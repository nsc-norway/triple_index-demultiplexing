#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


using namespace std;
using namespace boost::iostreams;
#include "bounded_levenshtein_distance.cpp"

class BarcodeAndSpacer {
    public:
    string id, barcode, spacer, full;
    size_t total_length;
    BarcodeAndSpacer(string id, string barcode, string spacer)
        : id(id), barcode(barcode), spacer(spacer), full(barcode+spacer) {
    }
};

class Sample {
public:
    BarcodeAndSpacer barcode[2];
    string name;
    unsigned long n_reads = 0, n_perfect_barcode = 0;
    filtering_ostream *out_r1 = nullptr, *out_r2 = nullptr;

    Sample(const string name, BarcodeAndSpacer bc0, BarcodeAndSpacer bc1) :
        name(name), barcode({bc0, bc1}) {}

    ~Sample() {
        delete out_r1;
        delete out_r2;
    }

    bool openFiles(string prefix) {
        file_descriptor_sink fds1(prefix + name + "_R1.fq.gz"), fds2(prefix + name + "_R2.fq.gz");
        out_r1 = new filtering_ostream();
        out_r1->push(gzip_compressor());
        out_r1->push(fds1);
        out_r2 = new filtering_ostream();
        out_r2->push(gzip_compressor());
        out_r2->push(fds2);
        return fds1.is_open() && fds2.is_open();
    }

};

vector<BarcodeAndSpacer> getBarcodesAndSpacers(const string& barcode_file) {
    ifstream bc(barcode_file);
    if (!bc) {
        cerr << "Error: Unable to open barcode file " << barcode_file << endl;
        exit(1);
    }
    vector<BarcodeAndSpacer> result;
    string line;
    while (bc) {
        getline(bc, line);
        if (!line.empty()) {
            stringstream ss_line(line);
            string id, barcode, spacer;
            ss_line >> id; ss_line >> barcode; ss_line >> spacer;
            result.emplace_back(BarcodeAndSpacer(id, barcode, spacer));
        }
    }
    return result;
}

void getSamples(vector<Sample>& samples, const vector<BarcodeAndSpacer>& barcodes, const string& sample_sheet) {
    ifstream ssheet(sample_sheet);
    int linenum = 1;
    while (ssheet) {
        string line;
        getline(ssheet, line);
        if (!line.empty()) {
            stringstream linestream(line);
            string name, bc1id, bc2id;
            linestream >> name; linestream >> bc1id; linestream >> bc2id;
            if (!linestream && !linestream.eof()) {
                cerr << "Incomplete data on sample sheet line " << linenum << "." << endl;
                exit(1);
            }
            const BarcodeAndSpacer * bc1 = nullptr, * bc2 = nullptr;
            for (const BarcodeAndSpacer & bctest : barcodes) {
                if (bc1id == bctest.id) bc1 = &bctest;
                if (bc2id == bctest.id) bc2 = &bctest;
            }
            if (!bc1) {
                cerr << "Error: Barcode 1 ID " << bc1id << " not found, for sample " << name
                     << " line " << linenum << "." << endl;
                exit(1);
            }
            if (!bc2) {
                cerr << "Error: Barcode 2 ID " << bc2id << " not found, for sample " << name 
                     << " line " << linenum << "." << endl;
                exit(1);
            }
            samples.emplace_back(Sample(name, *bc1, *bc2));
        }
        linenum++;
    }
}

int mismatches(const string& templ, const string& test) {
    unsigned int mismatch = 0;
    size_t len = templ.length();
    for (int i = 0; i < len; ++i) {
        if (templ[i] != test[i]) mismatch++;
    }
    return mismatch;
}

int main(int argc, char* argv[]) {

    if (argc < 6) {
        cerr << "usage:\n " << argv[0]
             << " BARCODE_FILE SAMPLE_SHEET INPUT_R1 INPUT_R2 OUTPUT_PREFIX \\\n"
             << "             [BARCODE_MISMATCHES_PER_READ=L1 \\\n"
             << "             [ALIGNMENT_MISMATCHES=BARCODE_MISMATCHES_PER_READ]]" << endl;
        return 1;
    }

    // First parse options: most are strings, except the mismatch parameters
    const string    barcode_file(argv[1]), sample_sheet(argv[2]),
                    input_file_r1(argv[3]), input_file_r2(argv[4]),
                    output_prefix(argv[5]);

    // Barcode mismatches can be either Ln for Levenshtein distance n, or n
    // for Hamming distance n.
    unsigned int barcode_mismatches = 1, alignment_mismatches;
    bool use_levens = true;
    if (argc >= 7) {
        if (argv[6][0] == 'L' || argv[6][0] == 'l') {
            barcode_mismatches = stoul(argv[6]+1);
            use_levens = barcode_mismatches > 0;
        }
        else {
            barcode_mismatches = stoul(argv[6]);
            use_levens = false;
        }
    }
    // Alignment distance can optionally be set higher than barcode mismatch,
    // to allow the alignment of more errors in spacer sequences.
    if (argc >= 8) {
        alignment_mismatches = stoul(argv[7]);
    }
    else {
        alignment_mismatches = barcode_mismatches;
    }
    
    // Read named barcodes from a single tab-separated file. This file should also
    // contain the heterogeneity spacer sequences.
    vector<BarcodeAndSpacer> barcodes = getBarcodesAndSpacers(barcode_file);
    if (barcodes.empty()) {
        cerr << "Error: Barcode file is empty." << endl;
        return 1;
    }
    // Read samples from sample sheet. Each sample record refers to the names of the
    // barcode sequences.
    vector<Sample> samples;
    getSamples(samples, barcodes, sample_sheet);
    if (samples.empty()) {
        cerr << "Error: No samples found." << endl;
        return 1;
    }
    unsigned int bc0len = samples[0].barcode[0].barcode.length();
    unsigned int bc1len = samples[0].barcode[1].barcode.length();

    // Check that all samples have sample length for barcode 1 and barcode 2
    if (!all_of(samples.begin(), samples.end(),
                [bc0len,bc1len](const Sample& s) {
                    return s.barcode[0].barcode.length() == bc0len
                        && s.barcode[1].barcode.length() == bc1len;
                })) {
        cerr << "Error: Barcodes have different length." << endl;
        return 1;
    }

    // Make a list to map barcode1 to a list of samples with that barcode1
    vector<BarcodeAndSpacer> barcode0_list;
    vector<list<Sample*> > barcode0_sample_list_list;
    for (Sample& s : samples) {
        bool found = false;
        for (int i=0; i<barcode0_list.size(); ++i) {
            if (s.barcode[0].barcode == barcode0_list[i].barcode) {
                found = true;
                barcode0_sample_list_list[i].push_back(&s);
                break;
            }
        }
        if (!found) {
            barcode0_list.push_back(s.barcode[0]);
            barcode0_sample_list_list.push_back(list<Sample*>(1, &s));
        }
    }

    // Input files
    filtering_istream input_r1_ifs, input_r2_ifs;
    file_descriptor_source r1raw(input_file_r1), r2raw(input_file_r2);
    input_r1_ifs.push(gzip_decompressor());
    input_r1_ifs.push(r1raw);
    input_r2_ifs.push(gzip_decompressor());
    input_r2_ifs.push(r2raw);
    istream & input_r1 = input_r1_ifs;
    istream & input_r2 = input_r2_ifs;
    
    // Open output files for all samples
    for (Sample& sample : samples) {
        if (!sample.openFiles(output_prefix)) {
            cerr << "Failed to open output file for sample " << sample.name << endl;
            exit(1);
        }
    }
    // and Undetermined
    filtering_ostream und_r1, und_r2;
    und_r1.push(gzip_compressor());
    und_r1.push(file_descriptor_sink(output_prefix + "Undetermined_R1.fq.gz"));
    und_r2.push(gzip_compressor());
    und_r2.push(file_descriptor_sink(output_prefix + "Undetermined_R2.fq.gz"));

    // Print information on startup
    cerr << "\nDemultiplexing " << samples.size() << " samples...\n\n";
    cerr << " Allowed barcode mismatches: " << barcode_mismatches << '\n';
    cerr << " String distance:            ";
    if (use_levens) cerr << "Levenshtein" << '\n';
    else            cerr << "Hamming" << '\n';
    cerr << " Alignment string distance:  " << alignment_mismatches;
    cerr << '\n' << endl;

    // Buffer for 1 FASTQ record (4 lines) from each file R1 and R2
    string data[2][4];
    unsigned long undetermined_reads = 0;
    int i;
    unsigned long n_total_reads = 0;
    while (input_r1 && input_r2 /* && n_total_reads < 1 */) {

        // Read one record from each FQ
        getline(input_r1, data[0][0]);
        if (!input_r1) {
            // Make a dummy read on R2 if error / EOF on R1, to sync it up
            getline(input_r2, data[1][0]);
            break;
        }
        for (i=1; i<4; ++i)
            getline(input_r1, data[0][i]);
        for (i=0; i<4; ++i)
            getline(input_r2, data[1][i]);
        
        bool undetermined = true;
        size_t data0_len = data[0][1].length();
        size_t data1_len = data[1][1].length();
        if (data0_len >= bc0len && data1_len >= bc1len) {
            int bc_mismatches[2], n_trim_r[2];
            int bc0_match_index = -1;

            // Try to match barcode 1
            for (i = 0; bc0_match_index == -1 && i < barcode0_list.size(); ++i) {
                BarcodeAndSpacer & bcs = barcode0_list[i];
                if (use_levens) {
                    auto result = mismatch_and_alignment(barcode_mismatches + 1,
                            alignment_mismatches + 1,
                            bc0len, bcs.full, data[0][1]);
                    bc_mismatches[0] = result.first;
                    n_trim_r[0] = result.second;
                }
                else {
                    bc_mismatches[0] = mismatches(bcs.barcode, data[0][1]);
                    n_trim_r[0] = bcs.full.length();
                }
                if (bc_mismatches[0] <= barcode_mismatches) bc0_match_index = i;
            }

            if (bc0_match_index != -1) {
                // Loop over samples to match bc2
                for (Sample* pSample : barcode0_sample_list_list[bc0_match_index]) {
                    Sample & sample = *pSample;
                    BarcodeAndSpacer & bcs = sample.barcode[1];
                    if (use_levens) {
                        auto result = mismatch_and_alignment(barcode_mismatches + 1,
                                alignment_mismatches + 1, bc1len, bcs.full, data[1][1]);
                        bc_mismatches[1] = result.first;
                        n_trim_r[1] = result.second;
                    }
                    else {
                        bc_mismatches[1] = mismatches(bcs.barcode, data[1][1]);
                        n_trim_r[1] = bcs.full.length();
                    }

                    if (bc_mismatches[1] <= barcode_mismatches) {
                        // Sample barcode match -- output trimmed FQ to sample's output file
                        // (in case of Levenshtein distance, n_trim_r[i] may be zero if the
                        // alignment to the spacer seequence failed)
                        sample.n_reads++;
                        if (bc_mismatches[0] + bc_mismatches[1] == 0) sample.n_perfect_barcode++;

                        (*sample.out_r1) << data[0][0] << '\n';
                        (*sample.out_r1) << data[0][1].substr(n_trim_r[0]) << '\n';
                        (*sample.out_r1) << data[0][2] << '\n';
                        (*sample.out_r1) << data[0][3].substr(n_trim_r[0]) << '\n';

                        (*sample.out_r2) << data[1][0] << '\n';
                        (*sample.out_r2) << data[1][1].substr(n_trim_r[1]) << '\n';
                        (*sample.out_r2) << data[1][2] << '\n';
                        (*sample.out_r2) << data[1][3].substr(n_trim_r[1]) << '\n';
                        
                        undetermined = false;
                        break;
                    }
                }
            }
        }
        // No match found, output untrimmed FQ
        if (undetermined) {
            undetermined_reads++;
            for (i=0; i<4; ++i) und_r1 << data[0][i] << '\n';
            for (i=0; i<4; ++i) und_r2 << data[1][i] << '\n';
        }
        if (++n_total_reads % 1000000 == 0) {
            cerr << "Processed " << n_total_reads << " reads." << endl;
        }
    }
    if (!input_r1.eof() || !input_r2.eof()) {
        cerr << "Error: Input error. Input valid flags: R1: " << (bool)input_r1
             << ", R2: " << (bool)input_r2 << endl;
        return 1;
    }
    else {
        cerr << "\nCompleted demultiplexing " << n_total_reads << " PE reads.\n" << endl;
        cout << "SAMPLE_NAME\tNUM_READS\tPCT_READS\tPCT_PERFECT_BARCODE\n";
        cout << "------------------------------------------------------\n";
        for (Sample& sample : samples) {
            cout.precision(2);
            cout << sample.name << '\t' << sample.n_reads << '\t'
                << fixed
                << sample.n_reads * 100.0 / max(n_total_reads, 1ul) << '\t'
                << sample.n_perfect_barcode * 100.0 / max(sample.n_reads, 1ul)
                << '\n';
        }
        cout << "------------------------------------------------------\n";
        cout << "Undetermined\t" << undetermined_reads << '\t'
            << undetermined_reads * 100.0 / n_total_reads << "\t-\n";
    }
    return 0;
}


