#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/program_options.hpp>
#include "thread_source.hpp"
#include "bounded_levenshtein_distance.cpp"

#define BUFFER_SIZE 1024
using namespace std;
using namespace boost::iostreams;

class BarcodeAndSpacer {
    public:
    string id, barcode, spacer;
    size_t total_length;
    BarcodeAndSpacer(string id, string barcode, string spacer)
        : id(id), barcode(barcode), spacer(spacer) {
        total_length = barcode.length() + spacer.length();
    }
};

class Sample {
public:
    BarcodeAndSpacer barcode[2];
    string name;
    unsigned long n_reads = 0, n_perfect_barcode = 0;
    filtering_ostream *out_r1 = nullptr, *out_r2 = nullptr;

    Sample(const string name, BarcodeAndSpacer bc1, BarcodeAndSpacer bc2) :
        name(name), barcode({bc1, bc2}) {}

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

vector<Sample> getSamples(const vector<BarcodeAndSpacer>& barcodes, const string& sample_sheet) {
    ifstream ssheet(sample_sheet);
    vector<Sample> samples;
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
    return samples;
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
        cerr << "usage: " << argv[0]
             << " BARCODE_FILE SAMPLE_SHEET INPUT_R1 INPUT_R2 OUTPUT_PREFIX [MISMATCHES_PER_READ=L1]" << endl;
        return 1;
    }

    const string    barcode_file(argv[1]), sample_sheet(argv[2]),
                    input_file_r1(argv[3]), input_file_r2(argv[4]),
                    output_prefix(argv[5]);

    unsigned int allow_mismatches = 1;
    bool use_levens = true;
    if (argc == 7) {
        if (argv[6][0] == 'L' || argv[6][0] == 'l') {
            use_levens = true;
            allow_mismatches = stoul(argv[6]+1);
        }
        else {
            use_levens = false;
            allow_mismatches = stoul(argv[6]);
        }
    }
    
    vector<BarcodeAndSpacer> barcodes = getBarcodesAndSpacers(barcode_file);
    if (barcodes.empty()) {
        cerr << "Error: Barcode file is empty." << endl;
        return 1;
    }
    vector<Sample> samples = getSamples(barcodes, sample_sheet);
    if (samples.empty()) {
        cerr << "Error: No samples found." << endl;
        return 1;
    }
    unsigned int bc1len = samples[0].barcode[0].barcode.length();
    unsigned int bc2len = samples[0].barcode[1].barcode.length();

    // Check that all samples have sample length for barcode 1 and barcode 2
    if (!all_of(samples.begin(), samples.end(),
                [bc1len,bc2len](const Sample& s) {
                    return s.barcode[0].barcode.length() == bc1len
                        && s.barcode[1].barcode.length() == bc2len;
                })) {
        cerr << "Error: Barcodes have different length." << endl;
        return 1;
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
    //stream<thread_source> input_r1(input_r1_ifs);
    //stream<thread_source> input_r2(input_r2_ifs);
    //input_r1->start();
    //input_r2->start();
    
    // Open output files
    for (Sample& sample : samples) {
        if (!sample.openFiles(output_prefix)) {
            cerr << "Failed to open output file for sample " << sample.name << endl;
            exit(1);
        }
    }
    filtering_ostream und_r1, und_r2;
    und_r1.push(gzip_compressor());
    und_r1.push(file_descriptor_sink(output_prefix + "Undetermined_R1.fq.gz"));
    und_r2.push(gzip_compressor());
    und_r2.push(file_descriptor_sink(output_prefix + "Undetermined_R2.fq.gz"));

    cerr << "\nDemultiplexing " << samples.size() << " samples...\n\n";
    cerr << " Allowed mismatches: " << allow_mismatches << '\n';
    cerr << " String distance:    ";
    if (use_levens) cerr << "Levenshtein" << '\n' << endl;
    else            cerr << "Hamming" << '\n' << endl;

    // Buffer for 1 FASTQ record (4 lines) from each file R1 and R2
    string data[2][4];
    unsigned long undetermined_reads = 0;
    int i;
    unsigned long n_total_reads = 0;
    while (input_r1 && input_r2) {
        // One record from each FQ
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
        
        // List of pairs (mismatches, index) of samples matching barcode1. index is an
        // index into the main list of samples.
        list<pair<unsigned int,int>> bc1_matching_sample;
        bool undetermined = true;
        if (data[0][1].length() >= bc1len && data[1][1].length() >= bc2len) {
            for (Sample& sample : samples) {
                unsigned int bc_mismatches[2];
                if (use_levens) {
                    bc_mismatches[0] = bounded_levenshtein_distance(allow_mismatches + 1,
                            bc1len, sample.barcode[0].barcode, bc1len, data[0][1]);
                }
                else {
                    bc_mismatches[0] = mismatches(sample.barcode[0].barcode, data[0][1]);
                }
                if (bc_mismatches[0] <= allow_mismatches) {
                    if (use_levens) {
                        bc_mismatches[1] = bounded_levenshtein_distance(allow_mismatches + 1,
                                bc2len, sample.barcode[1].barcode, bc2len, data[1][1]);
                    }
                    else {
                        bc_mismatches[1] = mismatches(sample.barcode[1].barcode, data[1][1]);
                    }

                    if (bc_mismatches[1] <= allow_mismatches) {
                        size_t n_trim_r[2];
                        // Determine how much to trim for each of R1, R2
                        for (int i=0; i<2; ++i) {
                            if (use_levens) {
                                if (bc_mismatches[i] < 0) {
                                    // Align spacer only
                                    n_trim_r[i] = aligned_s2_length(
                                            sample.barcode[i].spacer.length() + 1,
                                            sample.barcode[i].spacer,
                                            data[i][1].substr(
                                                    sample.barcode[i].barcode.length()
                                                    )
                                            ) + sample.barcode[i].barcode.length();
                                }
                                else {
                                    // Align barcode + spacer sequence
                                    string bcsp = sample.barcode[i].barcode +
                                                        sample.barcode[i].spacer;
                                    n_trim_r[i] = aligned_s2_length(
                                            sample.barcode[i].spacer.length() +
                                                allow_mismatches + 1,
                                            bcsp,
                                            data[i][1]
                                            );
                                }
                            }
                            else {
                                n_trim_r[i] = sample.barcode[i].total_length;
                            }
                        }

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
                << sample.n_reads * 100.0 / n_total_reads << '\t'
                << sample.n_perfect_barcode * 100.0 / sample.n_reads << '\n';
        }
        cout << "------------------------------------------------------\n";
        cout << "Undetermined\t" << undetermined_reads << '\t'
            << undetermined_reads * 100.0 / n_total_reads << "\t-\n";
    }
    return 0;
}


