#include <iostream>
#include <algorithm>
#include <vector>
#include <mutex>
#include <atomic>
#include <queue>
#include <deque>
#include <sstream>
#include <fstream>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

const size_t BATCH_SIZE = 100; // reads

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
    unsigned long n_reads = 0, n_perfect_barcode = 0, n_spacer_fail = 0;
    string path1, path2;

    Sample(const string name, BarcodeAndSpacer bc0, BarcodeAndSpacer bc1) :
        name(name), barcode({bc0, bc1}) {}

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

class Batch {
    size_t num_in_batch;
    size_t trim[BATCH_SIZE][2];
    string data[BATCH_SIZE][4];
};

class DemultiplexingManager {

    const size_t QUEUE_MAX = 16;
    const size_t QUEUE_HIGH_LEVEL = 12;
    //const size_t QUEUE_LOW_LEVEL = 4;
    //
    CompressedFastqInput & inputs[2];
    mutex inputmx[2];
    queue<Batch> inputqs[2];
    atomic_bool input_busy[2] = {false, false};
   
    atomic_bool analysis_busy = false;
    mutex analysismx;
    Analysis & analysis;
    
    vector<queue<Batch>> outputqs;
    vector<mutex> outputmx;

    bool finish = false, error = false;

    public:
    DemultiplexingManager(CompressedFastqInputs& inputs[2]) 
        : inputs(inputs) {
    }

    void run(int thread_index) {
        while (!error && !finish) {
            // Try to start input processing
            for (int i=0; i<1; ++i) {
                // Select one of the two input queues, different threads
                // will try different queues first.
                int qindex = (i ^ thread_index) & 1;
                size_t qsize = inputqs[qindex].size();
                if (qsize < QUEUE_HIGH_LEVEL && 
                        !input[qindex].eof() &&
                        !input_busy[qindex]) {
                    if (runInputFunction(qindex)) continue;
                }
            }
            if (!analysis_busy) {
                if (runAnalysisFunction()) continue;
            }
            
        }
    }

    // Calls input function. Returns true if it was executed, false if
    // there was already an active input job for this queue.
    bool runInputFunction(int qindex) {
        size_t queue_fill;
        {
            unique_lock lk(inputmx[qindex]);
            if (input_busy[qindex]) return false;
            if (inputs[qindex].eof) return false;
            input_busy[qindex] = true;
            queue_fill = inputqs[qindex].size();
        }
        while (queue_fill < QUEUE_MAX && !inputs[qindex].eof()) {
            Batch* bat = new Batch;
            if (inputs[qindex].readBatch(bat)) {
                {
                    unique_lock lk(inputmx[qindex]);
                    inputqs[qindex].push(bat);
                    queue_fill = inputqs[qindex].size();
                }
            }
            else {
                cerr << "Error while reading file " << 
                    input[qindex].path << endl;
                error = true;
                return true;
            }
        }
        input_busy[qindex] = false;
    }

    // Calls analysis function and feeds results into output queues
    bool runAnalysisFunction() {
        unique_lock lk(analysismx, try_lock);
        if (!lk) return false;
        if (analysis_busy) return false;
        analysis_busy = true;
        while (!inputqs[0].empty() && !inputqs[1].empty()) {
            Batch* b[2];
            for (int i=0; i<2; ++i) {
                unique_lock lk2(inputmx[i]);
                b[i] = inputqs[i].front;
                inputqs[i].pop();
            }
            assert(b[0]->num_in_batch == b[1]->num_in_batch);
            vector<int> routings = analysis.analyseAndRoute(b);
            vector<int> addable;
            for (int i=0; i<outputqs.size(); ++i) {
                for (int j=0; j<routings.size(); ++j) {
                    if (routings[j] == i) {
                        addable.push_back(j);
                    }
                }
            }
            
        }
        analysis_busy = false;
    }
};

class CompressedFastqInput {
    filtering_istream input_stream;

    public:
    string path;

    CompressedFastqInput(const string& path) : path(path) {
        file_descriptor_source raw(path);
        input_stream.push(gzip_decompressor());
        input_stream.push(raw);
    }

    bool readBatch(Batch* bat) {
        int i;
        for (i=0; i<BATCH_SIZE; ++i) {
            getline(input_stream, bat->data[i][0]);
            if (input_stream) {
                for (j=1; j<4; ++j)
                    getline(input_stream, bat->data[i][j]);
            }
            else if (input_stream.eof()) {
                break;
            }
            else {
                return false;
            }
        }
        bat->num_in_batch = i;
        return true;
    }
    
    bool eof() {
        return input_stream.eof();
    }
}


class Analysis {

    const vector<Sample>& samples;

    vector<BarcodeAndSpacer> barcode0_list;
    vector<list<Sample*> > barcode0_sample_list_list;

    public:

    unsigned long undetermined_reads = 0;
    unsigned long n_total_reads = 0;


    Analysis(const vector<Sample>& samples) : samples(samples) {
        // Make a list to map barcode1 to a list of samples with that barcode1
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
    }

    vector<int> analyseAndRouteReads(Batch* bat) {
        vector<int> result(bat->num_in_batch, -1);

    }
}


class CompressedFastqOutput {
    filtering_ostream *out_stream = nullptr;
    string out_path;

    ~CompressedFastqOutput() {
        delete out_r1;
    }
    bool openFile(string path) {
        file_descriptor_sink fds(path);
        out_stream = new filtering_ostream();
        out_stream->push(gzip_compressor());
        out_stream->push(fds1);
        if (fds.is_open()) {
            out_path = path;
            return true;
        }
        else {
            delete out_stream;
            out_stream = nullptr;
            return false;
        }
    }
    void writeBuffers() {
    }
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

    // Input files
    CompressedFastqInput input1(input_file_r1), input2(input_file_r2);
    
    // Output files
    vector<CompressedFastqOutput> sample_outputs_r1(samples.size());
    vector<CompressedFastqOutput> sample_outputs_r2(samples.size());
    CompressedFastqOutput undetermined_output1, undetermined_output2;

    // Open output files for all samples
    for (int i=0; i<samples.size(); ++i) {
        string path1 = prefix + samples[i].name + "_R1.fq.gz";
        string path2 = prefix + samples[i].name + "_R2.fq.gz";

        if (!sample_outputs_r1.openFile(path)) {
            cerr << "Failed to open output file " << path << " for sample "
                 << sample.name << endl;
            exit(1);
        }
        if (!sample_outputs_r2.openFile(path)) {
            cerr << "Failed to open output file " << path << " for sample "
                 << sample.name << endl;
            exit(1);
        }
    }
    // and Undetermined
    if (! (undetermined_output1.openFile(output_prefix + "Undetermined_R1.fq.gz") && 
            undetermined_output2.openFile(output_prefix + "Undetermined_R2.fq.gz"))) {
        cerr << "Failed to open undetermined files." << endl;
        exit(1);
    }

    // Print information on startup
    cerr.precision(2);
    cerr << "\nDemultiplexing " << samples.size() << " samples...\n\n";
    cerr << " Allowed barcode mismatches: " << barcode_mismatches << '\n';
    cerr << " String distance:            ";
    if (use_levens) cerr << "Levenshtein" << '\n';
    else            cerr << "Hamming" << '\n';
    cerr << " Alignment string distance:  " << alignment_mismatches;
    cerr << '\n' << endl;

    // Buffer for 1 FASTQ record (4 lines) from each file R1 and R2
    string data[2][4];
    int i;
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
                        for (int i=0; i<2; ++i) {
                            if (n_trim_r[i] == 0) {
                                sample.n_spacer_fail++;
                                n_trim_r[i] = sample.barcode[i].full.length();
                            }
                        }
                        
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
            cerr << "Processed " << n_total_reads << " reads. Undetermined: "
                 << (undetermined_reads * 100.0 / n_total_reads) << " %."<< endl;
        }
    }
    if (!input_r1.eof() || !input_r2.eof()) {
        cerr << "Error: Input error. Input valid flags: R1: " << (bool)input_r1
             << ", R2: " << (bool)input_r2 << endl;
        return 1;
    }
    else {
            bool delete_files = n_reads == 0 && out_r1 && out_r2;
            if (delete_files) {
                unlink(path1.c_str());
                unlink(path2.c_str());
            }
        cerr << "\nCompleted demultiplexing " << n_total_reads << " PE reads.\n" << endl;
        cout << "SAMPLE_NAME\tNUM_READS\tPCT_READS\tPCT_PERFECT_BARCODE\tPCT_SPACER_FAIL\n";
        cout << "---------------------------------------------------------------\n";
        for (Sample& sample : samples) {
            cout.precision(2);
            cout << sample.name << '\t' << sample.n_reads << '\t'
                << fixed
                << sample.n_reads * 100.0 / max(n_total_reads, 1ul) << '\t'
                << sample.n_perfect_barcode * 100.0 / max(sample.n_reads, 1ul) << '\t'
                << sample.n_spacer_fail * 50.0 / max(sample.n_reads, 1ul)
                << '\n';
        }
        cout << "---------------------------------------------------------------\n";
        cout << "Undetermined\t" << undetermined_reads << '\t'
            << undetermined_reads * 100.0 / n_total_reads << "\t-\n";
    }
    return 0;
}


