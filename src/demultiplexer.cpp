#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <queue>
#include <deque>
#include <memory>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <chrono>
#include <sstream>
#include <fstream>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

const size_t BATCH_SIZE = 1024; // reads

using namespace std;
using namespace boost::iostreams;
#include "bounded_levenshtein_distance.cpp"

// General data classes for barcode and sample
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


// Infrastructure for thread pool and work queues

class Batch {
    public:
    size_t num_in_batch;
    std::array<string, 4> data[BATCH_SIZE];
};


class CompressedFastqInput {
    unique_ptr<filtering_istream> input_stream;

    public:
    string path;

    CompressedFastqInput(const string& path) : path(path) {
        input_stream.reset(new filtering_istream);
        file_descriptor_source raw(path);
        input_stream->push(gzip_decompressor());
        input_stream->push(raw);
    }

    bool readBatch(Batch* bat) {
        int i;
        for (i=0; i<BATCH_SIZE; ++i) {
            getline(*input_stream, bat->data[i][0]);
            if (input_stream) {
                for (int j=1; j<4; ++j)
                    getline(*input_stream, bat->data[i][j]);
            }
            else if (input_stream->eof()) {
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
        return input_stream->eof();
    }
};


class OutputJob {
    public:
    int trim;
    std::array<string, 4> data;
    OutputJob(int trim, std::array<string, 4>&& data) :
        trim(trim), data(data) {}
};


class CompressedFastqOutput {
    shared_ptr<filtering_ostream> out_stream;
    string out_path;

    public:
    bool openFile(string path) {
        out_stream.reset(new filtering_ostream);
        file_descriptor_sink fds(path);
        out_stream->push(gzip_compressor());
        out_stream->push(fds);
        if (fds.is_open()) {
            out_path = path;
            return true;
        }
        else {
            return false;
        }
    }

    void writeBuffers(queue<OutputJob>& outputs) {
        while (!outputs.empty()) {
           OutputJob& output = outputs.front();
           (*out_stream) << output.data[0] << '\n';
           (*out_stream) << output.data[1].substr(output.trim) << '\n';
           (*out_stream) << output.data[2] << '\n';
           (*out_stream) << output.data[3].substr(output.trim) << '\n';
           outputs.pop();
        }
    }
};


class AnalysisResult {
    // Contains an int, routing, defining an output sample. For undetermined, a value of
    // the number of samples is used (i.e. one higher than the highest sample index).
    // Also contains the number of characters to trim off each read.
    public:
    int routing;
    int trim[2];
};


int mismatches(const string& templ, const string& test) {
    unsigned int mismatch = 0;
    size_t len = templ.length();
    for (int i = 0; i < len; ++i) {
        if (templ[i] != test[i]) mismatch++;
    }
    return mismatch;
}

class Analysis {

    vector<Sample>& samples;

    vector<BarcodeAndSpacer> barcode0_list;
    vector<list<int> > barcode0_sample_list_list;
    std::array<unsigned int, 2> bclen;
    bool use_levens;
    unsigned int barcode_mismatches, alignment_mismatches;

    public:

    unsigned long undetermined_reads = 0;
    unsigned long n_total_reads = 0;


    Analysis(vector<Sample>& samples, std::array<unsigned int, 2> bclen, bool use_levens,
            unsigned int barcode_mismatches, unsigned int alignment_mismatches) :
        samples(samples), bclen(bclen), use_levens(use_levens),
        barcode_mismatches(barcode_mismatches), alignment_mismatches(alignment_mismatches)
    {
        // Make a list to map barcode1 to a list of samples with that barcode1
        for (int sample_index=0; sample_index<samples.size(); ++sample_index) {
            Sample& s = samples[sample_index];
            bool found = false;
            for (int i=0; i<barcode0_list.size(); ++i) {
                if (s.barcode[0].barcode == barcode0_list[i].barcode) {
                    found = true;
                    barcode0_sample_list_list[i].push_back(sample_index);
                    break;
                }
            }
            if (!found) {
                barcode0_list.push_back(s.barcode[0]);
                barcode0_sample_list_list.push_back(list<int>(1, sample_index));
            }
        }
    }

    vector<AnalysisResult> analyseAndRouteReads(std::array<Batch*, 2>& bat) {

        vector<AnalysisResult> results(bat[0]->num_in_batch);

        cout << "debug running analysis of a batch with " << bat[0]->num_in_batch << " reads" << endl;
        for (int idata=0; idata < bat[0]->num_in_batch; ++idata) {
            bool undetermined = true;
            if (bat[0]->data[idata][1].length() >= bclen[0] && bat[1]->data[idata][1].length() >= bclen[1]) {
                int bc_mismatches[2], n_trim_r[2];
                int bc0_match_index = -1;

                // Try to match barcode 1
                for (int i = 0; bc0_match_index == -1 && i < barcode0_list.size(); ++i) {
                    BarcodeAndSpacer & bcs = barcode0_list[i];
                    if (use_levens) {
                        auto result = mismatch_and_alignment(barcode_mismatches + 1,
                                alignment_mismatches + 1,
                                bclen[0], bcs.full, bat[0]->data[idata][1]);
                        bc_mismatches[0] = result.first;
                        n_trim_r[0] = result.second;
                    }
                    else {
                        bc_mismatches[0] = mismatches(bcs.barcode, bat[0]->data[idata][1]);
                        n_trim_r[0] = bcs.full.length();
                    }
                    if (bc_mismatches[0] <= barcode_mismatches) bc0_match_index = i;
                }

                if (bc0_match_index != -1) {
                    // Loop over samples to match bc2
                    for (int sample_index : barcode0_sample_list_list[bc0_match_index]) {
                        Sample & sample = samples[sample_index];
                        BarcodeAndSpacer & bcs = sample.barcode[1];
                        if (use_levens) {
                            auto result = mismatch_and_alignment(barcode_mismatches + 1,
                                    alignment_mismatches + 1, bclen[1], bcs.full,
                                    bat[1]->data[idata][1]);
                            bc_mismatches[1] = result.first;
                            n_trim_r[1] = result.second;
                        }
                        else {
                            bc_mismatches[1] = mismatches(bcs.barcode, bat[1]->data[idata][1]);
                            n_trim_r[1] = bcs.full.length();
                        }

                        if (bc_mismatches[1] <= barcode_mismatches) {
                            for (int i=0; i<2; ++i) {
                                if (n_trim_r[i] == 0) {
                                    sample.n_spacer_fail++;
                                    n_trim_r[i] = sample.barcode[i].full.length();
                                }
                            }
                            // Sample barcode match -- output trimmed FQ
                            sample.n_reads++;
                            if (bc_mismatches[0] + bc_mismatches[1] == 0) sample.n_perfect_barcode++;
                            results[idata].routing = sample_index;
                            results[idata].trim[0] = n_trim_r[0];
                            results[idata].trim[1] = n_trim_r[1];
                            undetermined = false;
                            break;
                        }
                    }
                }
            }
            // No match found, output untrimmed FQ
            if (undetermined) {
                undetermined_reads++;
                results[idata].routing = samples.size();
                results[idata].trim[0] = 0;
                results[idata].trim[1] = 0;
            }
            if (++n_total_reads % 10000 == 0) {
                cerr << "Processed " << n_total_reads << " reads. Undetermined: "
                     << (undetermined_reads * 100.0 / n_total_reads) << " %."<< endl;
            }
        }
        return results;
    }
};


class DemultiplexingManager {

    const size_t QUEUE_MAX = 16;
    const size_t QUEUE_HIGH_LEVEL = 12;
    const size_t OUTPUT_QUEUE_MAX = 4096;

    const size_t n_thread;
 
    std::array<CompressedFastqInput, 2>& inputs;
    mutex inputmx[2];
    queue<Batch*> inputqs[2];
    atomic_bool input_busy[2];
   
    atomic_bool analysis_busy;
    mutex analysismx;
    Analysis & analysis;
    
    std::array<vector<CompressedFastqOutput>, 2> outputs;
    const size_t noutqs, ngsamples; // ngsamples=number of generalised samples, includes undetermined
    vector<queue<OutputJob>> outputqs;
    vector<mutex> outputqmx, outputmx;

    mutex any_work_mutex;
    condition_variable cv_any_work;

    atomic_bool error;

    public:
    DemultiplexingManager(size_t n_thread,
                        std::array<CompressedFastqInput, 2>& inputs,
                        Analysis& analysis,
                        std::array<vector<CompressedFastqOutput>, 2>& outputs) 
        : n_thread(n_thread), inputs(inputs), analysis(analysis),
          outputs(outputs), noutqs(outputs[0].size()+outputs[1].size()),
          ngsamples(outputs[0].size()), outputqs(noutqs), outputqmx(noutqs),
          outputmx(noutqs)
     {
         input_busy[0] = false;
         input_busy[1] = false;
         analysis_busy = false;
         error = false;
    }

    // Spawns threads and runs the demultiplexing. The current thread is used as one of
    // the workers, so n_threads-1 threads are created. The function returns when
    // the analysis is complete.
    bool execute() {
        int i;
        //vector <thread> workers;
        //for (i=0; i<n_thread-1; ++i) {
        //    workers.emplace_back(thread(run, i));
        //}
        run(0);
    }

    // Main thread function is a loop that continues until all data have been processed.
    // It first tries to do an input task, then analysis, then any output task.
    // Whenever a task is completed, it starts from the top again.
    void run(int thread_index) {
        bool finish_local = false;
        while (!error && !finish_local) {
            // Try to start input processing
            for (int i=0; i<2; ++i) {
                // Select one of the two input queues, different threads
                // will try different queues first.
                int qindex = (i ^ thread_index) & 1;
                size_t qsize = inputqs[qindex].size();
                if (qsize < QUEUE_HIGH_LEVEL && 
                        !inputs[qindex].eof() &&
                        !input_busy[qindex]) {
                    if (runInputFunction(qindex)) continue;
                }
            }
            size_t output_max = 0;
            for (auto& q : outputqs) output_max = max(q.size(), output_max);
            if (!analysis_busy && output_max < OUTPUT_QUEUE_MAX) {
                if (runAnalysisFunction()) continue;
            }

            cout << "debug no analysis function to run now, going to output!" << endl;

            for (int i=0; i<noutqs; ++i) {
                // Use different starting point for each thread
                int qindex = (i + thread_index*(noutqs / n_thread)) % noutqs;
                if (!outputqs[qindex].empty()) {
                    if (runOutputFunction(qindex)) continue;
                }
            }

            // We havent 'continue'd in the above code: it means there was nothing to do!
            // Either we are done, or there is not enough parallel tasks available at this point.
            cout << "debug Hit the end of task selection loop, nothing to do now?" << endl;
 
            // Check if this is really the end: Is the analysis complete?
            // Set finish_local, or if any check fails, we will go another round.
            if (inputs[0].eof() && inputs[1].eof()) {
                bool producer_finished = false;
                {
                    unique_lock<mutex> lk0(inputmx[0]);
                    unique_lock<mutex> lk1(inputmx[1]);
                    unique_lock<mutex> lk2(analysismx);
                    bool any_input_busy = false, empty = false;
                    for (int i=0; i<2; ++i) {
                        empty = empty && inputqs[i].empty(); 
                        any_input_busy = any_input_busy || input_busy[i];
                    }
                    producer_finished = !analysis_busy && empty && !any_input_busy;
                    if (!analysis_busy && !any_input_busy && (inputqs[0].empty() ^ inputqs[1].empty())) { 
                        cerr << "Error: Input files have different number of reads!" << endl;
                        error = true;
                    }
                }
                if (producer_finished) {
                    // There won't be any new output tasks. So check that all output
                    // queues are empty.
                    bool all_empty = true;
                    for (int i=0; i<noutqs; ++i) {
                        unique_lock<mutex> lk(outputqmx[i]);
                        all_empty = all_empty && outputqs[i].empty();
                    }
                    if (all_empty) finish_local = true;
                }
            }
            if (!error) { // There appears to be more tasks to do. Wait until the signal or timeout.
                unique_lock<mutex> lk(any_work_mutex);
                cv_any_work.wait_for(lk, std::chrono::seconds(1));
            }
        }
    }

    // Calls input function. Returns true if it was executed, false if
    // there was already an active input job for this queue.
    bool runInputFunction(int qindex) {
        cout << "debug running input function" << endl;
        size_t queue_fill;
        {
            unique_lock<mutex> lk(inputmx[qindex]);
            if (input_busy[qindex]) return false;
            if (inputs[qindex].eof()) return false;
            input_busy[qindex] = true;
            queue_fill = inputqs[qindex].size();
        }
        while (queue_fill < QUEUE_MAX && !inputs[qindex].eof()) {
            Batch* bat = new Batch;
            if (inputs[qindex].readBatch(bat)) {
                unique_lock<mutex> lk(inputmx[qindex]);
                inputqs[qindex].push(bat);
                queue_fill = inputqs[qindex].size();
            }
            else {
                cerr << "Error while reading file " << 
                    inputs[qindex].path << endl;
                error = true;
                return true;
            }
        }
        input_busy[qindex] = false;
        notifyNewWork(false); // Notify one possible analysis job available
    }

    // Calls analysis function and feeds results into output queues
    bool runAnalysisFunction() {
        unique_lock<mutex> lk(analysismx, try_to_lock);
        if (!lk) return false;
        if (analysis_busy) return false;
        analysis_busy = true;
        cout << "debug running analysis " << endl;
        while (!inputqs[0].empty() && !inputqs[1].empty()) {
            std::array<Batch*, 2> b;
            for (int i=0; i<2; ++i) {
                unique_lock<mutex> lk2(inputmx[i]);
                b[i] = inputqs[i].front();
                inputqs[i].pop();
            }
            assert(b[0]->num_in_batch == b[1]->num_in_batch);
            vector<AnalysisResult> routings = analysis.analyseAndRouteReads(b);
            addResultsToQueues(b, routings);
            delete b[0];
            delete b[1];
        }
        cout << "debug completed analysis function" << endl;
        analysis_busy = false;
        notifyNewWork(true); // Notify all -- there may be many output jobs, and input jobs
    }

    // Route the elements in this batch to output queues
    void addResultsToQueues(std::array<Batch*, 2>& b, vector<AnalysisResult>& results) {
        // Collect all elements to add to each pair of output queues
        vector<int> addable[ngsamples];
        for (int j=0; j<results.size(); ++j) {
            addable[results[j].routing].push_back(j);
        }
        for (int i=0; i<(ngsamples); ++i) {
            // Process all results destined for queue # i
            // Now lock the queues and add the elements
            if (!addable[i].empty()) {
                for (int read=0; read<2; ++read) {
                    unique_lock<mutex> lk(getqmutex(read, i));
                    for (int result_index : addable[i]) {
                        // Move the data buffer for this read (0/1) into the queue for i
                        getq(read, i).push(OutputJob(
                                    results[result_index].trim[read], 
                                    move(b[read]->data[result_index])
                                    ));
                    }
                }
            }
        }
    }

    bool runOutputFunction(int qindex) {
        unique_lock<mutex> lk(outputmx[qindex], try_to_lock);
        if (!lk) return false;
        queue<OutputJob> jobs;
        {
            unique_lock<mutex> lk2(outputqmx[qindex]);
            if (outputqs[qindex].empty()) return false;
            swap(jobs, outputqs[qindex]);
        }
        notifyNewWork(false); // Notify one possible analysis job available
                              // (in case of output queue limit)
        getoutput(qindex).writeBuffers(jobs);
    }

    void notifyNewWork(bool all) {
        unique_lock<mutex> lk(any_work_mutex);
        if (all) cv_any_work.notify_all();
        else     cv_any_work.notify_one();
    }

    // Functions to map between
    //  - one-dimensional array of output buffers
    //  - two dimensional representation with read 0/1 and sample index
    queue<OutputJob>& getq(int read, int sample) {
        return outputqs[sample + read * ngsamples];
    }

    mutex& getqmutex(int read, int sample) {
        return outputqmx[sample + read * ngsamples];
    }

    CompressedFastqOutput& getoutput(int qindex) {
        return outputs[qindex / ngsamples][qindex % ngsamples];
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
    std::array<unsigned int, 2> bclen = {
        (unsigned int)samples[0].barcode[0].barcode.length(),
        (unsigned int)samples[0].barcode[1].barcode.length()
    };

    // Check that all samples have sample length for barcode 1 and barcode 2
    if (!all_of(samples.begin(), samples.end(),
                [bclen](const Sample& s) {
                    return s.barcode[0].barcode.length() == bclen[0]
                        && s.barcode[1].barcode.length() == bclen[1];
                })) {
        cerr << "Error: Barcodes have different length." << endl;
        return 1;
    }

    // Input files
    std::array<CompressedFastqInput, 2> inputs = {
        CompressedFastqInput(input_file_r1),
        CompressedFastqInput(input_file_r2)
    };

    size_t n_samples = samples.size();
    
    // Output files -- represent undetermined as the last output stream
    std::array<vector<CompressedFastqOutput>, 2> sample_outputs = {
        vector<CompressedFastqOutput>(n_samples+1),
        vector<CompressedFastqOutput>(n_samples+1)
    };

    // Open output files for all samples
    for (int i=0; i<n_samples; ++i) {
        for (int read=0; read<2; ++read) {
            string path = output_prefix + samples[i].name + "_R" + to_string(read+1) + ".fq.gz";
            if (!sample_outputs[read][i].openFile(path)) {
                cerr << "Failed to open output file " << path << " for sample "
                     << samples[i].name << endl;
                exit(1);
            }
        }
    }
    // and Undetermined
    if (! (sample_outputs[0][n_samples].openFile(output_prefix + "Undetermined_R1.fq.gz") && 
            sample_outputs[1][n_samples].openFile(output_prefix + "Undetermined_R2.fq.gz"))) {
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

    Analysis analysis(samples, bclen, use_levens, barcode_mismatches, alignment_mismatches);
    DemultiplexingManager manager(4, inputs, analysis, sample_outputs);
    bool success = manager.execute();

    if (!success) { // The Manager will print an error message
        cerr << "\nThe program exited due to an error." << endl;
        return 1;
    }
    else {
        for (Sample sample : samples) { // Deletes empty files, would be corrupted gzip files.
            bool delete_files = sample.n_reads == 0;
            if (delete_files) {
                unlink(sample.path1.c_str());
                unlink(sample.path2.c_str());
            }
        }
        cerr << "\nCompleted demultiplexing " << analysis.n_total_reads << " PE reads.\n" << endl;
        cout << "SAMPLE_NAME\tNUM_READS\tPCT_READS\tPCT_PERFECT_BARCODE\tPCT_SPACER_FAIL\n";
        cout << "---------------------------------------------------------------\n";
        for (Sample& sample : samples) {
            cout.precision(2);
            cout << sample.name << '\t' << sample.n_reads << '\t'
                << fixed
                << sample.n_reads * 100.0 / max(analysis.n_total_reads, 1ul) << '\t'
                << sample.n_perfect_barcode * 100.0 / max(sample.n_reads, 1ul) << '\t'
                << sample.n_spacer_fail * 50.0 / max(sample.n_reads, 1ul)
                << '\n';
        }
        cout << "---------------------------------------------------------------\n";
        cout << "Undetermined\t" << analysis.undetermined_reads << '\t'
            << analysis.undetermined_reads * 100.0 / analysis.n_total_reads << "\t-\n";
    }
    return 0;
}


