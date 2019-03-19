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
#include <boost/program_options.hpp>

const size_t BATCH_SIZE = 8192; // reads to input, analyse
const size_t QUEUE_MAX = 16;        // number of batches in input queue
const size_t QUEUE_HIGH_LEVEL = 12; // max number of batch before running input
const size_t OUTPUT_QUEUE_SUM = 65535; // max read in ouput queues
const std::streamsize GZIP_INPUT_BUFFER_SIZE = 16*1024;
const std::streamsize GZIP_OUTPUT_BUFFER_SIZE = 4*1024;

using namespace std;
using namespace boost::iostreams;
namespace po = boost::program_options;

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


class FastqInput {
    unique_ptr<istream> input_stream;

    public:
    string path;
    size_t input_sizes[4] = {1};

    FastqInput(const string& path, bool compressed) : path(path) {
        if (compressed) {
            filtering_istream* fis = new filtering_istream;
            input_stream.reset(fis);
            file_descriptor_source raw(path);
            // Construct with buffer size. 15 is the default value of the first parameter.
            fis->push(gzip_decompressor(15, GZIP_INPUT_BUFFER_SIZE));
            fis->push(raw);
        }
        else {
            input_stream.reset(new ifstream(path, ios::binary));
        }
    }

    bool readBatch(Batch* bat) {
        int i;
        for (i=0; i<BATCH_SIZE; ++i) {
            for (int j=0; j<4; ++j) bat->data[i][j].reserve(input_sizes[j]);
            getline(*input_stream, bat->data[i][0]);
            input_sizes[0] = max(input_sizes[0], bat->data[i][0].capacity());
            if (input_stream->good()) {
                for (int j=1; j<4; ++j) {
                    getline(*input_stream, bat->data[i][j]);
                    input_sizes[j] = max(input_sizes[j], bat->data[i][j].capacity());
                }
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


class FastqOutput {
    shared_ptr<ostream> out_stream;

    public:
    string path;

    bool openFile(string basepath, bool compressed) {
        if (compressed) {
            path = basepath + ".fastq.gz";
            filtering_ostream* fos = new filtering_ostream;
            out_stream.reset(fos);
            file_descriptor_sink fds(path);
            fos->push(gzip_compressor(zlib::default_compression, GZIP_OUTPUT_BUFFER_SIZE));
            fos->push(fds);
            return (fds.is_open());
        }
        else {
            path = basepath + ".fastq";
            out_stream.reset(new ofstream(path, ios::binary));
            return out_stream->good();
        }
    }

    void writeBuffers(vector<OutputJob>& outputs) {
        for(OutputJob & output : outputs) {
            (*out_stream) << output.data[0] << '\n';
            (*out_stream) << output.data[1].substr(output.trim) << '\n';
            (*out_stream) << output.data[2] << '\n';
            (*out_stream) << output.data[3].substr(output.trim) << '\n';
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
    bool use_levens, no_trim;
    unsigned int barcode_mismatches, alignment_mismatches;

    public:

    unsigned long undetermined_reads = 0;
    unsigned long n_total_reads = 0;


    Analysis(vector<Sample>& samples, std::array<unsigned int, 2> bclen, bool use_levens,
            unsigned int barcode_mismatches, unsigned int alignment_mismatches, bool no_trim) :
        samples(samples), bclen(bclen), use_levens(use_levens), no_trim(no_trim),
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
                            if (no_trim) {
                                results[idata].trim[0] = 0;
                                results[idata].trim[1] = 0;
                            }
                            else {
                                results[idata].trim[0] = n_trim_r[0];
                                results[idata].trim[1] = n_trim_r[1];
                            }
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
            if (++n_total_reads % 1000000 == 0) {
                cerr << "Processed " << n_total_reads << " reads. Undetermined: "
                     << (undetermined_reads * 100.0 / n_total_reads) << " %."<< endl;
            }
        }
        return results;
    }
};


class DemultiplexingManager {

    const size_t n_thread;
 
    std::array<FastqInput, 2>& inputs;
    mutex inputmx[2];
    queue<Batch*> inputqs[2];
    atomic_bool input_busy[2];
   
    atomic_bool analysis_busy;
    mutex analysismx;
    Analysis & analysis;
    
    std::array<vector<FastqOutput>, 2> outputs;
    const size_t noutqs, ngsamples; // ngsamples=number of generalised samples, includes undetermined
    vector<vector<OutputJob>> outputqs;
    vector<mutex> outputqmx, outputmx;

    mutex any_work_mutex;
    condition_variable cv_any_work;

    atomic_bool error;
    atomic_size_t n_outputs_queued;

    public:
    DemultiplexingManager(size_t n_thread,
                        std::array<FastqInput, 2>& inputs,
                        Analysis& analysis,
                        std::array<vector<FastqOutput>, 2>& outputs) 
        : n_thread(n_thread), inputs(inputs), analysis(analysis),
          outputs(outputs), noutqs(outputs[0].size()+outputs[1].size()),
          ngsamples(outputs[0].size()), outputqs(noutqs), outputqmx(noutqs),
          outputmx(noutqs)
     {
         input_busy[0] = false;
         input_busy[1] = false;
         analysis_busy = false;
         error = false;
         n_outputs_queued = 0;
    }

    // Spawns threads and runs the demultiplexing. The current thread is used as one of
    // the workers, so n_threads-1 threads are created. The function returns when
    // the analysis is complete.
    bool execute() {
        int i;
        vector <thread> workers;
        for (i=0; i<n_thread; ++i) {
            workers.emplace_back(thread(&DemultiplexingManager::run, this, i));
        }
        //run(0);
        for (thread& t : workers) t.join();
        return !error;
    }

    // Main thread function is a loop that continues until all data have been processed.
    // It first tries to do an input task, then analysis, then any output task.
    // Whenever a task is completed, it starts from the top again.
    void run(int thread_index) {
        bool finish_local = false;
        while (!error && !finish_local) {
            bool inputed = false;
            // Try to start input processing
            for (int i=0; i<2 && !inputed; ++i) {
                // Select one of the two input queues, different threads
                // will try different queues first.
                int qindex = (i ^ thread_index) & 1;
                size_t qsize = inputqs[qindex].size();
                if (qsize < QUEUE_HIGH_LEVEL && 
                        !inputs[qindex].eof() &&
                        !input_busy[qindex]) {
                    inputed = runInputFunction(qindex);
                }
            }
            if (inputed) continue;
            if (!analysis_busy && n_outputs_queued < OUTPUT_QUEUE_SUM) {
                if (runAnalysisFunction()) continue;
            }

            bool outputed = false;
            for (int i=0; i<noutqs && !outputed; ++i) {
                // Use different starting point for each thread
                int qindex = (i + thread_index*(noutqs / n_thread)) % noutqs;
                if (!outputqs[qindex].empty()) {
                    outputed = runOutputFunction(qindex);
                }
            }
            if (outputed) continue;

            // We havent 'continue'd in the above code: it means there was nothing to do!
            // Either we are done, or there is not enough parallel tasks available at this point.
 
            // Check if this is really the end: Is the analysis complete?
            // Set finish_local, or if any check fails, we will go another round.
            if (inputs[0].eof() && inputs[1].eof()) {
                bool producer_finished = false;
                {
                    unique_lock<mutex> lk0(inputmx[0]);
                    unique_lock<mutex> lk1(inputmx[1]);
                    unique_lock<mutex> lk2(analysismx);
                    bool any_input_busy = false, empty = true;
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
                    if (all_empty) {
                        finish_local = true;
                        notifyNewWork(true);
                    }
                }
            }
            if (!error) { // There appears to be more tasks to do. Wait until the signal or timeout.
                unique_lock<mutex> lk(any_work_mutex);
                cv_any_work.wait_for(lk, std::chrono::seconds(5));
            }
        }
    }

    // Calls input function. Returns true if it was executed, false if
    // there was already an active input job for this queue.
    bool runInputFunction(int qindex) {
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
                cerr << "Error while reading file " << inputs[qindex].path << endl;
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
                        getq(read, i).emplace_back(
                                    results[result_index].trim[read], 
                                    move(b[read]->data[result_index])
                                    );
                    }
                }
            }
        }
        n_outputs_queued.fetch_add((size_t)(results.size()*2));
    }

    bool runOutputFunction(int qindex) {
        unique_lock<mutex> lk(outputmx[qindex], try_to_lock);
        if (!lk) return false;
        vector<OutputJob> jobs;
        {
            unique_lock<mutex> lk2(outputqmx[qindex]);
            if (outputqs[qindex].empty()) return false;
            jobs.reserve(outputqs[qindex].capacity());
            swap(jobs, outputqs[qindex]);
        }
        n_outputs_queued.fetch_sub((size_t)jobs.size());
        notifyNewWork(false); // Notify one possible analysis job available
                              // (in case of output queue limit)
        getoutput(qindex).writeBuffers(jobs);
        return true;
    }

    void notifyNewWork(bool all) {
        unique_lock<mutex> lk(any_work_mutex);
        if (all) cv_any_work.notify_all();
        else     cv_any_work.notify_one();
    }

    // Functions to map between
    //  - one-dimensional array of output buffers
    //  - two dimensional representation with read 0/1 and sample index
    vector<OutputJob>& getq(int read, int sample) {
        return outputqs[sample + read * ngsamples];
    }

    mutex& getqmutex(int read, int sample) {
        return outputqmx[sample + read * ngsamples];
    }

    FastqOutput& getoutput(int qindex) {
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

    const string usage =
                "usage:\n " + string(argv[0]) + " [options] \\\n"
                "     BARCODE_FILE SAMPLE_SHEET \\\n" + 
                "     INPUT_FILE_R1 INPUT_FILE_R2 \\\n" + 
                "     OUTPUT_PREFIX \n\n";


    // Parse command line into these variables
    string    barcode_file, sample_sheet,
              input_file_r1, input_file_r2,
              output_prefix;
    unsigned int barcode_mismatches, num_threads;
    int alignment_mismatches;
    bool use_hamming, no_trim;

    unsigned int cores = thread::hardware_concurrency();

    cerr << "\ndemultiplexer " << VERSION << "\n" << endl;

    po::options_description visible("Allowed options");
    visible.add_options()
        ("barcode-mismatches,b", po::value<unsigned int>(&barcode_mismatches)->default_value(1),
            "Allowed mismatches in barcode.")
        ("alignment-mismatches,a", po::value<int>(&alignment_mismatches)->default_value(-1),
            "Allowed mismatches in alignment (default=barcode-mismatches+1).")
        ("use-hamming,H", po::bool_switch(&use_hamming),
            "Use Hamming distance instead of Levenshtein distance.")
        ("no-trim,n", po::bool_switch(&no_trim),
            "Disable trimming of spacers and barcodes.")
        ("threads,t", po::value<unsigned int>(&num_threads)->default_value(min(cores, 16u)),
            "Number of threads to use.")
        ("help,h", "Show this help message.")
    ;
    po::options_description positionals("Positional options(hidden)");
    positionals.add_options()
        ("BARCODE_FILE", po::value<string>(&barcode_file)->required(),
            "Barcode definition file")
        ("SAMPLE_SHEET", po::value<string>(&sample_sheet)->required(),
            "Sample sheet")
        ("INPUT_FILE_R1", po::value<string>(&input_file_r1)->required(),
            "Input file read 1")
        ("INPUT_FILE_R2", po::value<string>(&input_file_r2)->required(),
            "Input file read 2")
        ("OUTPUT_PREFIX", po::value<string>(&output_prefix)->required(),
            "Output prefix (for directory, use a trailing slash).")
    ;
    po::options_description all_options("Allowed options");
    all_options.add(visible);
    all_options.add(positionals);

    po::positional_options_description pos_desc;
    pos_desc.add("BARCODE_FILE", 1);
    pos_desc.add("SAMPLE_SHEET", 1);;
    pos_desc.add("INPUT_FILE_R1", 1);
    pos_desc.add("INPUT_FILE_R2", 1);
    pos_desc.add("OUTPUT_PREFIX", 1);

    po::variables_map vm;
    try {
        po::store(
                po::command_line_parser(argc, argv).options(all_options).positional(pos_desc).run(),
                vm
                );
        if (vm.count("help") > 0) {
            cerr << usage << visible << endl;
            return 0;
        }
        else {
            po::notify(vm);
        }
    }
    catch(po::error& e) 
    { 
      cerr << "ERROR: " << e.what() << "\n\n";
      cerr << usage << visible << endl; 
      return 1; 
    }

    bool use_levens = !use_hamming;
    if (alignment_mismatches == -1) alignment_mismatches = barcode_mismatches + 1;

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

    bool compressed =
            input_file_r1.size() > 3 &&
            input_file_r1.substr(input_file_r1.size() - 3) == ".gz" &&
            input_file_r2.size() > 3 &&
            input_file_r2.substr(input_file_r1.size() - 3) == ".gz";

    // Input files
    std::array<FastqInput, 2> inputs = {
        FastqInput(input_file_r1, compressed),
        FastqInput(input_file_r2, compressed)
    };

    size_t n_samples = samples.size();
    
    // Output files -- represent undetermined as the last output stream
    std::array<vector<FastqOutput>, 2> sample_outputs = {
        vector<FastqOutput>(n_samples+1),
        vector<FastqOutput>(n_samples+1)
    };

    // Open output files for all samples
    for (int i=0; i<n_samples; ++i) {
        for (int read=0; read<2; ++read) {
            string basepath = output_prefix + samples[i].name + "_R" + to_string(read+1);
            if (!sample_outputs[read][i].openFile(basepath, compressed)) {
                cerr << "Failed to open output file " << sample_outputs[read][i].path
                     << " for sample " << samples[i].name << endl;
                exit(1);
            }
        }
    }
    // and Undetermined
    if (! (sample_outputs[0][n_samples].openFile(output_prefix + "Undetermined_R1", compressed) && 
            sample_outputs[1][n_samples].openFile(output_prefix + "Undetermined_R2", compressed))) {
        cerr << "Failed to open output files for undetermined reads." << endl;
        exit(1);
    }

    // Print information on startup
    cerr.precision(1);
    cerr << fixed;
    cerr << "\nDemultiplexing " << samples.size() << " samples...\n\n";
    cerr << " Allowed barcode mismatches: " << barcode_mismatches << '\n';
    cerr << " String distance:            ";
    if (use_levens) cerr << "Levenshtein" << '\n';
    else            cerr << "Hamming" << '\n';
    cerr << " Alignment string distance:  " << alignment_mismatches << '\n';
    cerr << " Input/output compression:   " << (compressed ? "gzip" : "off") << '\n';
    cerr << " Threads:                    " << num_threads << '\n';
    cerr << endl;

    Analysis analysis(samples, bclen, use_levens, barcode_mismatches, alignment_mismatches, no_trim);
    DemultiplexingManager manager(num_threads, inputs, analysis, sample_outputs);
    bool success = manager.execute();

    if (!success) { // The Manager will print an error message
        cerr << "\nThe program exited due to an error." << endl;
        return 1;
    }
    else {
        for (int i=0; i<samples.size(); ++i) { // Deletes empty files, would be corrupted gzip files.
            bool delete_files = samples[i].n_reads == 0;
            if (delete_files) {
                unlink(sample_outputs[0][i].path.c_str());
                unlink(sample_outputs[1][i].path.c_str());
            }
        }
        cerr << "\nCompleted demultiplexing " << analysis.n_total_reads << " PE reads.\n" << endl;
        cout << "SAMPLE_NAME\tR1_BC\tR2_BC\tNUM_READS\tPCT_READS\tPCT_PERFECT_BARCODE\tPCT_SPACER_FAIL\n";
        for (Sample& sample : samples) {
            cout.precision(2);
            cout << fixed;
            cout << sample.name << '\t'
                << sample.barcode[0].id << '\t' << sample.barcode[1].id << '\t'
                << sample.n_reads << '\t'
                << sample.n_reads * 100.0 / max(analysis.n_total_reads, 1ul) << '\t'
                << sample.n_perfect_barcode * 100.0 / max(sample.n_reads, 1ul) << '\t'
                << sample.n_spacer_fail * 50.0 / max(sample.n_reads, 1ul)
                << '\n';
        }
        cout << "Undetermined\t-\t-\t" << analysis.undetermined_reads << '\t'
            << analysis.undetermined_reads * 100.0 / analysis.n_total_reads << "\t0\t0\n";
    }
    return 0;
}


