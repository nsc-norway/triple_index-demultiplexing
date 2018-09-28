#ifndef THREAD_INPUT_FILTER_INCLUDED
#define THREAD_INPUT_FILTER_INCLUDED

#include <thread>
#include <condition_variable>
#include <mutex>
#include <memory>
#include <boost/iostreams/categories.hpp>

/**
 * Filter to buffer input and perform input read operations in
 * another thread. read() calls to the istream are made in a 
 * different thread, so all filters/devices upstream of this one
 * will run in that thread.
 *
 * This only works if the istream is always the same, which I think
 * it usually is. Why is istream an argument to read() and not to
 * the constructor?
 */

#include <chrono>
using namespace std;

class thread_source {

    condition_variable cv_output_ready, cv_input_ready;
    mutex m;
    istream& wrapped;
    streamsize buffer_size;
    thread read_worker_thread; 
    vector<char> buffer[2];
    // Buffer size is the size of the array, buffer_content is the
    // number of bytes of data in it.
    size_t buffer_content[2], read_ptr = 0;
    int output_buffer_index = 0;
    bool terminate = false, eof = false, error= false, started = false;
    bool buffer_output[2] = {false, false};

public:
    typedef char                         char_type;
    typedef boost::iostreams::source_tag category;

    thread_source(istream& wrapped, streamsize buffer_size = 1024*1024)
        : wrapped(wrapped), buffer_size(buffer_size) {
    }

    thread_source(const thread_source& other)
        : wrapped(other.wrapped), buffer_size(other.buffer_size) {}

    ~thread_source() {
        terminate = true;
        cv_input_ready.notify_all();
        //cv_output_ready.notify_all(); // We're not threadsafe that way
        if (read_worker_thread.joinable()) read_worker_thread.join();
    };

    void start() {
        if (!started) {
            for (int i=0; i<2; ++i) {
                buffer[i].resize(buffer_size);
            }
            read_worker_thread = thread(&thread_source::readWorkerLoop, this);
            started = true;
        }
    }

    // Provide data to the user of this source
    streamsize read(char* c, streamsize n) {
        size_t num_read = 0;
        while (num_read < n) {
            if (!buffer_output[output_buffer_index]) {
                unique_lock<mutex> lk(m);
                if (eof && !buffer_output[0] && !buffer_output[1]) {
                    if (num_read > 0) return num_read;
                    else return -1;
                }
                cv_input_ready.wait(lk,
                        [&]{return buffer_output[output_buffer_index] || error || eof;});
            }
            if (error) {
                throw new ios_base::failure("Wrapped stream failed");
            }
            streamsize ncpy = min(n - num_read,
                    buffer_content[output_buffer_index] - read_ptr);
            memcpy(c + num_read, buffer[output_buffer_index].data() + read_ptr, ncpy);
            read_ptr += ncpy;
            num_read += ncpy;

            if (read_ptr == buffer_content[output_buffer_index]) { // Swap buffers?
                buffer_output[output_buffer_index] = false;
                cv_output_ready.notify_all();
                output_buffer_index ^= 1;
                read_ptr = 0;
            }
        }
        return num_read;
    }
    
    // Thread function to read data from wrapped source
    void readWorkerLoop() {
        int mybuffer = 0;
        while (!terminate && !eof && !error) {
            // Input into this buffer
            wrapped.read(buffer[mybuffer].data(), buffer_size);

            if (!wrapped) {
                if (wrapped.eof()) {
                    buffer_content[mybuffer] = wrapped.gcount();
                    eof = true;
                }
                else {
                    buffer_content[mybuffer] = 0;
                    error = true;
                }
            }
            else {
                buffer_content[mybuffer] = buffer_size;
            }

            // Swap buffers
            {
                unique_lock<mutex> lk(m);

                // Clear input flag on buffer and allow output
                if (buffer_content[mybuffer] != 0)
                    buffer_output[mybuffer] = true;
                // Notify, of error / eof / new data
                cv_input_ready.notify_all();

                // Move on to next buffer
                mybuffer ^= 1;

                // Wait until output is finished on this buffer
                cv_output_ready.wait(lk,
                        [&]{return !buffer_output[mybuffer] || terminate || eof;}
                        );
            }
        }
    }

};

#endif // #ifndef THREAD_INPUT_FILTER_INCLUDED
