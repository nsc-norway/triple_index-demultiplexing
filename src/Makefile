all: demultiplexer

CFLAGS += -O3 -std=c++11

demultiplexer: demultiplexer.cpp bounded_levenshtein_distance.cpp
	$(CXX) -o $@ demultiplexer.cpp $(CFLAGS) -pthread -lz -lboost_iostreams$(BOOST_LIB_SUFF)

clean:
	rm -f demultiplexer
