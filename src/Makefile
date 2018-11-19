all: demultiplexer

CFLAGS += -O3 -std=c++11

# On RHEL7 -- Works inside devtoolset-7
demultiplexer.static: demultiplexer.cpp bounded_levenshtein_distance.cpp
	$(CXX) -o $@ demultiplexer.cpp $(CFLAGS) -static -pthread -lboost_iostreams$(BOOST_LIB_SUFF) -lz

demultiplexer: demultiplexer.cpp bounded_levenshtein_distance.cpp
	$(CXX) -o $@ demultiplexer.cpp $(CFLAGS) -pthread -lboost_iostreams$(BOOST_LIB_SUFF) -lz

clean:
	rm -f demultiplexer
