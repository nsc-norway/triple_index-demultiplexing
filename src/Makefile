all: demultiplexer

CFLAGS += -O3 -std=c++11

# On RHEL7 -- Works inside devtoolset-7
demultiplexer.static: demultiplexer.cpp bounded_levenshtein_distance.cpp
	$(CXX) -o $@ demultiplexer.cpp $(CFLAGS) -lz /usr/lib64/libboost_iostreams.a

demultiplexer: demultiplexer.cpp bounded_levenshtein_distance.cpp
	$(CXX) -o $@ demultiplexer.cpp $(CFLAGS) -lz -lboost_iostreams$(BOOST_LIB_SUFF)

clean:
	rm -f demultiplexer
