all: demultiplexer

CFLAGS += -O3 -std=c++11
VERSION = $(shell git describe --tags --dirty)


# On RHEL7 -- Works inside devtoolset-7
demultiplexer.static: demultiplexer.cpp bounded_levenshtein_distance.cpp
	$(CXX) -o $@ demultiplexer.cpp $(CFLAGS) -static -pthread -lboost_program_options$(BOOST_LIB_SUFF) -lboost_iostreams$(BOOST_LIB_SUFF) -lz -DVERSION=\"${VERSION}\"

demultiplexer.debug: demultiplexer.cpp bounded_levenshtein_distance.cpp
	$(CXX) -o $@ demultiplexer.cpp -g -std=c++11 -pthread -lboost_program_options$(BOOST_LIB_SUFF) -lboost_iostreams$(BOOST_LIB_SUFF) -lz -DVERSION=\"${VERSION}\"

demultiplexer: demultiplexer.cpp bounded_levenshtein_distance.cpp
	$(CXX) -o $@ demultiplexer.cpp $(CFLAGS) -pthread -lboost_program_options$(BOOST_LIB_SUFF) -lboost_iostreams$(BOOST_LIB_SUFF) -lz -DVERSION=\"${VERSION}\"

clean:
	rm -f demultiplexer
