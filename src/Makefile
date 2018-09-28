all: demultiplexer

CFLAGS += -O3 -std=c++11

demultiplexer: demultiplexer.cpp
	$(CXX) -o $@ $^ $(CFLAGS) -pthread -lz -lboost_iostreams$(BOOST_LIB_SUFF)

clean:
	rm -f demultiplexer
