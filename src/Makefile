C++ = g++
CFLAGS = -c -g -O3 -pedantic -std=c++11 -fopenmp 

all: kmer_counter

kmer_counter: 	kmer2read_distr.o time.o taxonomy.o kraken_processing.o
		$(C++) -o kmer2read_distr time.o kmer2read_distr.o taxonomy.o kraken_processing.o -lgomp

clean:
		rm -f *.o

%.o:	%.cpp
		$(C++) $(CFLAGS) $*.cpp
