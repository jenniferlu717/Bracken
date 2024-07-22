CXXFLAGS := -c -g -O3 -pedantic -std=c++11

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CXX ?= g++
else
	CXX ?= clang++
endif

IS_CLANG := $(findstring clang++,$(CXX))
ifeq ($(IS_CLANG),clang++)
	CXXFLAGS += -Xpreprocessor -fopenmp
	LDFLAGS += -lomp
else
	CXXFLAGS += -fopenmp
	LDFLAGS += -lgomp
endif

all: kmer2read_distr

kmer2read_distr: kmer2read_distr.o ctime.o taxonomy.o kraken_processing.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $<

