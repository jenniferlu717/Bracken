#ifndef KRAKEN_PROCESSING_H
#define KRAKEN_PROCESSING_H

#include "kmer2read_headers.h"
#include "taxonomy.h"
#include "time.h"
void evaluate_kfile(string, string, const taxonomy *, const map<int, taxonomy *> *, map<string, int>, const int, const int);

int read_kfile(string, const map<string, int> *, map<int, string> *, map<int, string> *, map<int, int> *, map<int, string> *);

void convert_distribution(string, int, const map<int, string> *, const map<int, string> *, const map<int, int> *, const map<int, string> *, const taxonomy *, const map<int, taxonomy *> *, const int, const int);

int get_classification(vector<int> *, const taxonomy *, const map<int, taxonomy *> *);

//Each line = outputfile
//id2seqid value
//id2taxid value
//id2tandl value
//read distributions
void print_distribution(string, const string, const int, string, map<int, int> );  
#endif
