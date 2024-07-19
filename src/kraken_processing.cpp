/*********************************************************************
 * kraken_processing.cpp is used as part of the kmer2distr script
 * Copyright (C) 2016-2023 Jennifer Lu, jlu26@jhmi.edu
 *
 * This file is part of Bracken.
 * Bracken is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the license, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.*/
/************************************************************************
 * Jennifer Lu, jlu26@jhmi.edu
 * Updated: 2022/03/31
 */
#include <set>

#include "kraken_processing.h"

#include <limits>

struct TaxidInfo {
    int count;
    int score;
    std::set<int> active_ancestors;
    std::set<int> ancestors_seen;
    bool active_ancestors_updated;
};

// Second idea: Change the logic so that a kmer links to all the taxids that
// depend on it. That way instead of scanning the active kmers to find out who
// to update we will instead have a list of dependents that we know we have to
// update. Any dependent with a count of 0, we remove from the set because it
// has fallen off.
class KmerClassifier {
 public:
    KmerClassifier() {
        scores[0].count = 0;
        scores[1].count = 0;
    }

    int
    classify_kmers(int kmer_to_add, int kmer_to_remove,
                   const std::map<int, taxonomy *> *taxid2node) {
        if (kmer_to_remove >= 0) {
            scores[kmer_to_remove].count -= 1;
            scores[kmer_to_remove].score -= 1;

            if (scores[kmer_to_remove].count == 0) {
                active_kmers.erase(kmer_to_remove);
            }
        }

        bool new_kmer = false;

        if (kmer_to_add > 1) {
             if (active_kmers.find(kmer_to_add) == active_kmers.end()) {
                new_kmer = true;
                active_kmers.insert(kmer_to_add);
                scores[kmer_to_add].count = 0;
                scores[kmer_to_add].score = 0;
                scores[kmer_to_add].active_ancestors_updated = true;
            }

            if (new_kmer) {
                taxonomy *new_node = taxid2node->find(kmer_to_add)->second;
                auto comp = [] (taxonomy *a, taxonomy *b) {
                    return a->get_lvl_num() < b->get_lvl_num();
                };

                std::vector<taxonomy *> nodes;
                nodes.reserve(active_kmers.size() - 1);

                for (auto it = active_kmers.begin(); it != active_kmers.end(); it++) {
                    taxonomy *node = taxid2node->find(*it)->second;
                    if (node != new_node) {
                        nodes.push_back(node);
                    }
                }

                std::sort(nodes.begin(), nodes.end(), comp);

                while (nodes.size() != 0) {
                    auto node = nodes.back();
                    auto taxid = node->get_taxid();

                    if (node->get_lvl_num() < new_node->get_lvl_num()) {
                        break;
                    } else if (scores[taxid].ancestors_seen.find(kmer_to_add) != scores[taxid].ancestors_seen.end()) {
                        scores[taxid].score += 1;
                        scores[taxid].active_ancestors.insert(kmer_to_add);
                    } else {
                        while (node->get_parent() != NULL && node->get_lvl_num() > new_node->get_lvl_num()) {
                            node = node->get_parent();
                        }
                        if (node->get_taxid() == kmer_to_add) {
                            scores[taxid].score += 1;
                            scores[taxid].active_ancestors.insert(kmer_to_add);
                            scores[taxid].ancestors_seen.insert(kmer_to_add);
                        }
                    }
                    nodes.pop_back();
                }

                taxonomy *curr_n = new_node;
                auto score = scores[kmer_to_add];
                while (!nodes.empty()) {
                    auto node = nodes.back();
                    if (scores[kmer_to_add].ancestors_seen.find(node->get_taxid()) != scores[kmer_to_add].ancestors_seen.end()) {
                        scores[kmer_to_add].score += scores[node->get_taxid()].count;
                        scores[kmer_to_add].active_ancestors.insert(node->get_taxid());
                    } else {
                        curr_n = new_node;
                        while (curr_n->get_parent() != NULL /* && curr_n->get_lvl_num() < node->get_lvl_num()*/) {
                            curr_n = curr_n->get_parent();
                            if (curr_n->get_taxid() == nodes.back()->get_taxid()) {
                                scores[kmer_to_add].score += scores[curr_n->get_taxid()].count;
                                scores[kmer_to_add].active_ancestors.insert(curr_n->get_taxid());
                                scores[kmer_to_add].ancestors_seen.insert(curr_n->get_taxid());
                            }
                        }
                    }
                    nodes.pop_back();
                }
            }
        }

        scores[kmer_to_add].count += 1;
        scores[kmer_to_add].score += 1;

        if (active_kmers.size() == 0) {
            if (scores[1].count > 0) {
                return 1;
            } else {
                return 0;
            }
        }

        int max_score = 0;
        int max_taxid = 0;
        for (auto it = active_kmers.begin(); it != active_kmers.end(); it++) {
            if (!new_kmer || *it != kmer_to_add) {
                auto not_found = scores[*it].active_ancestors.end();
                auto item1 = scores[*it].active_ancestors.find(kmer_to_remove);
                auto item2 = new_kmer ? not_found : scores[*it].active_ancestors.find(kmer_to_add);

                if (item1 != not_found) {
                    scores[*it].score -= 1;
                    if (scores[kmer_to_remove].count == 0)
                        scores[*it].active_ancestors.erase(item1);
                }

                if (item2 != not_found) {
                    scores[*it].score += 1;
                }
            }

            if (scores[*it].score > max_score) {
                max_score = scores[*it].score;
                max_taxid = *it;
            } else if (scores[*it].score == max_score) {
                taxonomy *n1 = taxid2node->find(*it)->second;
                taxonomy *n2 = taxid2node->find(max_taxid)->second;
                //Get to the same level
                while (n1->get_lvl_num() > n2->get_lvl_num()) {
                    if (n1->get_parent() == NULL){
                        cerr << n1->get_taxid() << endl;
                    }

                    n1 = n1->get_parent();
                }
                while (n1->get_lvl_num() < n2->get_lvl_num()) {
                    if (n2->get_parent() == NULL){
                        cerr << n2->get_taxid() << endl;
                    }
                    n2 = n2->get_parent();
                }
                //Find LCA
                while (n1 != n2) {
                    n1 = n1->get_parent();
                    n2 = n2->get_parent();
                }
                max_taxid = n1->get_taxid();
            }
        }

        return max_taxid;
    }

    void reset() {
        scores[0].count = 0;
        scores[0].score = 0;
        scores[1].count = 0;
        scores[1].score = 0;

        for (auto it = active_kmers.begin(); it != active_kmers.end(); ++it) {
            scores[*it].count = 0;
            scores[*it].score = 0;
        }

        active_kmers.clear();
    }

 private:
 std::set<int> active_kmers;
 std::map<int, TaxidInfo> scores;
};

inline unsigned int fast_atou(const char *str)
{
    unsigned int val = 0;
    while(*str) {
        val = (val << 1) + (val << 3) + *(str++) - 48;
    }
    return val;
}

/*METHOD: Evaluate the kraken database file*/
void evaluate_kfile(string k_file, string o_file, const taxonomy *my_taxonomy, const map<int, taxonomy *> *taxid2node, const map<string, int> seqid2taxid, const int kmer_len, const int read_len){
    /*Parallel Variables*/

    FILE * kraken_file = fopen(k_file.c_str(),"r");
    int fd = fileno(kraken_file);
    struct stat sb;
    fstat(fd,&sb);
    size_t dataSize = sb.st_size;
    char * data = static_cast<char*>(mmap(NULL, dataSize, PROT_READ,MAP_PRIVATE,fd,0));
    size_t globalLineCounter = 0;

    int seqs_read = 0;
    /*Iterate over kraken file in parallel*/
    printf("\t>>STEP 3: CONVERTING KMER MAPPINGS INTO READ CLASSIFICATIONS:\n");
    printf("\t\t%imers, with a database built using %imers\n",read_len, kmer_len);
    cerr << "\t\t0 sequences converted...";
    //Open file to append
    ofstream outfile;
    outfile.open(o_file, ofstream::out);

    #pragma omp parallel
    {
        bool processMore = true;
        size_t lastReadPos = 0;
        size_t newLineCnt = 0;
        size_t i;
        string kraken_line;
        KmerClassifier classifier;

        //Get the line and process
        while(processMore) {
            //Get line number
            int currLineToProcess = __sync_fetch_and_add(&globalLineCounter,1);
            //Get line
            for (i = lastReadPos; i < dataSize; i++){
                newLineCnt += (data[i]=='\n');
                if(newLineCnt == currLineToProcess){
                    char * lineStart = &data[i + (data[i]=='\n' ? 1 : 0)];
                    size_t pos=0;
                    while(lineStart[pos] != '\n' && (i + pos) < dataSize){
                        pos++;
                    }
                    char *lineEnd = &lineStart[pos];
                    size_t len = lineEnd - lineStart;
                    kraken_line.resize(len);
                    kraken_line.assign(lineStart, len);
                    //Variables for things to save
                    string seqid = "";
                    int taxid = -1;
                    std::map<int, int> taxids_mapped;

                    //CALL METHOD TO PROCESS THE LINE
                    convert_line(kraken_line, &seqid2taxid, read_len, kmer_len, my_taxonomy, taxid2node, seqid, taxid, taxids_mapped, classifier);
                    //PRINT FOR LINE
                    #pragma omp atomic
                    seqs_read += 1;
                    #pragma omp critical
                    {
                        cerr << "\r\t\t" << seqs_read << " sequences converted (finished: ";
                        cerr << seqid << ")";
                        //Print read information
                        outfile << seqid << "\t";
                        outfile << taxid << "\t";
                        outfile << "" << "\t";
                        //Print distributions
                        for (auto it=taxids_mapped.begin(); it!=taxids_mapped.end(); ++it){
                            outfile << it->first << ":" << it->second << " ";
                        }
                        outfile << "\n";
                    }
                    lastReadPos = i+1;
                    break;
                }
            }
            if(i == dataSize) {
                processMore = false;
            }
        }
    }
    cerr << "\r\t\t" << globalLineCounter << " sequences converted\n";

}

// /***************************************************************************************/
// /*METHOD: CONVERT DISTRIBUTIONS INTO READ MAPPINGS - SEND TO PRINT*/
void convert_line(string line, const std::map<string,int> *seqid2taxid, const int read_len, const int kmer_len, const taxonomy *my_taxonomy, const std::map<int, taxonomy *> *taxid2node, string &seqid, int &taxid, std::map<int,int> &taxids_mapped, KmerClassifier &classifier){
    int pos1, pos2, pos3, pos4, pos5;
    pos1 = line.find("\t");
    pos2 = line.find("\t", pos1+1);
    pos3 = line.find("\t", pos2+1);
    pos4 = line.find("\t", pos3+1);
    pos5 = line.find("\n", pos4+1);
    //Extract seqid and taxid
    seqid = line.substr(pos1 + 1, pos2 - pos1 - 1);
    taxid = seqid2taxid->find(seqid)->second;
    char *curr_ks = &line[pos4 + 1];

    /*Initialize variables for getting read mappings instead of kmer mappings */
    int n_kmers = read_len - kmer_len + 1;

    //Saving values
    vector<std::pair<int, int>> all_kmers;
    size_t count_kmers = 0;
    //Iterate through all of the kmer pairs
    int mid, end;
    int pair_taxid, pair_count;
    string curr_pair, pair_tstr;
    size_t len = line.size() - pos4;
    size_t n_alloc = 0;
    curr_ks[len - 1] = ' ';
    for (int i = 0; i < len; i++) {
        if (curr_ks[i] == ' ')
          n_alloc++;
    }
    all_kmers.reserve(n_alloc);
    for (size_t i = 0; i < len; i++) {
        // Split up this pair into the taxid and the number of kmers
        char *data = &curr_ks[i];
        char *mid = (char *)memchr(data, (int)':', len - i);
        if (mid == NULL) break;
        char *end = (char *)memchr(data, (int)' ', len - i);
        *mid = '\0';
        *end  = '\0';
        pair_count = fast_atou(mid + 1);
        if (data[0] == 'A')
            pair_taxid = 0;
        else {
            pair_taxid = fast_atou(data);
        }
        // Add kmers to queue
        all_kmers.push_back(std::make_pair(pair_taxid, pair_count));
        count_kmers += 1;

        i = end - curr_ks;
    }

    //Process all mappings
    deque<int> curr_kmers;
    int mapped_taxid = -1;
    uint32_t prev_kmer = -1, next_kmer;
    uint32_t prev_taxid;
    uint32_t hash1, hash2;
    size_t hash;

    for (int k = 0; k < count_kmers; k++) {
        next_kmer = all_kmers[k].first;
        int count = all_kmers[k].second;
        for (int j = 0; j < count; j++) {
            curr_kmers.push_back(next_kmer);
            if (curr_kmers.size() == n_kmers) {
                if (prev_kmer == next_kmer) {
                    mapped_taxid = prev_taxid;
                } else {
                    mapped_taxid = classifier.classify_kmers(
                        next_kmer, prev_kmer, taxid2node);
                }
                //Save to map
                auto t_it = taxids_mapped.find(mapped_taxid);
                if (t_it == taxids_mapped.end()){
                    taxids_mapped[mapped_taxid] = 1;
                    //.insert(std::pair<int,int>(mapped_taxid, 1));
                } else {
                    t_it->second += 1;
                }
                prev_taxid = mapped_taxid;
                //Remove last element
                prev_kmer = curr_kmers.front();
                curr_kmers.pop_front();
            } else {
                classifier.classify_kmers(next_kmer, prev_kmer, taxid2node);
            }
        }
    }
    classifier.reset();
}
