/*********************************************************************
 * kraken_processing.cpp is used as part of the kmer2distr script
 * Copyright (C) 2016-2020 Jennifer Lu, jlu26@jhmi.edu
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
 * Updated: 2020/03/31
 */
#include "kraken_processing.h"

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
                    char * lineEnd = &lineStart[pos]; 
                    string kraken_line(lineStart, lineEnd - lineStart);
                    //Variables for things to save
                    string seqid = "";
                    int taxid = -1;
                    map<int, int> taxids_mapped;
                    //CALL METHOD TO PROCESS THE LINE
                    convert_line(kraken_line, &seqid2taxid, read_len, kmer_len, my_taxonomy, taxid2node, seqid, taxid, taxids_mapped); 
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
                        for (map<int, int>::iterator it=taxids_mapped.begin(); it!=taxids_mapped.end(); ++it){
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
/***************************************************************************************/
/*METHOD: CONVERT DISTRIBUTIONS INTO READ MAPPINGS - SEND TO PRINT*/
void convert_line(string line, const map<string,int> *seqid2taxid, const int read_len, const int kmer_len, const taxonomy *my_taxonomy, const map<int, taxonomy *> *taxid2node, string &seqid, int &taxid, map<int,int> &taxids_mapped){ 
   
    //Find delimiter indices
    int pos1, pos2, pos3, pos4, pos5;
    pos1 = line.find("\t");
    pos2 = line.find("\t", pos1+1);
    pos3 = line.find("\t", pos2+1);
    pos4 = line.find("\t", pos3+1);
    pos5 = line.find("\n", pos4+1);
    //Extract seqid and taxid 
    seqid = line.substr(pos1+1, pos2-pos1-1);
    taxid = seqid2taxid->find(seqid)->second;
    string curr_ks = line.substr(pos4+1, pos5-pos4-1);

    /*Initialize variables for getting read mappings instead of kmer mappings */
    int n_kmers = read_len - kmer_len + 1;
    
    //Saving values 
    vector<int> all_kmers;
    int count_kmers = 0;
    //Iterate through all of the kmer pairs 
    int mid, end; 
    string buf; 
    std::stringstream ss(curr_ks);
    int pair_taxid, pair_count;
    string curr_pair, pair_tstr;
    while(ss >> buf) {
        //Extract the string of this pair
        curr_pair = buf;
        //Split up this pair into the taxid and the number of kmers
        mid = curr_pair.find(":");
        end = curr_pair.find("\n"); 
        pair_tstr = curr_pair.substr(0,mid);
        pair_count = atoi(curr_pair.substr(mid+1,end-mid-1).c_str());
        if (pair_tstr == "A")
            pair_taxid = 0;
        else
            pair_taxid = atoi(pair_tstr.c_str());
        //Add kmers to queue
        for (int j = 0; j < pair_count; j++) {
            all_kmers.push_back(pair_taxid);
            count_kmers += 1;
        }
    }
    //Process all mappings
    vector<int> curr_kmers; 
    int mapped_taxid;
    int prev_kmer, next_kmer; 
    int prev_taxid; 
    for (int k = 0; k < count_kmers; k++) {
        next_kmer = all_kmers[k];
        curr_kmers.push_back(next_kmer);
        if (curr_kmers.size() == n_kmers) {
            if (prev_kmer == next_kmer) {
                mapped_taxid = prev_taxid;
            } else {
                mapped_taxid = get_classification(&curr_kmers, my_taxonomy, taxid2node);
            } 
            //if (mapped_taxid != 0) {
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
            prev_kmer = curr_kmers[0];
            curr_kmers.erase(curr_kmers.begin());
        } 
    }
}


/***************************************************************************************/
/*METHOD: Process the array of the kmer distributions*/ 
int get_classification(vector<int> *curr_kmers, const taxonomy *my_taxonomy, const map<int, taxonomy *> *taxid2node) {
    //Initialize values
    map<int, int> taxid2kmers; 
    int count_ts = 0;
    int curr_t; 
    int save_t; 
    int found_root = 0;
    //Figure out number of taxids in current list 
    for (int i = 0; i < curr_kmers->size(); i++) {
        //Get the first kmer out 
        curr_t = curr_kmers->at(i);  
        //Check for unclassified 
        if (curr_t == 0){
            continue;
        } else if (curr_t == 1){
            found_root = 1;
            continue;   
        }
        //Save taxids
        auto x_it = taxid2kmers.find(curr_t);
        if (x_it == taxid2kmers.end()) {
            //Not found in map
            count_ts += 1;
            save_t = curr_t;
            taxid2kmers[curr_t] = 1;
        } else {
            x_it->second += 1;
        }
    }
    //Determine what to send back
    if (count_ts == 0) {
        if (found_root != 0){
            return 1;
        } else {
            return 0;
        }
    } else if (count_ts == 1) {
        return save_t; 
    } else {
        //MUST MAKE THIS FASTER SOMEHOWWW
        //can i sort by node numbers? while array has something -- pop if visit.....
        //Score each taxid -- make faster by ignoring taxids already visited? 
        map<int, int> taxid2scores; 
        taxonomy * curr_n;
        for (map<int, int>::iterator it=taxid2kmers.begin(); it!=taxid2kmers.end(); ++it){
            taxid2scores[it->first] = it->second;
            curr_n = taxid2node->find(it->first)->second;
            while(curr_n->get_parent() != NULL) {
                curr_n = curr_n->get_parent();
                if (taxid2kmers.find(curr_n->get_taxid()) != taxid2kmers.end()) {
                    taxid2scores[it->first] += taxid2kmers[curr_n->get_taxid()];
                }
            }
        }
        //Find the maximum score 
        int max_score = 0;
        int max_taxid = 0;
        for (map<int, int>::iterator it=taxid2scores.begin(); it!=taxid2scores.end(); ++it){
            if (it->second > max_score){
                //Get max scoring
                max_score = it->second;
                max_taxid = it->first;
            } else if (it->second == max_score){
                //Break ties 
                taxonomy * n1 = taxid2node->find(it->first)->second;
                taxonomy * n2 = taxid2node->find(max_taxid)->second;
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
}
