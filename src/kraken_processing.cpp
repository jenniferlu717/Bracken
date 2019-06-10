/*********************************************************************
 * kraken_processing.cpp is used as part of the kmer2distr script
 * Copyright (C) 2016-2018 Jennifer Lu, jlu26@jhmi.edu
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
 * Updated: 2018/09/06
 */
#include "kraken_processing.h"

/*METHOD: Evaluate the kraken database file*/
void evaluate_kfile(string k_file, string o_file, const taxonomy *my_taxonomy, const map<int, taxonomy *> *taxid2node, const map<string, int> seqid2taxid, const int kmer_len, const int read_len){
    /*Read the file and get maps:
     * map of number seqid to the kmer distribution  
     * map of number seqid to taxid 
     */
    map<int, string> id2kmers;
    map<int, int> id2taxid;
    map<int, string> id2seqid;
    map<int, string> id2tandl;
    int num_reads = read_kfile(k_file, &seqid2taxid, &id2seqid, &id2kmers, &id2taxid, &id2tandl);
    /*For each seqid, in parallel, convert kmer distribution to read distribution*/
    convert_distribution(o_file, num_reads, &id2seqid, &id2kmers, &id2taxid, &id2tandl, my_taxonomy, taxid2node, kmer_len, read_len);
}
/***************************************************************************************/
/*METHOD: READ THE FILE AND SAVE THE SEQID TO KMERS AND TAXID*/
int read_kfile(string k_file, const map<string, int> *seqid2taxid, map<int, string> *id2seqid, map<int, string> *id2kmers, map<int, int> *id2taxid, map<int, string> *id2tandl)
{
    /*Initialize Variables*/
    int s_count = 0;
    int curr_taxid; 
    int pos1, pos2, pos3, pos4, pos5;
    string line; 
    string curr_seqid;
    string curr_kmers;
    string curr_ataxid, curr_len; 
    string ataxid_length;
    /*Read through the file*/ 
    ifstream krakenfile (k_file);
    if (krakenfile.is_open()){
        printf("\t>>STEP 3: READING DATABASE.KRAKEN FILE\n");
        printf("\t\t0 sequences read...");
        while(getline(krakenfile, line)) {
            s_count += 1;
            if (s_count % 10 == 0) 
                printf("\r\t\t%i sequences read...", s_count);
            //Find delimiter indices
            pos1 = line.find("\t");
            pos2 = line.find("\t", pos1+1);
            pos3 = line.find("\t", pos2+1);
            pos4 = line.find("\t", pos3+1);
            pos5 = line.find("\n", pos4+1);
            //Extract seqid and taxid 
            curr_seqid = line.substr(pos1+1, pos2-pos1-1);
            curr_ataxid = line.substr(pos2+1, pos3-pos2-1);
            curr_len = line.substr(pos3+1, pos4-pos3-1);
            curr_kmers = line.substr(pos4+1, pos5-pos4-1);
            curr_taxid = seqid2taxid->find(curr_seqid)->second;
            ataxid_length = curr_ataxid + "\t" + curr_len; 
            //Save in maps
            id2seqid->insert(std::pair<int, string>(s_count, curr_seqid)); 
            id2kmers->insert(std::pair<int, string>(s_count,curr_kmers)); 
            id2taxid->insert(std::pair<int, int>(s_count, curr_taxid));
            id2tandl->insert(std::pair<int, string>(s_count, ataxid_length));
        }
        printf("\r\t\t%i total sequences read\n", s_count);
    } else {
        printf("  cannot open %s", k_file.c_str());  
        exit(1);
    }
    return s_count; 
}
/***************************************************************************************/
/*METHOD: CONVERT DISTRIBUTIONS INTO READ MAPPINGS - SEND TO PRINT*/
void convert_distribution(string o_file, int s_count, const map<int, string> *id2seqid, const map<int, string> *id2kmers, const map<int, int> *id2taxid, const map<int, string> *id2tandl, const taxonomy *my_taxonomy, const map<int, taxonomy *> *taxid2node, const int kmer_len, const int read_len){
    /*Initialize variables for getting read mappings instead of kmer mappings */
    int n_kmers = read_len - kmer_len + 1;
    int seqs_read = 0;
    /*Iterate over taxid2kmers in parallel*/
    printf("\t>>STEP 4: CONVERTING KMER MAPPINGS INTO READ CLASSIFICATIONS:\n");
    printf("\t\t%imers, with a database built using %imers\n",read_len, kmer_len);
    cerr << "\t\t" << seqs_read << " sequences converted...";
    int i;
    #pragma omp parallel for
    for(i = 1; i <= s_count; i++) {
        #pragma omp critical
        {
            cerr << "\r\t\t" << seqs_read << " sequences converted...(up next: ";
            cerr << id2seqid->find(i)->second << ")";
        }
        //Get values to parse here
        string curr_ks = id2kmers->find(i)->second;
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
        map<int,int> taxids_mapped; 
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
                } else {
                    t_it->second += 1;
                }
                prev_taxid = mapped_taxid;
                //Remove last element
                prev_kmer = curr_kmers[0];
                curr_kmers.erase(curr_kmers.begin());
            } 
        }
        //Update User
        #pragma omp atomic
        seqs_read += 1;
        #pragma omp critical
        {
            //cerr << "\r\t\t" << seqs_read << " sequences converted...";
            print_distribution(o_file, id2seqid->find(i)->second, id2taxid->find(i)->second, id2tandl->find(i)->second, taxids_mapped);
        }
    }
    cerr << "\r\t\t" << seqs_read << " sequences converted...\n";
}


/***************************************************************************************/
/*METHOD: Process the array of the kmer distributions*/ 
int get_classification(vector<int> *curr_kmers, const taxonomy *my_taxonomy, const map<int, taxonomy *> *taxid2node) {
    //Initialize values
    map<int, int> taxid2kmers; 
    int count_ts = 0;
    int curr_t; 
    int save_t; 
    //Figure out number of taxids in current list 
    for (int i = 0; i < curr_kmers->size(); i++) {
        //Get the first kmer out 
        curr_t = curr_kmers->at(i);  
        //Check for unclassified 
        if (curr_t == 0){
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
        return 0;
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
                    n1 = n1->get_parent();
                }
                while (n1->get_lvl_num() < n2->get_lvl_num()) {
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
/***************************************************************************************/
/*METHOD - PRINT DISTRIBUTION*/ 
void print_distribution(string o_file, const string seqid, const int taxid, string ataxid_length, map<int, int> distribs){
    //Outstream 
    ofstream outfile;
    //Open file to append
    outfile.open(o_file, ofstream::app);
    //Print read information
    outfile << seqid << "\t"; 
    outfile << taxid << "\t"; 
    outfile << ataxid_length << "\t"; 
    //Print distributions
    for (map<int, int>::iterator it=distribs.begin(); it!=distribs.end(); ++it){
        outfile << it->first << ":" << it->second << " ";
    }
    outfile << "\n";
}

