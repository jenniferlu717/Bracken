/*********************************************************************
 * kmer2read_distr.cpp is the main function for this package
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
 * Updated: 2018/09/06
 */

#include "kmer2read_headers.h"
#include "ctime.h"
#include "taxonomy.h"
#include "kraken_processing.h"

#define MAX_LEN_PRINT 50

/*General Function Declarations*/
void parse_command_line(int argc, char **argv);
void usage(int exit_code=0);

/*Function Declarations*/ 
void construct_taxonomy(const string, taxonomy *);
void get_seqid2taxid(string, map<string, int> *); 
/*Variables - Remains Constant*/
int num_threads = 1; 
int kmer_len = 31;
int read_len = 100; 
string taxid_file = "";
string seqid_file = "";
string kraken_file = ""; 
string output_file = "";
/*Other Program variables*/
map<string, int> seqid2taxid;
taxonomy *my_taxonomy = new taxonomy();
map<int, taxonomy *> taxid2node;
/*Main Driver Program*/
int main(int argc, char *argv[]) {
    /*set Default number of threads*/    
    omp_set_num_threads(1);
    
    /*Parse command line*/
    printf("\t>>STEP 0: PARSING COMMAND LINE ARGUMENTS\n");
    parse_command_line(argc, argv);  
    printf("\t\tTaxonomy nodes file: %s\n", taxid_file.c_str());
    printf("\t\tSeqid file:          %s\n", seqid_file.c_str());
    printf("\t\tNum Threads:         %i\n", num_threads);
    printf("\t\tKmer Length:         %i\n", kmer_len);
    printf("\t\tRead Length:         %i\n", read_len);
    
    //Time Vals
    struct timeval ta, tb, tresult; 
    gettimeofday (&ta, NULL); 
    /*Construct taxonomy*/
    get_seqid2taxid(seqid_file, &seqid2taxid);
    construct_taxonomy(taxid_file, my_taxonomy);
    evaluate_kfile(kraken_file, output_file, my_taxonomy, &taxid2node, seqid2taxid, kmer_len, read_len);
    gettimeofday( &tb, NULL);
    timeval_subtract(&tresult, &tb, &ta);
    int minutes = int (tresult.tv_sec / 60);
    tresult.tv_sec = tresult.tv_sec % 60;  
    printf("\tTime Elaped: %i minutes, %i seconds, %0.5f microseconds\n", minutes, tresult.tv_sec, tresult.tv_usec);
    printf("\t=============================\n");
}


/* METHOD: Process command line arguments. */
void parse_command_line(int argc, char **argv) {
    int opt;
    int intval;
    /*Help Message*/
    if (argc > 2 && strcmp(argv[1], "-h") == 0)
        usage(0);
    
    /*Set arguments*/
    static struct option all_options[] = { 
        {"seqid2taxid", required_argument, 0, 'a'},
        {"taxonomy",    required_argument, 0, 'b'},
        {"kraken",      required_argument, 0, 'c'},
        {"output",     required_argument, 0, 'o'},
        {"help",        no_argument, 0, 'h'},
        {"threads",     required_argument, 0, 't'},
        {"kmerlen",     required_argument, 0, 'k'},
        {"readlen",     required_argument, 0, 'l'},
        {0, 0}
        };
    /*Process arguments*/
    int option_index = 0;
    struct stat buffer;
    while ((opt = getopt_long(argc, argv, "h:l:k:t:o:a:b:c:", all_options, &option_index)) != -1) {
        switch(opt) {
            case 'h':
                usage(0);
                break;
            case 'b':
                /*taxonomy FOLDER*/
                taxid_file = optarg;
                if (taxid_file.back() == '/') {
                    taxid_file = taxid_file + "nodes.dmp";
                } else {
                    taxid_file = taxid_file + "/nodes.dmp";
                }
                /*check if file exists*/
                /*if (stat(taxid_file.c_str(), &buffer) != 0) {
                    printf("  %s does not exist",taxid_file.c_str());
                    usage(1);
                }*/
                break;
            case 'a':
                /*seqid2taxid file*/
                seqid_file = optarg;
                /*check if file exists*/
                /*
                if (stat(seqid_file.c_str(), &buffer) != 0) {
                    printf("  %s does not exist", seqid_file.c_str());
                    usage(1);
                }*/
                break;
            case 'o':
                /*Output file name*/
                output_file = optarg;
                break;
            case 'c':
                /*database.kraken file*/
                kraken_file = optarg;
                /*check if file exists*/
                /*
                if (stat(kraken_file.c_str(), &buffer) != 0) {
                    printf("  %s does not exist", kraken_file.c_str());
                    usage(1);
                }*/
                break;
            case 'l':
                /*check negative kmer length*/
                /*do not allow kmer lengths <= 1*/
                read_len = atoi(optarg);
                if (read_len <= 1) {
                    errx(1, "  read lengths must be >= 1\n");
                    usage(1);
                }
                break;
            case 'k':
                /*check negative kmer length*/
                /*do not allow kmer lengths <= 1*/
                kmer_len = atoi(optarg);
                if (kmer_len <= 1) {
                    errx(1, "  kmer lengths must be >= 1\n");
                    usage(1);
                }
                break;
            case 't':
                intval = atoi(optarg);
                /*check negative number of threads*/
                if (intval <= 0) {
                    errx(1, "  can't use nonpositive threads");
                    usage(1);
                }
                if (intval > omp_get_num_procs()) {
                    errx(1, "  thread count exceeds number of processors");
                    usage(1);
                }
                /*set number of threads*/
                num_threads = intval;
                omp_set_num_threads(num_threads);
                break;
            default:
                usage(1);
                break;
        }
    }
    /*Check mandatory options*/
    if (taxid_file == "") {
        printf("  Must specify --taxonomy folder!\n");
        usage(1);
    } else if (seqid_file == "") {
        printf("  Must specify --seqid2taxid file!\n");
        usage(1);
    } else if (kraken_file == "") {
        printf("  Must specify --kraken file! (database.kraken file)\n");
        usage(1);
    } else if (output_file == "") {
        printf("  Must specify --output file!\n");
        usage(1);
    }
    /*check if files exists*/
    //taxid_file = "taxonomy/nodes.dmp";
    ifstream test1(taxid_file.c_str());
    ifstream test2(seqid_file.c_str());
    ifstream test3(kraken_file.c_str());
    if (!test1.is_open()) {
        printf("  %s does not exist",taxid_file.c_str());
        usage(1);
    } else if (!test2.is_open()) {
        printf("  %s does not exist",seqid_file.c_str());
        usage(1);
    } else if (!test3.is_open()) {
        printf("  %s does not exist",kraken_file.c_str());
        usage(1);
    }
}

/* METHOD: Print usage message and exit. */
void usage(int exit_code) {
    if (exit_code == 1) {
        printf("  For usage, please run: \n");
        printf("     kmer2read_distr --help\n");  
        exit(exit_code);
    } 
    cerr << "--------------------------------------------------------------------------" << endl;
    cerr << "Usage: kmer2read_distr [options]" << endl << endl
        << "  *Required Parameters:" << endl
        << "     --seqid2taxid FILE     seqid2taxid file generated during " << endl
        << "                            Kraken database building process" << endl
        << "     --taxonomy FOLDER      taxonomy folder containing the nodes.dmp file" << endl
        << "                            (typically downloaded with the Kraken taxonomy)" << endl
        << "     --kraken FILE          kraken file of all classifications of all library" << endl
        << "                            sequences (typically database.kraken)" << endl
        << "     --output FILE          name of an output file to print read distributions to" << endl
        << "                            (suggested name: databaseXmers.kraken_cnts)" << endl
        << "  *Optional Parameters" << endl
        << "     -k NUM                 kmer length used to build Kraken database" << endl
        << "                            (default = 31)" << endl
        << "     -l NUM                 read length (evaluate every l-length read)" << endl
        << "                            (default = 100)" << endl
        << "     -t NUM                 number of threads" << endl
        << "                            (default = 1)" << endl
        << "  User must specify --seqid2taxid, --taxonomy, --kraken, and --output options" 
        << endl;
    cerr << "---------------------------------------------------------------------------" << endl;
    cerr << endl;
    exit(exit_code);
}

/*METHOD: Use the nodes.dmp to construct the taxonomy!*/ 
//MUST MAKE TAXONOMY HEADER/STRUCTURE
//copy method of building from python folder
void construct_taxonomy(const string t_file, taxonomy *my_taxonomy) {
    /*Initialize variables*/
    int pos1, pos2, pos3;
    int n_count = 0;
    string line;
    int curr_taxid = 0;
    int curr_parent = 0;
    taxonomy *curr_node;
    string rank = "";
    map<int, map<int, taxonomy *>> parent2kids; 
    /*Read through file line by line*/
    ifstream nodefile (t_file);
    if (nodefile.is_open()){
        printf("\t>>STEP 2: READING NODES.DMP FILE\n");
        printf("\t\t0 nodes read");
        while(getline(nodefile, line)) {
            n_count += 1;
            if (n_count % 1000 == 0) 
                printf("\r\t\t%i nodes read", n_count);
            //Find delimiter indices
            pos1 = line.find("\t|\t");
            pos2 = line.find("\t|\t", pos1+1);
            pos3 = line.find("\t|\t", pos2+1);
            //Extract taxid, parent, and rank information
            curr_taxid = atoi(line.substr(0, pos1).c_str());
            curr_parent = atoi(line.substr(pos1+3, pos2-pos1-3).c_str()); 
            rank = line.substr(pos2+3, pos3-pos2-3); 
            //Set Node information
            if (curr_taxid == 1) {
                curr_node = new taxonomy(curr_taxid, rank);
                my_taxonomy = curr_node;
                curr_node->set_lvl_num(1);
                taxid2node[curr_taxid] = curr_node;
            } else {
                if (taxid2node.find(curr_parent) == taxid2node.end()) {
                    //Parent does not exist
                    //Create node and insert into both vectors
                    curr_node = new taxonomy(curr_taxid, rank);
                    taxid2node[curr_taxid] = curr_node;
                    parent2kids[curr_parent][curr_taxid] = curr_node;
                } else {
                    //Parent does exist
                    curr_node = new taxonomy(curr_taxid, rank, taxid2node[curr_parent]);
                    taxid2node[curr_taxid] = curr_node;
                    taxid2node[curr_parent]->add_child(curr_node);
                }
            }
        }
        printf("\r\t\t%i total nodes read\n", n_count);
        nodefile.close();
    } else {
        printf("  cannot open %s", t_file.c_str());  
        usage(1); 
    }
    // printf("\t>>STEP 2.2: LINKING PARENTS/CHILDREN TAXONS\n");
    /*Traverse through parent array and make sure all parents/children are linked*/ 
    map<int, taxonomy *> children;
    taxonomy * parent; 
    for (auto const& pair : parent2kids) {
        curr_parent = pair.first;
        parent = taxid2node[curr_parent];
        children = pair.second;
        for(auto const& child : children) {
            curr_taxid = child.first;
            taxid2node[curr_taxid]->add_parent(parent);
            parent->add_child(taxid2node[curr_taxid]);
        }
    }
    /*Traverse through entire tree and set level numbers */
    //printf("\t>>STEP 2.3: SETTING LEVELS\n");
    n_count = 1;
    //printf("\t\t%i node updated", n_count);
    queue<taxonomy *> curr_nodes;
    //Start with root'schildren
    my_taxonomy->set_lvl_num(1);
    vector<taxonomy *> retrieve_nodes = my_taxonomy->get_children();
    for(taxonomy * add_n: retrieve_nodes) { 
        curr_nodes.push(add_n);
    }
    //Iterate over all children in order
    while (!curr_nodes.empty()) {
        n_count += 1;
        //if (n_count % 1000 == 0) 
        //    printf("\r\t\t%i nodes updated", n_count);
        curr_node = curr_nodes.front();
        curr_nodes.pop();
        curr_node->set_lvl_num(curr_node->get_parent()->get_lvl_num() + 1);
        //Add this node's children
        retrieve_nodes = curr_node->get_children();
        for(taxonomy * add_n: retrieve_nodes) {
            curr_nodes.push(add_n);
        }
    }
    //printf("\r\t\t%i total nodes updated\n", n_count);
}

/*METHOD: Create map of seqids to taxonomy ids from the seqid2taxid file*/
void get_seqid2taxid(string s_file, map<string, int> *seqid2taxid) {
    /*Read through file line by line*/
    int s_count = 0;
    int curr_taxid, pos1, pos2;
    string curr_seqid, line;
    ifstream mapfile (s_file);
    if (mapfile.is_open()){
        printf("\t>>STEP 1: READING SEQID2TAXID MAP\n");
        printf("\t\t0 sequences read");
        while(getline(mapfile, line)) {
            s_count += 1;
            if (s_count % 1000 == 0) 
                printf("\r\t\t%i sequences read", s_count);
            //Find delimiter indices
            pos1 = line.find("\t");
            pos2 = line.find("\n");
            //Extract seqid and taxid 
            curr_seqid = line.substr(0, pos1);
            curr_taxid = atoi(line.substr(pos1+1, pos2-pos1).c_str()); 
            //Save in map
            seqid2taxid->insert(std::pair<string,int>(curr_seqid, curr_taxid));
        }
        printf("\r\t\t%i total sequences read\n", s_count);
    } else {
        printf("  cannot open %s", s_file.c_str());  
        usage(1); 
    }
}

