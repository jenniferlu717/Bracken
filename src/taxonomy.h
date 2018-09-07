#ifndef TAXONOMY_H
#define TAXONOMY_H
#include "kmer2read_headers.h"

/* Class defining the  Taxonomy.
 * This class specifies the methods and variables
 * used for adding to/accessing the taxonomy tree node
 */
class taxonomy{
    public:
        /*Constructor and Destructor*/
        taxonomy();
        taxonomy(int, string);
        taxonomy(int, string, taxonomy *);
        ~taxonomy();
        /*Methods for accessing this node's information*/
        int get_taxid() const;
        int get_lvl_num() const;
        string get_lvl_type() const; 
        taxonomy* get_parent() const;
        vector<taxonomy *> get_children() const;
        /*Methods for manipulating the tree*/
        void add_parent(taxonomy *);
        void add_child(taxonomy *);
        void set_lvl_num(int);
        /*Other methods for comparisons*/
        void get_lca(taxonomy *); 
    private:
        int taxid;
        int lvl_num;
        string lvl_type;
        vector<taxonomy *> children;
        taxonomy *parent;  
};
#endif
