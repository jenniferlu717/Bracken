#include "taxonomy.h"
using std::vector;
using std::string;

/*Null constructor*/
taxonomy::taxonomy() {
    this->taxid = -1;
    this->lvl_type = "N";
    this->parent = NULL;
}
/*Constructor without Parent*/
taxonomy::taxonomy(int taxid, string level_type) {
    this->taxid = taxid;
    this->lvl_type = level_type;
    this->parent = NULL;
}
/*Constructor for the Taxonomy Node*/
taxonomy::taxonomy(int taxid, string level_type, taxonomy *parent) {
    this->taxid = taxid;
    this->lvl_type = level_type;
    this->parent = parent;
}
/*Destructor */ 
taxonomy::~taxonomy() {
}

/*Accessor Methods*/

int taxonomy::get_taxid() const {
    return this->taxid; 
}

int taxonomy::get_lvl_num() const {
    return this->lvl_num;
}

string taxonomy::get_lvl_type() const {
    return this->lvl_type;
}

taxonomy* taxonomy::get_parent() const{
    return this->parent;
}

vector<taxonomy *> taxonomy::get_children() const {
    return this->children;
}

/*Methods for manipulating the tree*/
void taxonomy::add_parent(taxonomy *parent) {
    this->parent = parent;
}

void taxonomy::add_child(taxonomy *new_child) {
    this->children.push_back(new_child);
}

void taxonomy::set_lvl_num(int num) {
    this->lvl_num = num;
}
