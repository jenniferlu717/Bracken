/*********************************************************************
 * taxonomy.h is used as part of the kmer2distr script
 * Copyright (C) 2016-2019 Jennifer Lu, jlu26@jhmi.edu
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
