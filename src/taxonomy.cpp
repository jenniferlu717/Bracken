/*********************************************************************
 * taxonomy.cpp is used as part of the kmer2distr script
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
 * Updated: 2018/09/06
 */

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
