/*********************************************************************
 * taxonomy.cpp is used as part of the kmer2distr script
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

#include "taxonomy.h"
using std::vector;
using std::string;

/*Null constructor*/
taxonomy::taxonomy() {
    this->taxid = -1;
    this->lvl_type = "N";
    this->parent = NULL;
    this->lvl_num = 0;
}
/*Constructor without Parent*/
taxonomy::taxonomy(int taxid, string level_type) {
    this->taxid = taxid;
    this->lvl_type = level_type;
    this->parent = NULL;
    this->lvl_num = 0;
}
/*Constructor for the Taxonomy Node*/
taxonomy::taxonomy(int taxid, string level_type, taxonomy *parent) {
    this->taxid = taxid;
    this->lvl_type = level_type;
    this->lvl_num = 0;
    this->parent = parent;
}
/*Destructor */ 
taxonomy::~taxonomy() {
}
