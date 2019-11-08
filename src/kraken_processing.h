/*********************************************************************
 * kraken_processing.h is used as part of the kmer2distr script
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
#ifndef KRAKEN_PROCESSING_H
#define KRAKEN_PROCESSING_H

#include "kmer2read_headers.h"
#include "taxonomy.h"
#include "time.h"
void evaluate_kfile(string, string, const taxonomy *, const map<int, taxonomy *> *, map<string, int>, const int, const int);

int read_kfile(string, const map<string, int> *, map<int, string> *, map<int, string> *, map<int, int> *);

void convert_distribution(string, int, const map<int, string> *, const map<int, string> *, const map<int, int> *, const taxonomy *, const map<int, taxonomy *> *, const int, const int);

int get_classification(vector<int> *, const taxonomy *, const map<int, taxonomy *> *);

#endif
