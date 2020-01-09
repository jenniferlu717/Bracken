/*********************************************************************
 * kmer2read_headers.h provides the necessary headers for the kmer2distr script
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

#include <omp.h>
#include <map>
#include <deque>
#include <string>
#include <string.h>
#include <vector>
#include <queue>
#include <algorithm>

#include <getopt.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <unistd.h>
using std::atoi; 
using std::min;
using std::max;
using std::endl;
using std::cerr;
using std::cout;
using std::deque; 
using std::queue;
using std::vector;
using std::map;
using std::ifstream;
using std::ofstream;
using std::string;
