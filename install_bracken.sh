#!/bin/bash

#####################################################################
#install_bracken.sh checks for all dependencies
#Copyright (C) 2016-2022 Jennifer Lu, jlu26@jhmi.edu
#
#This file is part of Bracken.
#
#Bracken is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or 
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#####################################################################

set -e

VERSION="2.5.3"


if [ ! -z "$@" ];
then
    if [ "$@" = "-v" ]; 
        then echo install_bracken.sh v${VERSION} && exit 0; 
        fi
fi

#Move to Bracken installation directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
cd $DIR
#Allow build script to be run without sh preceding it 
chmod +x "bracken-build"
chmod +x "bracken"

cd src/ && make

echo "Bracken installation complete."
echo 
echo "For simplification, you can symlink the following files"
echo "into a directory in your PATH:"
echo "    bracken"
echo "    bracken-build" 
