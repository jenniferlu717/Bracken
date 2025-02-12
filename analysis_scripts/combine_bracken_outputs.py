#! /usr/bin/env python
"""
#####################################################################
# combine_bracken_output.py combines multiple Bracken output files for comparison
# Copyright (C) 2016-2023 Jennifer Lu, jlu26@jhmi.edu

# This file is part of Bracken.

# Bracken is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the license, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

#####################################################################
# Jennifer Lu, jlu26@jhmi.edu
# 02/21/2018
#
# This program takes multiple Bracken OUTPUT files and creates a single
# text tab-delimined file with estimated abundance across all samples
# for comparison.
#
# Input files:
#   - Bracken output files (not Bracken report files)
#       [NOTE: All samples being compared must have the same level estimated]
#
# Optional Input Parameter:
#   - Sample Names for each sample separated by commas.
#     Number of sample names must match the number of samples
#     If not given, basename of sample files will be used for column headers
#       [e.g. --names S1,S2 will yield column headers S1_num S1_frac S2_num S2_frac]
#
# Output File format (tab-delimited, for species abundance estimation)
#   - Name of species
#   - Taxonomy ID for that species
#   - Taxonomy Level specified for
#   - Sample #1: Number of reads estimated for that species
#   - Sample #1: Fraction of total reads in the sample estimated for this species
#       [NOTE: Fractions are of total reads classified, unclassified reads not accounted for]
#   - Sample #2: Number of reads
#   - Sample #2: Fraction of total reads
#   - ...etc.
#
# Methods:
#   - main
#####################################################################
"""

import argparse
import os
import sys
from time import gmtime, strftime


def get_arguments():
    """Parse arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--files",
        dest="files",
        nargs="+",
        type=str,
        required=True,
        help="Bracken output files to combine.",
    )
    parser.add_argument(
        "-n",
        "--names",
        dest="names",
        default="",
        required=False,
        help=(
            "Names for each input file - to be used in column headers of output "
            "[separate names with commas]"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="Name of output file with combined Bracken results.",
    )
    return parser.parse_args()


def get_sample_names(args):
    """Rename samples or use file names"""
    all_samples = []
    if len(args.names) == 0:
        all_samples = [os.path.basename(f) for f in args.files]
    else:
        all_samples = [args.names.split(",")]
    return all_samples


def read_files(args):
    """Read all input_files"""
    # Initialize variables
    sample_counts = {}  # species :: sample1: counts, samples2: counts
    all_samples = get_sample_names(args)
    total_reads = {sample: 0 for sample in all_samples}

    # Read each file information in
    level = ""
    i = 0
    for f in args.files:
        # Print update
        curr_name = all_samples[i]
        i += 1
        sys.stdout.write(f"Processing Output File {f}:: Sample {curr_name}\n")
        # Iterate through file
        header = True
        with open(f, "r", encoding="utf-8") as i_file:
            for line in i_file:
                # Header line
                if header:
                    header = False
                    continue
                # Process line
                # name, taxid, taxlvl, kreads, areads, estreads, frac
                [name, taxid, taxlvl, _, _, estreads, _] = line.strip().split(
                    "\t"
                )
                estreads = int(estreads)
                # Error Checks
                if name not in sample_counts:
                    sample_counts[name] = {}
                    sample_counts[name][taxid] = {}
                elif taxid != list(sample_counts[name].keys())[0]:
                    sys.exit(
                        f"Taxonomy IDs not matching for species {name}: "
                        f"({taxid}\t{sample_counts[name].keys()[0]})"
                    )
                if len(level) == 0:
                    level = taxlvl
                elif level != taxlvl:
                    sys.exit("Taxonomy level not matching between samples")
                # Save counts
                total_reads[curr_name] += estreads
                sample_counts[name][taxid][curr_name] = estreads
            # Close file
            i_file.close()

    return sample_counts, total_reads, all_samples, level


def process_and_write_output(args, sample_counts, total_reads, all_samples, level):
    """Do the dict processing and file writing"""
    with open(args.output, "w", encoding="utf-8") as o_file:
        o_file.write("name\ttaxonomy_id\ttaxonomy_lvl")
        for name in all_samples:
            o_file.write(f"\t{name}_num\t{name}_frac")
        o_file.write("\n")
        # Print each sample
        for name in sample_counts:
            # Print information for classification
            taxid = list(sample_counts[name].keys())[0]
            o_file.write(f"{name}\t{taxid}\t{level}")
            # Calculate and print information per sample
            for sample in all_samples:
                if sample in sample_counts[name][taxid]:
                    num = sample_counts[name][taxid][sample]
                    perc = float(num) / float(total_reads[sample])
                    o_file.write(f"\t{num}\t{perc:0.5f}")
                else:
                    o_file.write("\t0\t0.00000")
            o_file.write("\n")
        o_file.close()



def main():
    """Main method"""

    args = get_arguments()

    # Start program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write(f"PROGRAM START TIME: {time}\n")

    sample_counts, total_reads, all_samples, level = read_files(args)

    # Print output file header
    process_and_write_output(args, sample_counts, total_reads, all_samples, level)
    # End program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write(f"PROGRAM END TIME: {time}\n")


if __name__ == "__main__":
    main()
