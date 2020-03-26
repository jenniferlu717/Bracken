# Bracken 2.5.3 abundance estimation
For Bracken news, updates, and instructions: https://ccb.jhu.edu/software/bracken/ 

Bracken's peer-reviewed paper (published Jan 2, 2017): https://peerj.com/articles/cs-104/

# Installation
Bracken is a companion program to Kraken 1 or Kraken 2
While Kraken classifies reads to multiple levels in the taxonomic tree, 
Bracken allows estimation of abundance at a single level using those classifications
(e.g. Bracken can estimate abundance of species within a sample). 

Prior to installing Bracken, please install Kraken:
Kraken can be downloaded from here: http://ccb.jhu.edu/software/kraken/

## Easy Bracken Installation: 
    sh install_bracken.sh
## Hard Bracken Installation:
    cd src/ && make
    Add bracken/bracken-build and scripts in src/ to your PATH 

# Bracken 2.5.3 Changes
Bracken 2.5.3 has small changes in options to allow for 
non-traditional abundance estimation (e.g. S1, G1, K7) and 
allows specification for naming the kraken-style report 
with bracken read counts. 

# Bracken 2.5 Changes
Bracken 2.5 has a 30x faster build-time. Previously, 1 million 
database sequences would take hours to process by Bracken. Now, 
the same process requires less than 10 minutes (with 16 threads).
All output is identical to previous output. 

# Bracken 2.0 Changes
Bracken 2.0 does not introduce any changes in the main algorithm. 
[Output from Bracken 1.0 is identical to output from Bracken 2.0]
However, additional scripts are provided to allow easier/faster 
installation and execution of the Bracken code. 

Users can either run the higher-level scripts now provided, OR 
follow the individual steps described.

# RUNNING BRACKEN: EASY VERSION
Steps 0/1 are run once per database. If you would like to generate
Bracken files for multiple read lengths, repeat Step 1 specifying the same
database but different read lengths. The script will skip any step already complete. 

If you run Kraken using one of the pre-built MiniKraken databases, you can find 
corresponding Bracken files [here](https://ccb.jhu.edu/software/bracken/). 
Do not run bracken-build with MiniKraken. 

## Step 0: Build a Kraken 1.0 or Kraken 2.0 database
        kraken-build --db ${KRAKEN_DB} --threads ${THREADS}
        kraken2-build --db ${KRAKEN_DB} --threads ${THREADS} 
   
   * `${KRAKEN_DB}` is the path to a built Kraken database which also must contain:
        * the taxonomy/nodes.dmp file
        * and library sequences `*.fna`, `*.fa`, or `*.fasta` in the `library` directory.

## Step 1: Generate the Bracken database file (databaseXmers.kmer_distrib)  
   It is highly encouraged for users to run the following scripts with 10-20 threads.
     [if run single-threaded, kraken/kraken2 and kmer2read_distr will take hours-days]
   Please note that the flags for this script are single lettered
      
   If Kraken 1.0 or Kraken 2.0 is included in your PATH, run the following
        
        bracken-build -d ${KRAKEN_DB} -t ${THREADS} -k ${KMER_LEN} -l ${READ_LEN}

   Otherwise, direct the program using "-x" to the installation/location of the ./kraken or ./kraken2 scripts
        
        bracken-build -d ${KRAKEN_DB} -t ${THREADS} -k ${KMER_LEN} -l ${READ_LEN} -x ${KRAKEN_INSTALLATION}

            `${KRAKEN_DB}`  = location of the built Kraken 1.0 or Kraken 2.0 database
            `${THREADS}`    = number of threads to use with Kraken and the Bracken scripts
            `${KMER_LEN}`   = length of kmer used to build the Kraken database 
                                    Kraken 1.0 default kmer length = 31
                                    Kraken 2.0 default kmer length = 35
                                    Default set in the script is 35. 
            `${READ_LEN}`   = the read length of your data 
                                    e.g., if you are using 100 bp reads, set it to `100`. 

## Step 2: Run Kraken 1.0 or Kraken 2.0 AND Generate a report file 
   Kraken 1.0 requires a 2-step process to generate the report file needed by Bracken
        
        kraken --db ${KRAKEN_DB} --threads ${THREADS} ${SAMPLE}.fq > ${SAMPLE}.kraken
        kraken-report --db ${KRAKEN_DB} ${SAMPLE}.kraken > ${SAMPLE}.kreport 

   Kraken 2.0 requires the addition of the --report flag 
        
        kraken2 --db ${KRAKEN_DB} --threads ${THREADS} --report ${SAMPLE}.kreport ${SAMPLE}.fq > ${SAMPLE}.kraken

## Step 3: Run Bracken for Abundance Estimation
        
        bracken -d ${KRAKEN_DB} -i ${SAMPLE}.kreport -o ${SAMPLE}.bracken -r ${READ_LEN} -l ${LEVEL} -t ${THRESHOLD}


# RUNNING BRACKEN: HARD VERSION
## Step 0: Build a Kraken 1.0 or Kraken 2.0 database
        kraken-build --db ${KRAKEN_DB} --threads ${THREADS} 
        kraken2-build --db ${KRAKEN_DB} --threads ${THREADS}
   
   * `${KRAKEN_DB}` is the path to a built Kraken database which also must contain:
        * the taxonomy/nodes.dmp file
        * and library sequences `*.fna`, `*.fa`, or `*.fasta` in the `library` directory.

## Step 1: Generate the Bracken database file (databaseXmers.kmer_distrib)  
   * It is highly encouraged for users to run the following scripts with 20 threads.
### Step 1a: Search all library input sequences against the database
Run the following scripts WITHIN the Kraken database folder: 

        kraken --db=${KRAKEN_DB} --threads=10 <( find -L library \(-name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) -exec cat {} + )  > database.kraken
        kraken2 --db=${KRAKEN_DB} --threads=10 <( find -L library \(-name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) -exec cat {} + )  > database.kraken

### Step 1b: Compute classifications for each perfect read from one of the input sequences

        /src/kmer2read_distr --seqid2taxid ${KRAKEN_DB}/seqid2taxid.map --taxonomy ${KRAKEN_DB}/taxonomy --kraken database.kraken --output database${READ_LEN}mers.kraken
            -k ${KMER_LEN} -l ${READ_LEN} -t ${THREADS}
            
            `${KRAKEN_DB}`  = location of the built Kraken 1.0 or Kraken 2.0 database
            `${THREADS}`    = number of threads to use [recommended: 20]
            `${KMER_LEN}`   = length of kmer used to build the Kraken database 
                                    Kraken 1.0 default kmer length = 31
                                    Kraken 2.0 default kmer length = 35
                                    [default: 35]
            `${READ_LEN}`   = the read length of your data 
                                    e.g., if you are using 100 bp reads, set it to `100`. 

### Step 1c: Generate the kmer distribution file
The kmer distribution file is generated using the following command line:

    python generate_kmer_distribution.py -i database${READ_LEN}mers.kraken -o database${READ_LEN}mers.kmer_distrib
    
## Step 2: Run Kraken 1.0 or Kraken 2.0 AND Generate a report file 

Kraken 1.0 requires a 2-step process to generate the report file needed by Bracken
        
        kraken --db ${KRAKEN_DB} --threads ${THREADS} ${SAMPLE}.fq > ${SAMPLE}.kraken
        kraken-report --db ${KRAKEN_DB} ${SAMPLE}.kraken > ${SAMPLE}.kreport 

Kraken 2.0 requires the addition of the --report flag 
        
        kraken2 --db ${KRAKEN_DB} --threads ${THREADS} --report ${SAMPLE}.kreport ${SAMPLE}.fq > ${SAMPLE}.kraken

## Step 3: Run Bracken for Abundance Estimation
Given the expected kmer distribution for genomes in a kraken database along
with a kraken report file, the number of reads belonging to each species (or
genus) is estimated using the estimate_abundance.py file, run with the
following command line:

    python estimate_abundance.py -i ${SAMPLE}.kreport -k database${READ_LEN}mers.kmer_distrib -l ${CLASSIFICATION_LVL} -t ${THRESHOLD} -o ${BRACKEN_OUTPUT_FILE}.bracken

The following required parameters must be specified:
- `${SAMPLE}`.kreport - the kraken report generated for a given dataset 
- database`${READ_LEN}`mers.kmer_distrib - the file generated by generate_kmer_distribution.py 
- `{BRACKEN_OUTPUT_FILE}`.bracken - the desired name of the output file to be generated by the code

The following optional parameters may be specified:
- `${CLASSIFICATION_LVL}` - Default = `'S'`. This specifies that abundance
    estimation will calculate estimated reads for each species. Other possible
    options are K (kingdom level), P (phylum), C (class), O (order), F (family),
    and G (genus). 
- `${THRESHOLD}` - Default = `10`. For species classification, any species
    with <= 10 (or otherwise specified) reads will not receive any additional reads
    from higher taxonomy levels when distributing reads for abundance estimation.
    If another classification level is specified, thresholding will occur at
    that level.  

# Output Kraken-Style Bracken Report

By default, this script will also recreate the report file using the new Bracken numbers. 
 1. The new report file will be found in the same folder as the original report file, with "bracken" included in the name.
 2. Levels below the estimate-level will not be printed. 
 3. Any levels whose reads were below the threshold will not be included
 4. Percentages will be re-calculated for the remaining levels
 5. Unclassified reads will not be included in the report.  

# Example abundance estimation
The following sample input and output files are included in the sample_data/ folder: 
    `sample_test.report` - Kraken report file generated from the kraken-report command. 
    `sample_kmer_distr_75mers.txt` - example kmer distribution file generated by generate_kmer_distribution.py
    `sample_output_species_abundance.txt` - Bracken species abundance estimation for sample_test.report 
    `sample_output_bracken.report` - Kraken report style file with all reads redistributed to the species level

Due to size constraints, the following files are not included in the sample_data/ folder:
    `sample_test.kraken` - Kraken output file used to generate the Kraken report file
    `database.kraken` - Initial Kraken classification of every genome
    `database75mers.kraken_cnts` - Counting of kmer abundances

The following commands were used to generate each individual file: 

1. ```
    kraken --db${KRAKEN_DB} --threads=10 sample.fa > sample_test.kraken
    kraken-report --db=${KRAKEN_DB} sample_test.kraken > sample_test.report 
   ```
2. ```
    kraken --db=${KRAKEN_DB} --fasta_input --threads=10 <( find -L library -name "*.fna" -o -name "*.fa" -o -name "*.fasta" -exec cat {} + ) > database.kraken 
    perl count-kmer-abubndances.pl --db=${KRAKEN_DB} --read-length=75 database.kraken > database75mers.kraken_cnts
   ```
3. ```
    python generate_kmer_distribution.py -i database75mers.kraken_cnts -o sample_kmer_distr_75mers.txt
   ```
4. ```
    python estimate_abundance.py -i sample_test.report -k sample_kmer_distr_75mers.txt -l S -t 10 -o sample_output_species_abundance.txt 
   ```
# Copyright and licensing
Copyright (C) 2020 Jennifer Lu, jlu26@jhmi.edu

Bracken is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the license, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
    
# Author information
Jennifer Lu (jlu26@jhmi.edu, ccb.jhu.edu/people/jennifer.lu)

Florian Breitwieser (fbreitw1@jhu.edu, ccb.jhu.edu/people/florian)

Last Updated On: 03/26/2020
