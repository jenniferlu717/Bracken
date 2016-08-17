#!/bin/env perl
#count-kmer-abundance.py counts the kmer abundance for genomes classified by Kraken. 
#Copyright (C) 2016 Florian P. Breitwieser, fbreitw1@jhu.edu

#This file is part of Bracken.

#Bracken is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use List::MoreUtils qw/pairwise/;
use File::Basename;
use Parallel::ForkManager;
use Getopt::Long;

my $read_length = 75;
my $n_threads = 1;
my $db_path = ".";
my $opt_help = 0;
my $opt_version = 0;

my $prog = basename($0);
my $version = "0.1";
my $USAGE = "Usage: $prog [options] <database.kraken>

Options:
  --db NAME               Path to Kraken DB (default: $db_path)
  --threads NUM           Number of threads (default: $n_threads)
  --read-length NUM       Read length (default: $read_length)
  --help                  Print this message
  --version               Print version information
";

GetOptions ("db=s" => \$db_path,
            "l|read-length=i" => \$read_length,
            "p|threads=i" => \$n_threads,
            "v|version"  => \$opt_version,
            "h|help" => \$opt_help
) 
or usage(1);

sub usage {
  my $exitval = shift;
  print STDERR $USAGE;
  exit $exitval;
}

sub version {
  print STDERR "$prog version $version\n";
  exit 0;
}

usage(0) if $opt_help;
version() if $opt_version;

if (!defined $ARGV[0] || !-f $ARGV[0]) {
  print STDERR "Give database.kraken file as first command line argument.\n";
  print STDERR "See '$prog --help' for more information.\n";
  exit 1;
}

my $last_real_taxid;
my $last_classified_taxid;
my $current_hashref;
my $current_cntref;
my $n=$read_length-31+1;


my $depth_map = get_taxid_to_rankid();
$depth_map->{"A"} = 0;
$depth_map->{0} = 0;
my $seqid_to_taxid_map = get_seqid_mapping();

my $current_line = 0;
my $n_lines = `cat $ARGV[0] | wc -l`;
chomp $n_lines;

print STDERR "Going through Kraken file ...  ";
my ($seqid,$classified_taxid);

++$|;

my $counti = 0;
my $countb = 0;
my $pm = new Parallel::ForkManager($n_threads);
while (my $line = <>) {
    printf STDERR "\r%3.2f%%; line %10s / %s",(++$current_line/$n_lines*100),$current_line,$n_lines;
    $pm->start and next; # do the fork
    my (undef,$seqid,$classified_taxid,$length,$classification_string) = split(/\t/, $line);
    defined $classification_string or die "ERROR: $line ?\n";
    if ($length < $read_length) {
        $pm->finish; next;
    }
    my $real_taxid;
    if ($seqid =~ /^kraken:taxid/) {
        $seqid =~ /kraken:taxid\|([0-9]+)[^\t]*/;
        $real_taxid = $1;
    } else {
        if (defined $seqid_to_taxid_map->{$seqid}) {
            $real_taxid = $seqid_to_taxid_map->{$seqid};
        }     
    }

    if (defined $real_taxid) {
        
        my %real_taxid_cnts;
    
        my @classified_taxids;
        $classified_taxids[$n-1] = 0;
        
        my ($classified_taxids_len,$classified_taxids_i) = (0,0);
        my ($max_depth_idx,$current_max_depth,$current_classified_taxid) = (0,0,0);
    
        while (($classification_string =~ m/([0-9]+):([0-9]+)/g)) {
            my ($taxid,$count) = ($1,$2);
            ++$countb;
            my $taxid_depth = $depth_map->{$taxid};
            #if ($taxid eq 'A') {
            #    @classified_taxids = ();
            #    $classified_taxids[$n-1] = 0;
            #    ($classified_taxids_len,$classified_taxids_i) = (0,0);
            #    ($max_depth_idx,$current_max_depth,$current_classified_taxid) = (0,0,0);
            #    next;
            #}
            defined $taxid_depth or die "Couldn't find taxid $taxid in nodes.dmp (count $count) in line $current_line: ".substr($line,0,50)." ...\n";
            
            ## fill up the string if classified_taxids_len < n
            ## if count is greater than zero, overwrite the positions and update the count
            while ($count > 0) {
                $classified_taxids[$classified_taxids_i] = $taxid;
                --$count; ++$classified_taxids_len;
                next if $classified_taxids_len < $n;
    
                if ($taxid_depth > $current_max_depth) {
                    $current_classified_taxid = $taxid;
                    $max_depth_idx = $classified_taxids_i;
                    $current_max_depth = $taxid_depth;
                } elsif ($taxid_depth < $current_max_depth) {
                    if ($max_depth_idx == $classified_taxids_i) {
                        for (my $i=0; $i <= $#classified_taxids; ++$i) {
                            if ($depth_map->{$classified_taxids[$i]} > $depth_map->{$classified_taxids[$max_depth_idx]}) {
                                $max_depth_idx = $i;
                            }
                        }
                    }
                    $current_classified_taxid = $classified_taxids[$max_depth_idx];
                    die "depth for taxid $current_classified_taxid at idx $max_depth_idx is undefined" unless defined $depth_map->{$current_classified_taxid};
                    $current_max_depth = $depth_map->{$current_classified_taxid};
                } else {
                    $max_depth_idx = $classified_taxids_i;
                }
                ++$real_taxid_cnts{$current_classified_taxid};
                ++$counti;
                ++$classified_taxids_i;
                $classified_taxids_i -= $n if $classified_taxids_i >= $n;
            }
        }
        print $seqid,"\t",$real_taxid,"\t",$classified_taxid,"\t",$length,"\t",join(" ",map{ $_.":".$real_taxid_cnts{$_} } sort {$real_taxid_cnts{$b} <=> $real_taxid_cnts{$a}} keys %real_taxid_cnts),"\n";
    } else {
        print "$seqid\tno mapping\t$classified_taxid\t$length\t0\n";
    }
    $pm->finish; # do the exit in the child process

}

$pm->wait_all_children;

print STDERR "Done: $counti; $countb.\n";

sub get_seqid_mapping {
  my %seqid_to_taxid;
  open (my $S, "<", "$db_path/seqid2taxid.map")
      or die "Can't open seqid2taxid.map file - make sure the Kraken database path is defined correctly. $!\n";
  while (<$S>) {
    chomp;
    my ($seqid,$taxid) = split(/\t/);
    $seqid_to_taxid{$seqid} = $taxid;
  }
  close($S);
  return \%seqid_to_taxid;
}

sub get_taxid_to_rankid {
  my $child_lists;
  my %depth_map;
  #if (-f "$prefix/taxonomy/child_lists.perlhash" && -f "$prefix/taxonomy/rank_map.perlhash") {
  #  $rank_map = retrieve("$prefix/taxonomy/rank_map.perlhash");
  #} else {
  ++$|;
  print STDERR "Reading nodes.dmp ... "; 
    open NODES, "<", "$db_path/taxonomy/nodes.dmp"
      or die "Can't open nodes file - make sure the Kraken database path is defined correctly. $!\n";
    while (<NODES>) {
      my ($node_id, $parent_id) = split(/\t\|\t/);
      $child_lists->{$parent_id} ||= [];
      push @{ $child_lists->{$parent_id} }, $node_id;
    }
    close NODES;
  print STDERR "Done\n";
  print STDERR "Calculating taxon depths ... "; 
  get_rank_depth($child_lists, \%depth_map, 1, 1);
  print STDERR "Done\n";
  return \%depth_map;
}

sub get_rank_depth {
    my ($child_lists, $depth_map, $depth, $taxon) = @_;
    $depth_map->{$taxon} = $depth;
    return unless defined $child_lists->{$taxon};
    foreach my $child (@{ $child_lists->{$taxon} }) {
        next if $child == $taxon;
        die "Depth already defined for $child!" if defined $depth_map->{$child};
        get_rank_depth($child_lists, $depth_map, $depth+1, $child);
    }
}
