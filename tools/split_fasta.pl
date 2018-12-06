#!/usr/bin/perl

######################## Copyright and License #########################
#
# splitfasta.pl - Splits fasta file into multiple files of equal size.
# Copyright 2012 Fabian Kilpert, Bank Beszteri, Adam Monier
#
# This file is part of PhyloAssigner.
#  
# PhyloAssigner is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# PhyloAssigner is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PhyloAssigner.  If not, see <http://www.gnu.org/licenses/>.
#
########################################################################


use strict;
use warnings;

use Getopt::Long;
use File::Spec;
use File::Basename;


my $usage = "splitfasta.pl [options] seq_file 
    Split fasta file into multiple files of equal size.

    seq_file  : fasta file
    
    Options:
    -n number      : number of output files (default: 2)
    -v             : verbose mode
    
";


#- default values ------------------------------------------------------
my $parts = 2;      # number of output files
my $verbose;        # verbose mode


#- command line handling -----------------------------------------------

#Options
GetOptions ('number=i' => \$parts, 
            'verbose' => \$verbose) or die $usage;

#Arguments
my $infile = shift or die $usage;
$infile = File::Spec->rel2abs($infile);
my $dirname = dirname($infile);      #main dir of infile
my $basename = basename($infile);       #basename of infile

my $basename_pure;
my $suffix;
if ($basename =~ m/(.*)(\.\w*)$/) {
    $basename_pure = $1;
    $suffix = $2;
}

#- Handle output files -------------------------------------------------

my @out_files;
my @out_file_handlers;

print "Writing sequences to $parts output files:\n" if $verbose;

#for (my $j=0; $j < $parts; $j++) {
for(0..$parts-1) {
    my $part_name = sprintf ("${dirname}/${basename_pure}.part%.2i${suffix}", $_);    # make output file name
    push (@out_files, $part_name);
    open (my $fh_out, '>', $part_name) or die "Couldn't open file \"$part_name\": $!";
    push (@out_file_handlers, $fh_out);
    print $part_name, "\n" if $verbose;
}


#- file input ----------------------------------------------------------

open(my $fh, '<', $infile) or die "Couldn't open file \"$infile\": $!";

my $header;
my $sequence;
my $i = 1;
my $first_header = 1;
my $cur = 0;            # current active file number
my $cur_fh;             # current active file handler

while(my $line = <$fh>) {
    $line =~ s/\r\n/\n/g;      #replace CRLF by LF only
    chomp $line;
    
    #Get header
    if ($line =~ m/^>(.+)/) {
        if ($first_header) {
            $header = $1;
            
            $cur_fh = $out_file_handlers[$cur];
            print {$cur_fh} '>', $header, "\n";
            
            $first_header = 0;
            ++$i;
        }
        else {
            $sequence = &format_sequence_for_output($sequence);
            print {$cur_fh} $sequence;
            $sequence = '';
            
            if ($cur >= scalar @out_files - 1) {
                $cur = 0;
                $cur_fh = $out_file_handlers[$cur];
            }
            else {
                $cur += 1;
                $cur_fh = $out_file_handlers[$cur];
            }
            
            $header = $1;
            print {$cur_fh} '>', $header, "\n";
            ++$i;
        }
    }
    # Get sequence
    else {
        $line =~ s/\s//g; 
        $sequence .= $line;
    }
}

$sequence = &format_sequence_for_output($sequence);

print {$cur_fh} $sequence;

close $fh;  #close input file


#- Close all output files ----------------------------------------------

foreach my $ofh (@out_file_handlers) {
    close $ofh;    
}

print "Sequences processed: ",  $i-1, "\n" if $verbose;


#- FASTA output formatting subroutine ----------------------------------

sub format_sequence_for_output {
    my $sequence = shift @_;
    my $width = 50;
    my $formatted_sequence;
    
    for (my $pos=0; $pos < length $sequence; $pos += $width) {
        my $subsequence = substr($sequence, $pos, $width);
        $formatted_sequence .= $subsequence."\n";
    }
    return $formatted_sequence;
}
