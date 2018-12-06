#!/usr/bin/perl

######################## Copyright and License #########################
#
# parse_hmmer3_outfile.pl - Read HMMER3 outfile, write FASTA
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

use Cwd;

use Bio::SeqIO;


my $usage = "parse_hmmer3_outfile.pl fasta_file
Parse FASTA alignement (AFA) outfile from hmmalign (HMMER3)
(replace . by -, make every sequence character uppercase)
and output FASTA format (to stdout).\n\n";


#Arguments
my $infile = shift @ARGV or die $usage;   #hmmalign output file
$infile = Cwd::abs_path($infile);


# create one SeqIO object to read in
my $seqin = Bio::SeqIO->new(
                            -file   => "<$infile",
                            -format => "FASTA",
                            );


# write each entry in the input to the output
while (my $inseq = $seqin->next_seq) {
    #header
    print ">", $inseq->id, "\n";
    
    #sequence
    my $seq = $inseq->seq;
    $seq =~ s/\./-/g;     # replace "." by "-"
    $seq =~ tr/a-z/A-Z/;  # UPPERCASE
    print $seq, "\n";
}


