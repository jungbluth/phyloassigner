#!/usr/bin/perl

######################## Copyright and License #########################
#
# fasta2stockholm.pl - Convert from FASTA to STOCKHOLM format
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


my $usage = "fasta2stockholm.pl fasta_file
    Convert file from FASTA into STOCKHOLM format.\n\n";


#Arguments
my $infile = shift @ARGV or die $usage;
$infile = Cwd::abs_path($infile);

#print $infile, "\n";


#- read infile ---------------------------------------------------------
open (my $fh, "<", $infile) or die "Could not open $infile\n$!\n";

my @seqs;
my %seq;
my $header = '';
my $sequence = '';
my $i = 0;
my $first_header = 1;

while (my $line = <$fh>) {
    $line =~ s/\s*$//g;     #replace one or more linefeeds
    
    #Get header
    if ($line =~ m/^>(.+)/) {
        ++$i;
        
        if ($first_header) {
            $header = $line;
            $header =~ s/^>*//g;
            
            #print $i, "\n";
            #print $header, "\n";
            
            $first_header = 0;
        }
        else { 
            #print $sequence, "\n";
            
            $seq{$header} = $sequence;
            push(@seqs, $header);
            
            $header = $line;
            $header =~ s/^>*//g;
            
            #print $i, "\n";
            #print $header, "\n";
            
            $sequence = '';
        }
    }
    else {
        $line =~ s/\s//g;
        $sequence .= $line;
    }

}
close $fh;

$seq{$header} = $sequence;
push(@seqs, $header);



#- print out STOCKHOLM format ------------------------------------------

print "# STOCKHOLM 1.0\n";

foreach (@seqs) {
    printf "%s %s\n", $_, $seq{$_};
}

print "//\n";

