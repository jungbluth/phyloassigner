#!/usr/bin/perl

######################## Copyright and License #########################
#
# dos_to_unix.pl - Convert DOS (CRLF) to UNIX (LF)
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



my $usage = "dos_to_unix.pl [options] file 

    Check fasta file for issues, e.g. id duplicates, unallowed characters.

    file              : input file
    
    Options:
    -o outfile        : output file (may be the same as input file)
    -v                : verbose
    
";


my $verbose;
my $outfile;


#- Get options and arguments -------------------------------------------
#Options
GetOptions ('v' => \$verbose,
            'o=s' => \$outfile)
            or die $usage;

#Arguments
my $infile = shift or die $usage;


# generate outfile name
$outfile = $infile.".qids" if not $outfile;
$outfile = File::Spec->rel2abs($outfile);
print "Outfile: $outfile\n" if $verbose;
-e $infile or die "Could NOT open file! $!\n";


#=begin comment


#- Read file -----------------------------------------------------------
open (my $in_fh, '<', $infile) or die "Unable to open file: \"$infile\" $!\n";
my @lines = <$in_fh>;
close $in_fh;


#- Process lines -------------------------------------------------------
foreach (@lines) {
    $_ =~ s/\s+$//;      #delete one or more linefeeds
}

#print "@lines";


#- Write back to file --------------------------------------------------
open (my $out_fh, '>', $outfile) or die "Unable to open file: \"$outfile\" $!\n";
print {$out_fh} @lines;
close $out_fh;


#=end comment


