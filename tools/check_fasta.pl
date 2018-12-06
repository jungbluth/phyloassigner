#!/usr/bin/perl

######################## Copyright and License #########################
#
# check_fasta.pl - Check FASTA file for not allowed characters
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
use Bio::SeqIO;


#- Default values ------------------------------------------------------
my $usage = "check_fasta.pl [options] seq_file 

    Check fasta file for issues, e.g. id duplicates, unallowed characters.

    seq_file          : fasta file
    
    Options:
    -v                : verbose
    
";


my $verbose;


#- Get options and arguments -------------------------------------------
#Options
GetOptions ('v' => \$verbose)
            or die $usage;
#print "Correction mode...\n" if $correct_issues;

#Arguments
my $file = shift or die $usage;

-e $file or die "Could NOT open file! $!\n";


#-----------------------------------------------------------------------
# Input seq stream
my $seq_in = Bio::SeqIO->new(-file => "<$file",
                            -format => 'fasta');
                            
my @seq_list;    # list of all sequences (sequence objects) in the file
my @bad_id_list;
my @bad_seq_list;
my %id_dic;
my @duplicate_id_list;


# Store all seqs in list
while (my $seq = $seq_in->next_seq()) {
    push (@seq_list, $seq);
}

# Info about data
(printf "\nSequences in file: %i\n" , scalar @seq_list) if $verbose;


#- Check for duplicate IDs ---------------------------------------------
foreach my $seq (@seq_list) {
    if (exists $id_dic{$seq->id}) { 
        $id_dic{$seq->id} += 1; 
        push (@duplicate_id_list, $seq); }
    else { $id_dic{$seq->id} = 1; }
}

#output
if (@duplicate_id_list) {
    print "\nWARNING! Duplicate IDs:\n";
    
    while (my ($k, $v) = each %id_dic) {
        printf "%s : %s\n", $k, $v if $v != 1;
    }
}


#- Check for non-alphanumeric characters in id -------------------------
foreach (@seq_list) {
    push(@bad_id_list, $_) if ($_->id =~ /[^\w]/);
}

#output
if (@bad_id_list) {
    printf "\nWARNING! Non-alphanumerical characters detected in %i IDs:\n", scalar @bad_id_list;
    foreach (@bad_id_list) {
        (my $h = $_->id) =~ s/[\w]//g;
        printf "%s : %s\n", $_->id, $h;
        #print $_->id, "\n";
    }
}


#- Check for non-standard sequence characters --------------------------
foreach (@seq_list) {
    push(@bad_seq_list, $_) if ($_->seq =~ /[^AGCTUMRWSYKBDHV\-NX?]/gi);
}

#output
if (@bad_seq_list) {
    
    printf "\nWARNING! Non-IUPAC characters detected in %i sequences:\n", scalar @bad_seq_list;
    foreach (@bad_seq_list) {
        (my $s = $_->seq) =~ s/[AGCTUMRWSYKBDHV\-NX?]//gi;
        printf "%s : %s\n", $_->id, $s;
    }
}

if (not @bad_id_list and not @bad_seq_list and not @duplicate_id_list) {
    print "Sequence data OK\n";
}
else { exit 2 }


exit;
