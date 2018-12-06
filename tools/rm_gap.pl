#!/usr/bin/perl -w

######################## Copyright and License #########################
#
# rm_gap.pl - Remove gaps
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
use Bio::AlignIO;

my $inf = shift or die "Infile?\n";
my $from = shift or die "First aln position to include (1 - aln length)?\n";
my $to  = shift or die "Last aln position to include (1 - aln length)?\n";


my %seq;
my $reader = Bio::AlignIO->new( -file => $inf, -format => 'fasta');
my $aln = $reader->next_aln;

die "First aln position should be between 1 and " . $aln->length . "\n" if $from  < 1 or $from > $aln->length;
die "Last aln position should be between 1 and " . $aln->length . "\n" if $to  < 1 or $to > $aln->length;

my @all_gaps = 0..($to-$from);

for my $s ($aln->each_seq) {
   my @seq_chars = (split //, $s->seq)[($from-1)..($to-1)];
   $seq{$s->id} = \@seq_chars;
   my $cc = $#all_gaps;
   while ($cc >= 0) {
      unless ($seq_chars[$all_gaps[$cc]] eq '-' or $seq_chars[$all_gaps[$cc]] eq '.') {
         splice @all_gaps, $cc, 1;
      }
      $cc --;
   }
}

my %gaps;
map {$gaps{$_} =1} @all_gaps;
for my $id (keys %seq) {
   print ">$id\n";
   for my $ii (0..($to-$from)) {
      unless ($gaps{$ii}) {
         if ($seq{$id}->[$ii] eq '.') {
	    print '-';
	 } elsif ($seq{$id}->[$ii] eq 'U' or $seq{$id}->[$ii] eq 'u') {
	    print 'T';
	 } else {
	    print $seq{$id}->[$ii];
	 }
      }
   }
   print "\n";
}
