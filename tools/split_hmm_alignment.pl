#!/usr/bin/perl -w

######################## Copyright and License #########################
#
# split_hmm_alignment.pl - Split in query and reference; output FASTA
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

my $inf = shift or die "Infile?\n";
my $q_file = shift or die "Query file?\n";

my %q;
if ($q_file ne 'NONE') {
   open (my $ifh, "<", $q_file) or die $!;
   while (<$ifh>) {
      chomp;
      next if /^\s?$/;
      $q{$_} ++;
   }
} else {
   (my $qname = basename($inf)) =~ s/\..*//;
   $q{$qname} = 1;
}
open (my $qfh, ">", $inf . ".q.fasta") or die $!;
open (my $rfh, ">", $inf . ".ref.fasta") or die $!;
$/ = ">";
open (my $in, "<", $inf) or die $!;
my $disc = <$in>;
my %had; ## to make sure we don't feed duplicates to pplacer, as it dies if we do
while (my $seq = <$in>){
   chomp $seq;
   my ($id, $s) = split /\n/, $seq;
   next if $had{$id};
   $had{$id} ++;
   my $out = $q{$id} ? $qfh : $rfh;
   print $out ">$seq";
}

close $in;
close $qfh;
close $rfh;
