#!/usr/bin/perl -w

######################## Copyright and License #########################
#
# clean_arb_tree.pl - read NEWICK tree from ARB, write cleaned tree
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

## script to clean up a tree file exported from arb (with node names)
## also takes care of non-unique node names

use Bio::TreeIO;

my $arb_file = shift or die "arb tree file?\n";
my $tmp_file = $arb_file . ".tmp";

## clean infile of (most) spaces
open (my $in, "<", $arb_file) or die $!;
open (my $tout, ">", $tmp_file) or die $!;
while (<$in>) {
   s/\s*\)/)/g;
   s/\s*\(/(/g;
   s/\n//;
   print $tout $_;
}
close $in;
close $tout;

my $a_in = Bio::TreeIO->new(-file=>$tmp_file, -format =>'newick');
my $arb_tree = $a_in->next_tree;

map {my $n=$_->id; $n =~ s/\s//g; $_->id($n)} grep {$_->is_Leaf} $arb_tree->get_nodes;

my %had;
for my $a_node (grep {! $_->is_Leaf} $arb_tree->get_nodes) {
   next unless $a_node->id;
   next if $a_node->id =~ /^\s*$/;
   (my $label = $a_node->id) =~ s/\s+$//;
   $label =~ s/^'//;
   $label =~ s/'$//;
   if ($had{$label}) {
      #if ($label =~ /^\'.*\'$/) {
      #   $label =~ s/'$//;
      #   $label = $label . "_" . $had{$label}++ . "'";
      #} else {
         $label = $label . "_" . $had{$label}++;
      #}
   } else {
      $had{$label} ++;
   }
   $a_node->id($label);
   my @taxa = grep {$_->is_Leaf} $a_node->get_all_Descendents;
   
   #print "$label\n";
}
my $a_out = Bio::TreeIO->new(-file=>">$arb_file.unodes", -format=>'newick');
$a_out->write_tree($arb_tree);
unlink $tmp_file;
