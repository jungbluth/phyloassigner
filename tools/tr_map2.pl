#!/usr/bin/perl -w

######################## Copyright and License #########################
#
# tr_map2.pl - Output the '.mapping' file for PhyloAssigner
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

use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;

my $raxml_file = shift or die "RAxML/pplacer tree?\n";
my $arb_file = shift or die "arb tree (with node names, after clean-up)?\n";
my $type = shift or die "Does the first tree come from raxml or pplacer?\n";
my $outgroup = shift || 'Archaea';

die "Type of tree file should be either 'raxml' or 'pplacer'\n" unless ($type eq 'raxml' or $type eq 'pplacer');

my $r_in = Bio::TreeIO->new(-file=>$raxml_file, -format =>'newick');
my $rax_tree = $r_in->next_tree;

if ($type eq 'raxml') {
   # pre-process RAxML tree: move raxml node IDs to a tag from node names
   for my $rn ($rax_tree->get_nodes) {
      next if $rn->internal_id == $rax_tree->get_root_node->internal_id;
      my $id_in = $rn->id;
      warn "Need RAxML IDs of nodes in the form _I1_, _I2_ at the end of node IDs.\nWe have: $id_in\n" unless $id_in =~ /_(I\d+)_$/;
      my $rax_id = $1;
      $rn->add_tag_value("raxid", $rax_id);
      if ($rn->is_Leaf) {
         $id_in =~ s/_I\d+_$//;
         $rn->id($id_in);
      }
   }
} else  {
   # pre-process pplacer tree
   my $max_id = 0;
   for my $rn ($rax_tree->get_nodes) {
      my $id_in = $rn->id;
      my $pp_id;
      if ($rn->is_Leaf) {
         warn "Need pplacer IDs of nodes in the form 1\@node_a, 2\@node_b.\nWe have: $id_in\n" unless $id_in =~ /^(\d+)\@/;
         $pp_id = $1;
         $id_in =~ s/^\d+\@//;
         $rn->id($id_in);
      } else {
         $pp_id = $id_in;
      }
      $rn->add_tag_value("raxid", $pp_id);
      $max_id = $pp_id if $pp_id > $max_id;
   }
   $rax_tree->get_root_node->set_tag_value("raxid", $max_id + 1);
}

my $a_in = Bio::TreeIO->new(-file=>$arb_file, -format=>'newick');
my $arb_tree = $a_in->next_tree;

#   trial to re-root RAxML tree to fit arb tree - can also do this eg in dendroscope
my $a_archaea = $arb_tree->find_node(-id=>$outgroup);
my @archaea;
#~ =begin comment1
if ($a_archaea->is_Leaf) {
   push @archaea, $a_archaea->id;
} else {
   @archaea = map {$_->id} grep {$_->is_Leaf} $a_archaea->get_all_Descendents;
}
my $tmp_root;
for my $nnn (grep {$_->is_Leaf} $rax_tree->get_nodes){
   unless (grep {$_ eq $nnn->id} @archaea) {
      $tmp_root = $nnn;
      last;
   }
}
$rax_tree->reroot($tmp_root);
my @r_arch_nodes=();
for my $a (@archaea) {
   die "Did not find node with ID $a in raxml tree\n" unless $rax_tree->find_node(-id=>$a);
   push @r_arch_nodes, $rax_tree->find_node(-id=>$a);
}

my $ar_lca;
if (scalar @r_arch_nodes > 1) {
   $ar_lca = $rax_tree->get_lca(@r_arch_nodes);
} else {
   $ar_lca = $r_arch_nodes[0];
}
$rax_tree->reroot($ar_lca);
my $tt_out = Bio::TreeIO->new(-file=>">$raxml_file.rooted", -format=>'newick');
$tt_out->write_tree($rax_tree);


my %had;
for my $n ($rax_tree->get_nodes) {
   #next if $n->internal_id == $rax_tree->get_root_node->internal_id;
   my @clade=();
   my $node_id = $n->get_tag_values("raxid") ;  
   next unless defined $node_id;
   print STDERR "Node ID $node_id not unique in pplacer tree!\n" if $had{$node_id};
   $had{$node_id}++;
   if ($n->is_Leaf) {
      @clade = ($n->id);
   } else {
      map {push @clade, $_->id if $_->is_Leaf} $n->get_all_Descendents;
   }
   map {s/\_I\d+\_$//} @clade; 
   # find the same node in arb tree:
   my @arb_clade = ();
   for my $cn (@clade) {
      die "Did not find node $cn in arb tree.\n" unless my $cnn =  $arb_tree->find_node(-id=>$cn);
      push @arb_clade, $cnn;
   }
   die "Mapping nodes between trees failed\n" unless scalar @arb_clade == scalar @clade;
   my $arb_node;
   if (scalar @clade > 1) {
      $arb_node = $arb_tree->get_lca(@arb_clade);
   } else {
      $arb_node = $arb_clade[0];
   }
   my $taxonomy = "";
   for my $ln ($arb_tree->get_lineage_nodes($arb_node)) {
      $taxonomy = $taxonomy . ";" . $ln->id if $ln->id;
   }
   if (! $arb_node->is_Leaf and defined $arb_node->id) {
      $taxonomy = $taxonomy . ";" . $arb_node->id;
   }
   $taxonomy =~ s/^\;//;
   print $node_id . "\t";
   print $taxonomy ; #"\t";
   #print join ", ", @clade;
   if (scalar @clade == 1) {
      print "\t" . $clade[0];
   }
   print "\n";
}
#~ =end comment1
#~ =cut
