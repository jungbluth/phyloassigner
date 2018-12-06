#!/usr/bin/perl

######################## Copyright and License #########################
#
# assign_labels.pl - Labels from LCA identification
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


use warnings;
use strict;

use Cwd;
use Getopt::Long;

use JSON;       #JSON - JSON (JavaScript Object Notation) encoder/decoder
use Bio::TreeIO;


my $usage = "assign_labels.pl [options] mapping_file placement_file
    
    Parse pplacer output (JSON), identify the last common ancestor (LCA)
    of best placements and output to table.
    
    mapping_file    : Text file mapping of node IDs to taxon names
    placement_file  : Output file from pplacer 1.1 (JSON file)
    
    Options:
    --leafnodes     : Include mapping to leaf nodes
    --pp            : Use posterior probabilities instead of likelihood 
                      weights for LCA identification
    --lca           : Cutoff cumulative likelihood weight / posterior 
                      probability value. Top scoring nodes are 
                      collected for each sequence classified until this 
                      cumulative cutoff is reached, and the last common 
                      ancestor of these nodes is written to the output. 
                      (Default: 0.9)
    -v              : verbose mode

            ";


#Option variables
my $pp;               # Bayesian posterior probability (PP) mode
my $lca = 0.9;        # LCA cumulative threshold
my $verbose;          # verbose mode
my $leafnodes;        # allow mapping to leaf nodes


#- Get options and arguments -------------------------------------------
#Options
GetOptions ('pp' => \$pp,
            'lca=f' => \$lca,
            'verbose' => \$verbose,
            'leafnodes' => \$leafnodes),
            or die $usage;

#Arguments
my $mapping_file = Cwd::abs_path(shift @ARGV) or die $usage;
-e $mapping_file or die "ERROR: Mapping file NOT found!!!\n\n", $usage;  

my $placement_file = Cwd::abs_path(shift @ARGV) or die $usage;
-e $placement_file or die "ERROR: Placement file NOT found!!!\n\n", $usage;


#some output on variables
print "Mapping file: $mapping_file\n" if $verbose;
print "Placement file: $placement_file\n" if $verbose;
if ($verbose) {
    if ($leafnodes) {print "Leaf node mapping: Yes\n"} else {print "Leaf node mapping: No\n"}
}
print "Verbose mode!\n" if $verbose;
print "LCA cutoff: $lca\n" if $verbose;
if ($pp) {print "Mode: Bayesian posterior probability (PP)\n" if $verbose}
else {print "Mode: Maximum Likelihood (ML)\n" if $verbose}



#- read mapping file ---------------------------------------------------

my %mapping;    #Mapping of tree node number to taxon string

open (my $m_fh, "<", $mapping_file) or die $!;
while (<$m_fh>) {
    if (length $_ > 0 and not $_ =~ /^\#/) {
        $_ =~ s/\s+$//;         #delete one or more linefeeds
        my  @columns = split /\t/, $_;
        $mapping{$columns[0]} = [ \$columns[1], \$columns[2] ];
    }
}
close $m_fh;

#~ # some output
#~ while ( my($key, $value) = each %mapping ) {
    #~ if ( ${$value->[1]} ) { print $key, " ", ${$value->[0]}, "\t", ${$value->[1]}, "\n"; }
    #~ else { print $key, " ", ${$value->[0]}, "\n"; }
#~ }



#- read placement file (JSON) ------------------------------------------

my $json_text;
my $version;
my $tree_string;
my $metadata_invocation;
my @fields;
my @results;

my %seq_dic;
my @seq_list;


open (my $p_fh, "<", $placement_file ) or die $!;
while (<$p_fh>) {
    $json_text .= $_;
}
close $p_fh;

my $href_json = decode_json( $json_text );



#- go through all main JSON elements -----------------------------------
while ( my($key, $value) = each %{$href_json} ) {

    #get version
    if ($key eq 'version') {
        $version = $href_json->{version};
        #print $version, "\n";
    }
    
    #get tree
    if ($key eq 'tree') {
        $tree_string = $href_json->{tree};
        #print $tree_string, "\n";
    }
    
    #get metadata
    if ($key eq 'metadata') {
        while ( my($m_key, $m_value) = each %{$href_json->{metadata}} ) {
            if ($m_key eq 'invocation') {
                $metadata_invocation = $m_value;
                #print $metadata_invocation, "\n";
            }
        }
    }

    #get fields
    if ($key eq 'fields') {
        @fields = @{$href_json->{fields}};
        #print "@fields\n";
    }
    
    #get placements
    if ($key eq 'placements') {
        foreach my $href_placements (@{$href_json->{placements}}){
            
            my $name;
            
            #names
            my $aref_n = $href_placements->{n};
            foreach (@{$aref_n}) {
                $name = $_;
                push (@seq_list, $name);
            
            
                #results
                my $aref_p = $href_placements->{p};
                
                my @seq_results_list;
                
                foreach ( @{$aref_p} ) {
                    my %seq_results_dic;
                    
                    @results = @{$_};
            
                    # e.g.: edge_num, likelihood, like_weight_ratio, distal_length, pendant_length, (post_prob, marginal_prob)
                    for(0..(scalar @fields)-1) { 
                        $seq_results_dic{$fields[$_]} = $results[$_];
                    }
                    push (@seq_results_list, \%seq_results_dic);
                }
                $seq_dic{$name} = \@seq_results_list;
            }
        }
    }
}
print "JSON file successfully loaded!\n\n" if $verbose;



#- make tree usable using Bioperl --------------------------------------

# (UGLY HACK!) In this special tree BioPerl is unable to parse node numbering from newick 
# By replacing {} by [] bioperl can recognize the node numbering by the bootsrap (UGLY HACK!)
$tree_string =~ s/{/[/g; 
$tree_string =~ s/}/]/g;

# open tree string as a file (make usable file handle)
open (my $fake_tree_fh, "+<", \$tree_string) or die $!;   # +< : read and write

#make newick tree a bioperl object
my $bio_tree = new Bio::TreeIO(-fh => $fake_tree_fh, -format => 'newick');
my $tree = $bio_tree->next_tree;

#simplify nodes
foreach ($tree->get_nodes) {
    #perldoc Bio::Tree::Node
    
    #print "ID: ", $_->id, "\t bootsrap: ", $_->bootstrap, "\n";
    $_->id($_->bootstrap);      #In this special tree BioPerl stores the IDs in the 'bootstrap' variable!!! (UGLY HACK!)

}


#- Generate output -----------------------------------------------------

#Description line
printf "#%s, %s, %s, %s, %s, %s, %s\n", "sequence ID", "best node", "best taxon", "best weight", "LCA node", "LCA taxon", "LCA weight";

# output JSON data
foreach my $s(@seq_list) {
    
    #output variables
    my $seq_id = $s;
    my $best_node;
    my $best_taxon;
    my $best_weight = 0;
    my $lca_node;
    my $lca_taxon;
    my $lca_weight = 0;
 
    my @placements = @{$seq_dic{$s}};
    my @node_list = ();

    foreach (@placements) {
        my $weight_ratio;
        
        #get relative weight ratio (like_weight_ratio OR post_prob)
        if ( $pp ) { 
            $weight_ratio = ${$_}{'post_prob'} or die "ERROR: No Bayesian posterior probability (PP) mode supported by current file!!!\n\n";
        } 
        else { 
            $weight_ratio = ${$_}{'like_weight_ratio'};
        }
        
        my $edge_num = ${$_}{'edge_num'};
        
        #get best placement
        if ( $weight_ratio > $best_weight ) {
            $best_weight = $weight_ratio;
            $best_node = $edge_num;
            
            $best_taxon = ${$mapping{$best_node}->[0]} || "NA";
            
            #leafnodes
            $best_taxon .= ";" . ${$mapping{$best_node}->[1]} if ${$mapping{$best_node}->[1]} and $leafnodes and $best_taxon ne "NA"; 
        }
        
        #LCA cutoff
        if ( $lca_weight < $lca ) {
            $lca_weight += $weight_ratio;
            push @node_list, $edge_num;
        }
    }
    
    #If there are less then two nodes left...
    if ( scalar @node_list < 2 ) {
        $lca_node = $best_node;
        $lca_taxon = $best_taxon;
    } else {
        #print "### node_list: @node_list\n";
        $lca_node = $tree->get_lca( map {$tree->find_node($_)} @node_list )->id;
        $lca_node =~ s/\[(\d*)\]/$1/;    #this is an ugly hack to correct an error in BioPerl, which gives incorrect result for base node!!!
        
        #print "### LCA: ", $lca_node, "\n";
        
        #$lca_taxon = $mapping{$lca_node};
        
        $lca_taxon = ${$mapping{$lca_node}->[0]} || "NA";
        
        #leafnodes
        $lca_taxon .= ";" . ${$mapping{$lca_node}->[1]} if ${$mapping{$lca_node}->[1]} and $leafnodes and $lca_taxon ne "NA"; 
    }
    
    #print "seq_id", $seq_id, "\n";
    $seq_id = "-" if not defined $seq_id;
    $best_node = "-" if not defined $best_node;
    $best_taxon = "-" if not defined $best_taxon;
    $best_weight = "-" if not defined $best_weight;
    $lca_node = "-" if not defined $lca_node;
    $lca_taxon = "-" if not defined $lca_taxon;
    $lca_weight = "-" if not defined $lca_weight;
    
    #final output line
    print   $seq_id, "\t",
            $best_node, "\t",
            $best_taxon, "\t",
            $best_weight, "\t",
            $lca_node, "\t",
            $lca_taxon, "\t",
            $lca_weight, "\n";
}


