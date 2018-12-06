#!/usr/bin/perl

######################## Copyright and License #########################
#
# extract_num_ref_tree.pl - read tree from pplacer outfile, write NEWICK
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


my $usage = "extract_num_ref_tree.pl [options] placement_file
    
    Parse pplacer output (JSON) and modify tree.
    
    placement_file  : Output file from pplacer 1.1 (JSON file)
    
    Options:
    -v              : verbose mode

            ";


#Option variables
my $verbose;          # verbose mode


#- Get options and arguments -------------------------------------------
#Options
GetOptions ('verbose' => \$verbose),
            or die $usage;

#Arguments
my $placement_file = Cwd::abs_path(shift @ARGV) or die $usage;
-e $placement_file or die "ERROR: Placement file NOT found!!!\n\n", $usage;

#some output on variables
print "Placement file: $placement_file\n" if $verbose;
print "Verbose mode!\n" if $verbose;



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

#~ #print $tree_string, "\n";
#~ 
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
    
    # root node
    if ( $_ == $tree->get_root_node ) {
        my $root_id = $_->id;
        $root_id =~ s/\[(\d*)\]/$1/;
        $_->id($root_id);
        $_->branch_length(0);
    }
    # other nodes
    else {
        
        $_->id($_->bootstrap."@".$_->id);      #In this special tree BioPerl stores the IDs in the 'bootstrap' variable!!! (UGLY HACK!)
    }
    #print "ID: ", $_->id, "\t bootsrap: ", $_->bootstrap, "\n";
}

my $tree_as_string = $tree->as_text('newick');
print $tree_as_string, "\n";

