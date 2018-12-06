#!/usr/bin/perl

######################## Copyright and License #########################
#
# label_count.pl - Count labels in '.tab' file, make summary
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



my $usage = "label_count.pl tab_file
    Collect total number of taxonomic labels from a tab-file.

    tab_file  : PhyloAssigner tab-file
    
    Options:
    -o file        : output file
    -v             : verbose mode
    
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
$outfile = $infile.".counts" if not $outfile;
$outfile = File::Spec->rel2abs($outfile);
print "Outfile: $outfile\n" if $verbose;
-e $infile or die "Could NOT open file! $!\n";


#- Read file -----------------------------------------------------------
open (my $in_fh, '<', $infile) or die "Unable to open file: \"$infile\" $!\n";
my @lines = <$in_fh>;
close $in_fh;



#- Process lines -------------------------------------------------------
my $i = 0;      # line counter
my %lca_label_dic;  # counting all found lca labels in a dic
my %best_label_dic;  # counting all found best labels in a dic

foreach (@lines) {
    $_ =~ s/\s+$//;         #delete one or more linefeeds
    
    if (length $_ > 0 and not $_ =~ /^\#/) {
        #print $_, "\n";
        ++$i;
        
        my @column_list = split(/\t/, $_);
        #print scalar @column_list, "\n";
        
        if (scalar @column_list!=7) {
            die "ERROR! Some lines do not have 7 columns!\n$_\n\n";
        } 
        else {
            my $label;      # a single label
            
            #- lca column --------------------------------------------------
            $label = &parse_label($column_list[5]);  #get the last label in column
            #print $i, " ", $label, "\n";
            
            # get absolute numbers by filling a dictionary/hash
            if (exists $lca_label_dic{$label}) {$lca_label_dic{$label} += 1}
            else {$lca_label_dic{$label} = 1}
            
            
            #- lca column --------------------------------------------------
            $label = &parse_label($column_list[2]);  #get the last label in column
            #print $i, " ", $label, "\n";
            
            # get absolute numbers by filling a dictionary/hash
            if (exists $best_label_dic{$label}) {$best_label_dic{$label} += 1}
            else {$best_label_dic{$label} = 1}
        }
    }
}

print "\n" if $verbose;;

open (my $fh_out, '>', $outfile) or die "Couldn't open file \"$outfile\": $!";

my $output_line;



#- best placements -----------------------------------------------------
my $best_total = 0;
foreach (values %best_label_dic) {
    $best_total += $_;
}

$output_line = "Best placements (".$best_total." total):\n";
foreach my $k (sort (keys(%best_label_dic))) {
    $output_line .= $k. "\t".$best_label_dic{$k}."\n";
}
print $output_line if $verbose;
print {$fh_out} $output_line;

print "\n" if $verbose;
print {$fh_out} "\n";


#- LCA placements ------------------------------------------------------
my $lca_total = 0;
foreach (values %lca_label_dic) {
    $lca_total += $_;
}

$output_line = "LCA placements (".$lca_total." total):\n";
foreach my $k (sort (keys(%lca_label_dic))) {
    $output_line .= $k. "\t".$lca_label_dic{$k}."\n";
}
print $output_line if $verbose;
print {$fh_out} $output_line;

print "\n" if $verbose;

close $fh_out;


#= Subroutines =========================================================

# parse the label column
sub parse_label {
    my $cur_column = shift @_;
    #print $cur_column, "\n";
    
    my @label_list = split(/;/, $cur_column);
    my $label;
    
    if (scalar @label_list == 1) { $label = shift @label_list; }
    else { $label = pop @label_list; }
    
    return $label;
}

