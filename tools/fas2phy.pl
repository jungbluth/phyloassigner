#!/usr/bin/perl -w

######################## Copyright and License #########################
#
# fas2phy - Convert FASTA to PHYLIP format
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
use File::Basename;
use Bio::AlignIO;

my $id_max_len = 120;   #maximum length of sequence IDs (it's PHYLIP!)

my $inf = shift or die "Infile?\n";
my $outf = shift or die "Outfile?\n";
(my $q = basename($inf)) =~ s/\..*//;
my $i = Bio::AlignIO->new(-file => $inf, -format=>'fasta');
my $o = Bio::AlignIO->new(-file => ">$outf", -format=>'phylip',
-idlength=>$id_max_len);
my $a = $i->next_aln;
$a->set_displayname_flat();
$o->write_aln($a);
