#!/usr/bin/perl

######################## Copyright and License #########################
#
# setupdb - Setup of a PhyloAssigner reference database
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

use Cwd;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::TreeIO;
use File::Basename;
use File::Copy;
#use File::Spec
use File::Spec::Functions qw( catfile );
use File::Path qw( make_path rmtree );
use Getopt::Long;


#- Variables -----------------------------------------------------------

my $program_info = "
        -------------------------------------------------------------
        setupdb - Setup of a PhyloAssigner reference database
        This program is part of PhyloAssigner v0.2.8.16
        Copyright (C) 2012 Fabian Kilpert, Bank Beszteri, Adam Monier
        This program comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions. Use option '--license' for details.
        -------------------------------------------------------------
";

print $program_info , "\n";


# usage
my $usage = "Usage:  setupdb.pl [options] ref_tree ref_alignment outgroup

    ref_tree            : reference NEWICK tree with clade names
    ref_alignment       : alignment FASTA file (e.g. from silva-arb)
    outgroup            : outgroup (name must exist in ref_tree!!!)
    
    [options]
    -o dir              : output directory
    -v                  : verbose mode
    --nocleanup         : disable cleanup (default: cleanup enabled)
    --nopack            : disable packing of .phyloassignerdb to .tar.gz
    --start position    : first position (bp) of required alignment slice
    --end position      : last position (bp) of required alignment slice
    --hmmerdir dir      : dir conatining HMMER3 binaries (DEFAULT: ./)
    --pplacerdir dir    : dir conatining pplacer 1.1 binaries (DEFAULT: ./)
    --phymldir dir      : dir containing PhyML3 binaries (DEFAULT: ./)
    
    -h                  : show this help
    --help              : show extended help (README)
    --license           : show GNU GPL v3 (COPYING)

";


## main ##
my $phyloassigner_dir = dirname(File::Spec->rel2abs($0));       #main dir of PhyloAssigner
my $tools_dir = catfile($phyloassigner_dir, "tools");           #tools dir of PhyloAssigner
my $cwd = getcwd;

## options ##
my $slice_start;        # first position of required alignment slice
my $slice_end;          # last position of required alignment slice
my $verbose;            # verbose mode
#my $dbname;             # user specified output name
my $hmmer_dir;          # Hmmer directory
my $phyml_dir;          # PhyML directory
my $pplacer_dir;        # Pplacer directory
my $no_cleanup;         # disable cleanup
my $no_pack;            # DONT pack output folder to archive (*.tar.gz)
my $help;               # flag to show help (usage)
my $help_ext;           # flag to show extended help (README)
my $license;            # flag to show license

## arguments ##
my $in_tree;            # input reference tree (with clade labels)
my $in_alignment;       # input reference alignment (fasta)
my $in_outgroup;        # user defined outgroup, e.g. "Alveolata/Dinophyceae"

## other ##
my $output_dir;         # output directory
my $shell_cmd;          # hold command line arguments
my @args;               # collect command line arguments 
my @garbage;            # collect files and folders for eventual clean up
my $archive;            # Final packed archive .phyloassigerdb.tar.gz
my $basename;           # basename of input file
my $base;               # basename of input file without extension


#- Get options and arguments -------------------------------------------

## options ##
GetOptions ('start=i' => \$slice_start,
            'end=i' => \$slice_end,
            'output_dir=s' => \$output_dir,
            'verbose' => \$verbose,
            'hmmerdir=s' => \$hmmer_dir,
            'pplacerdir=s' => \$pplacer_dir,
            'phymldir=s' => \$phyml_dir,
            'nocleanup' => \$no_cleanup,
            'nopack' => \$no_pack,
            'h' => \$help,
            'help' => \$help_ext,
            'license' => \$license)
            or die $usage;

# help
die $usage if $help;
if ( $help_ext ) {
    $shell_cmd = "more ".catfile($phyloassigner_dir, "README");
    system("$shell_cmd") == 0 or die "Unable to show README!\n";
    exit;
}

#license
if ( $license ) {
    $shell_cmd = "more ".catfile($phyloassigner_dir, "COPYING");
    system("$shell_cmd") == 0 or die "Unable to show COPYING!\n";
    exit;
}

die $usage if not defined @ARGV;


# pplacer1.1 dir
if ( $pplacer_dir ) {
    # if --pplacerdir
    -e catfile($pplacer_dir, "pplacer") or die "ERROR! Check --pplacerdir option! pplacer NOT found in:\n$pplacer_dir\n\n", $usage;
} else {
        $pplacer_dir = catfile( $phyloassigner_dir, "binaries") if -e catfile( $phyloassigner_dir, "binaries/pplacer");
        $pplacer_dir = get_PATH('pplacer') if get_PATH('pplacer');
        $pplacer_dir = $cwd if -e catfile($cwd, "pplacer");
}
die "ERROR! pplacer NOT found! Set path e.g. by using --pplacerdir option!!!\n\n" if not $pplacer_dir;


# HMMER3 dir
if ( $hmmer_dir ) {
    # if --hmmerdir
    -e catfile($hmmer_dir, "hmmbuild") or die "ERROR! Check --hmmerdir option! hmmbuild NOT found in:\n$hmmer_dir\n\n", $usage;
} else {
        $hmmer_dir = catfile( $phyloassigner_dir, "binaries") if -e catfile( $phyloassigner_dir, "binaries/hmmbuild");
        $hmmer_dir = get_PATH('hmmbuild') if get_PATH('hmmbuild');
        $hmmer_dir = $cwd if -e catfile($cwd, "hmmbuild");
}
die "ERROR! hmmer NOT found! Set path e.g. by using --hmmerdir option!!!\n\n" if not $hmmer_dir;


# phyml dir
if ( $phyml_dir ) {
    # if --phymldir
    -e catfile($phyml_dir, "phyml") or die "ERROR! Check --phymldir option! phyml NOT found in:\n$phyml_dir\n\n", $usage;
} else {
        $phyml_dir = catfile( $phyloassigner_dir, "binaries") if -e catfile( $phyloassigner_dir, "binaries/phyml");
        $phyml_dir = get_PATH('phyml') if get_PATH('phyml');
        $phyml_dir = $cwd if -e catfile($cwd, "phyml");
}
die "ERROR! phyml NOT found! Set path e.g. by using --phymldir option!!!\n\n" if not $phyml_dir;


## arguments ##
print "WARNING: Too many arguments in command line!\n" if scalar @ARGV > 3;

#reference tree
if ( defined $ARGV[0] ) {
    $in_tree = $ARGV[0];
    -e $in_tree or die "ERROR: Reference tree file NOT found!!!\n$in_tree\n\n";
} else { die "ERROR: Reference tree file NOT specified!!!\n\n"; }

#reference alignment
if ( defined $ARGV[1] ) {
    $in_alignment = $ARGV[1];
    -e $in_alignment or die "ERROR: Reference alignment file NOT found!!!\n$in_alignment\n\n"; 
    #$basename = basename($in_alignment);                        #basename of reference alignment
    ( $base = basename($in_alignment) ) =~ s/\.[^.]+$//;  #basename without extension
} else { die "ERROR: Reference alignment file NOT specified!!!\n\n"; }

#outgroup
if ( defined $ARGV[2] ) {
    $in_outgroup = $ARGV[2];
} else { die "ERROR: Outgroup NOT specified!!!\n\n"; }


## other ##

# output dir (database in dir)
$output_dir = $base.'.phyloassignerdb' if not defined $output_dir;                         #default
$output_dir = $output_dir.'.phyloassignerdb' if ( $output_dir !~ m/\.phyloassignerdb$/ );   #add extension if necessary


## Path check ##
-d $phyloassigner_dir or die "ERROR: PhyloAssigner folder NOT found!!!\n$phyloassigner_dir\n\n";
-d $tools_dir or die "ERROR: Tools folder NOT found!!!\n$tools_dir\n\n";


## set variables ##
my $ali_len = &alignment_length($in_alignment);    #length of alignment sequences
my $ref_aln = catfile($output_dir, $base.".aln");
my $ref_aln_sto = catfile($output_dir, $base.".sto");
my $ref_aln_phy = catfile ($output_dir, $base.".phy");
my $ref_hmm = catfile($output_dir, $base.".hmm3");
( my $ref_tree = basename($in_tree) ) =~ s/\.[^.]+$//;
$ref_tree = catfile($output_dir, $ref_tree.".tre");
my $ref_stats = catfile($output_dir, $base.".phy_phyml_stats.txt");

my $phyloassigner_testdb_dir = catfile($output_dir, "TEST".'.phyloassignerdb');
my $test_seq = catfile($phyloassigner_testdb_dir, "TEST".".fas");
my $phyloassigner_output_dir = catfile($output_dir, "TEST" . '.output');
my $pplacer_map = catfile($output_dir, $base.".mapping");
my $jplace = catfile($phyloassigner_output_dir, "TEST.fas.aln.jplace");
my $jplace_tree = catfile($output_dir, basename($jplace).'.tree');


## some output ##
print "S E T T I N G S :\n\n";

printf "%-15s : %s\n" , "Ref tree", $in_tree;
printf "%-15s : %s\n" , "Ref alignment", $in_alignment;
printf "%-15s : %s\n" , "Ref outgroup", $in_outgroup;
printf "%-15s : %s\n" , "HMMER 3 dir", $hmmer_dir;
printf "%-15s : %s\n" , "PhyML dir", $phyml_dir;
printf "%-15s : %s\n" , "pplacer1.1 dir", $pplacer_dir;
printf "%-15s : %s\n" , "Main dir", $phyloassigner_dir if $verbose;
printf "%-15s : %s\n" , "Tools dir", $tools_dir if $verbose;
if (not $no_cleanup) { printf "%-15s : %s\n" , "Cleanup", "Yes"; } else { printf "%-15s : %s\n" , "Cleanup", "No"; }
if ($verbose) { printf "%-15s : %s\n" , "Verbose", "Yes"; } else { printf "%-15s : %s\n" , "Verbose", "No"; }
printf "%-15s : %s\n" , "Slice start", $slice_start if $slice_start;
printf "%-15s : %s\n" , "Slice end", $slice_end if $slice_start;
printf "%-15s : %s\n" , "Output dir", $output_dir;


#=======================================================================


#- Lenght of sequences in alignment for slice --------------------------
print "-" x 80, "\n";
print "Determing alignment length...\n";

print "Alignment length [bp]: $ali_len\n";
 
# if no slice specified
$slice_start = 1 if not defined $slice_start; 
$slice_end = $ali_len if not defined $slice_end; 

# check for invalid values specified by user
die "Error: invalid slice value (--from): $slice_start\n\n", $usage if ($slice_start < 1 or $slice_start >= $ali_len);
die "Error: invalid slice value (--to): $slice_end\n\n", $usage if ($slice_end > $ali_len or $slice_end <= $slice_start);

print "Alignment slice [bp]: $slice_start-$slice_end\n";



#- Remove output directory if already existing -------------------------
if (-d $output_dir) {
    print "-" x 80, "\n";
    print "Output directory $output_dir already exists, trying to continue run where we left
off.\n";
    print "If you would like to start a new run from scratch, remove or rename $output_dir
or specify an alternative output directory!\n";
    #print "New output will OVERWRITE data in output directory:\n$output_dir\n";
    #rmtree($output_dir) or die "Unable to remove existing output directory!\n";
} else {
    make_path($output_dir, {mode => 0744,}) or die "Unable to create output directory!\n";
}


#- Cleanup of ref alignment --------------------------------------------
unless (-e $ref_aln) {
   print "-" x 80, "\n";
   print "Processing reference alignment (e.g. from silva-arb)...\n";

   $shell_cmd = sprintf "perl %s %s %i %i > %s",
                            catfile($tools_dir, "rm_gap.pl"),
                            $in_alignment,
                            $slice_start,
                            $slice_end,
                            $ref_aln;
   print $shell_cmd, "\n" if $verbose;
   print "This could take a while. Please be patient...\n";
   system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
} else {
   print "Reference alignment has already been processed...\n";
}
-e $ref_aln or die "File does NOT exist!\n$ref_aln\n";

#- Convert alignment from FASTA to STOCKHOLM format --------------------
unless (-e $ref_aln_sto) {
   print "-" x 80, "\n";
   print "Converting alignment format from FASTA to STOCKHOLM format...\n";

   #(my $basename = basename($ref_aln)) =~ s/\.[^.]+$//, "\n";
   #~

   $shell_cmd = sprintf "perl %s %s %s %s",
                (catfile($tools_dir, 'fasta2stockholm.pl'),
                $ref_aln,
                '>', $ref_aln_sto);

   print $shell_cmd, "\n" if $verbose;
   system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
} else {
   print "Reference alignment has already been converted to Stockholm format...\n";
}
push(@garbage, $ref_aln_sto);
-e $ref_aln_sto or die "File does NOT exist!\n$ref_aln\n";
print $ref_aln, "\n";

#- Build hmm -----------------------------------------------------------
unless (-e $ref_hmm) {
   print "-" x 80, "\n";
   print "Building hmm (HMMER3)...\n";

   @args = (catfile($hmmer_dir, 'hmmbuild'),
            '-n', $base,
            $ref_hmm,
            $ref_aln_sto);
   print "@args\n" if $verbose;
   system(@args) == 0 or die "Unable to run:\n@args\n";
} else {
   print "Reference HMM is already there...\n";
}
-e $ref_hmm or die "File does NOT exist!\n$ref_hmm\n";

#- Convert alignment from fasta to phy ---------------------------------
unless (-e $ref_aln_phy) {
   print "-" x 80, "\n";
   print "Converting alignment format from FASTA to PHYLIP...\n";

   @args = ('perl', catfile($tools_dir, 'fas2phy.pl'),
        $ref_aln,
        $ref_aln_phy);
   print "@args\n" if $verbose;
   system(@args) == 0 or die "Unable to run:\n@args\n";
} else {
   print "Reference alignment has already been converted to PHYLIP format...\n";
}
push(@garbage, $ref_aln_phy);



#- Save cleaned up copy of reference tree ------------------------------
unless (-e $ref_tree){
   print "-" x 80, "\n";
   print "Saving cleaned up copy of reference tree...\n";

   #remove the comment at the beginning of the ref_tree file and save a copy in the outdir folder
   &process_tree_file($in_tree, $ref_tree) or die "Unable to generate ref tree file!\n";
   print $ref_tree, "\n" if $verbose;
} else {
   print "Reference tree has already been cleaned up...\n";
}
push(@garbage, $ref_tree);

# check for consistent sequence ids between alignment and tree:
print "-" x 80, "\n";
print "Checking consistent sequence IDs between alignment and tree...\n";
#my $problems = &check_aln_tree ($ref_aln_phy, $ref_tree) ;
my $problems = &check_aln_tree ($ref_aln, $ref_tree) ;
die "There are $problems inconsistent sequence IDs between reference alignment and tree.\n"
if ($problems > 0);

print "Sequence ID consistency check OK.\n";

#- Run PhyML -----------------------------------------------------------
unless (-e catfile($output_dir, $base.".phy_phyml_tree.txt") && -e $ref_stats) {
   print "-" x 80, "\n";
   print "Running PhyML...\n";

   $shell_cmd = sprintf "%s %s %s %s %s %s %s %s %s %s %s %s %s %s",
                'echo Y |',
                catfile($phyml_dir, 'phyml'),
                '-i', catfile($output_dir, $base.".phy"),
                '-m', 'GTR',
                '-f', 'e', 
                '-a', 'e',
                '-u', $ref_tree,
                '-o', 'lr';
   print $shell_cmd, "\n" if $verbose;
   system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
}
$ref_tree = catfile($output_dir, $base.".phy_phyml_tree.txt");
-e $ref_tree or die "File does NOT exist!\n$ref_tree\n";
-e $ref_stats or die "File does NOT exist!\n$ref_stats\n";



#- Make test-phyloassigerdb --------------------------------------------
print "-" x 80, "\n";
print "Preparing PhyloAssigner test run...\n";

# remove output dir if existing
if (-d $phyloassigner_testdb_dir) {
    print "-" x 80, "\n";
    print "New output will OVERWRITE data in output directory:\n$phyloassigner_testdb_dir\n";
    rmtree($phyloassigner_testdb_dir) or die "Unable to remove existing output directory!\n";
}
make_path($phyloassigner_testdb_dir, {mode => 0744,}) or die "Unable to create output directory!\n";
push(@garbage, $phyloassigner_testdb_dir);


#copy files 
copy ( $ref_aln, catfile($phyloassigner_testdb_dir, 'REF.aln') ) or die "Unable to copy alignment to subfolder!";
copy ( $ref_hmm, catfile($phyloassigner_testdb_dir, 'REF.hmm3') ) or die "Unable to copy profile HMM3 to subfolder!";
copy ( $ref_stats, catfile($phyloassigner_testdb_dir, 'REF.phy_phyml_stats.txt') ) or die "Unable to copy statistics file to subfolder!";
copy ( $ref_tree, catfile($phyloassigner_testdb_dir, 'REF.phy_phyml_tree.txt') ) or die "Unable to copy tree file to subfolder!";



#- Make test sequence from reference alignment -------------------------
print "-" x 80, "\n";
print "Making test sequence from reference alignment...\n";
# test file with only a single sequence (first sequence from the reference dataset)
print $ref_aln, "\n";
print $test_seq, "\n";
&make_testfile($ref_aln, $test_seq) or die "Unable to generate test sequence!\n";



#- Run PhyloAssigner for test -----------------------------------------------
print "=" x 80, "\n";
print "Running PhyloAssigner with test sequence...\n";

@args = ('perl', catfile($phyloassigner_dir, 'phyloassigner.pl'),
        '-o', $phyloassigner_output_dir,
        '--nopost',
        '--nocleanup',
        '--hmmerdir', $hmmer_dir,
        '--pplacerdir', $pplacer_dir,
        $phyloassigner_testdb_dir,
        $test_seq);
print "@args\n" if $verbose;
system(@args) == 0 or die "Unable to run:\n@args\n";
push(@garbage, $phyloassigner_output_dir);



#- Parse pplacer output (JSON) and modify tree -------------------------
print "=" x 80, "\n";
print "Parsing pplacer output (JSON) and modifying tree...\n";

-e $jplace or die "File does NOT exist!\n$jplace\n";

$shell_cmd = sprintf "%s %s %s > %s",
                'perl', catfile($tools_dir, 'extract_num_ref_tree.pl'),
                $jplace,
                $jplace_tree;
print $shell_cmd, "\n" if $verbose;
system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
print $jplace_tree, "\n" if $verbose;
push(@garbage, $jplace_tree);



#- Modify labeled reference tree ---------------------------------------
# *.tree -> *.tree.unodes
print "-" x 80, "\n";
print "Processing labeled reference tree...\n";

# make a copy of the labeled reference tree in working folder
copy ($in_tree, $output_dir) or die "Unable to copy labeled reference tree to subfolder!";

@args = ('perl', catfile($tools_dir, 'clean_arb_tree.pl'),
            catfile($output_dir, basename($in_tree)));
print "@args\n" if $verbose;
system(@args) == 0 or die "Unable to run:\n@args\n";
push(@garbage, catfile($output_dir, basename($in_tree.'.unodes')));



#- Mapping -------------------------------------------------------------
print "-" x 80, "\n";
print "Mapping...\n";

$shell_cmd = sprintf "perl %s %s %s %s \"%s\" > %s",
                        catfile($tools_dir, 'tr_map2.pl'),
                        $jplace_tree,
                        catfile($output_dir, basename($in_tree).'.unodes'),
                        "pplacer",
                        $in_outgroup,
                        catfile($output_dir, $base.'.mapping');
print $shell_cmd, "\n" if $verbose;
system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
-e $pplacer_map or die "File does NOT exist!\n$pplacer_map\n";
push(@garbage, $jplace_tree.'.rooted');



#- Clean up of temporary files and folders -----------------------------
if (not $no_cleanup) {
    print "-" x 80, "\n";
    print "Cleaning up files and folders no longer needed...\n";
    print "Now deleting:\n" if $verbose;
    
    foreach (@garbage) {
        rmtree ($_) or die "Unable to remove: $_";
        print $_, "\n" if $verbose;
    }
}



#- Packing of .phyloassignerdb to .tar.gz ------------------------------
if (not $no_pack) {
    print "-" x 80, "\n";
    print "Finally, packing of phyloassingerdb to archive (tar.gz)...\n";
    
    chdir dirname($output_dir);
    
    $archive = basename($output_dir).".tar.gz";
    
    $shell_cmd = sprintf "tar cfvz %s %s",
                            $archive,
                            basename($output_dir);
    print $shell_cmd, "\n" if $verbose;
    system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
    -e basename($output_dir).".tar.gz" or die "File does NOT exist!\n$archive\n";
    
    chdir $cwd;
    
    rmtree ($output_dir) or die "Unable to remove: $output_dir";
}



#- End of script. ------------------------------------------------------
print "-" x 80, "\n";
print "setupdb.pl finished successfully.\n\n";
if ( $archive and -e catfile(dirname($output_dir),$archive) ) {
    print "PhyloAssigner database (archive) stored to:\n".catfile(dirname($output_dir),$archive)."\n\n";
} else {
    print "PhyloAssigner database stored to:\n$output_dir\n\n";
}

#=begin comment
#=end comment
#=cut

#= Subroutines =========================================================

# Remove arb comment from tree file
sub process_tree_file {
    my $infile = shift @_;
    my $outfile = shift @_;

    my $comment ;
    my $tree_string = "";

    #read infile
    open (my $in_fh, "<", $infile) or die "Could not open $infile\n$!\n";
    while (my $line = <$in_fh>) {
        $comment = 1 if $line =~ /\[/;
        $tree_string .= $line unless $comment;
        $comment = 0 if $line =~ /\]$/;
    }
    close $in_fh;

    #write outfile
    open (my $out_fh, ">", $outfile) or die "Could not open $outfile for writing\n$!\n";
    print {$out_fh} $tree_string;
    close $out_fh;
    return 1;
} 


# Make test sequence from reference alignment
sub make_testfile {
    my $infile = shift @_;
    my $outfile = shift @_;
    
    my $seq_in = Bio::SeqIO->new('-file' => "<$infile", '-format' => 'fasta');
    my $seq = $seq_in->next_seq();
    $seq->id("TEST_".$seq->id);     #make an addition to sequence name
    
    my $seq_out = Bio::SeqIO->new('-file' => ">$outfile", '-format' => 'fasta');
    $seq_out->write_seq($seq);
    
    return 1;
}


# Get alignment length
sub alignment_length {
    my $file = shift @_;
    
    my $seq_in = Bio::SeqIO->new(-file => "<$file",
                            -format => 'fasta');
    
    my $len;    # sequence length
    my $is_first = 1;   # flag
    
    while (my $seq = $seq_in->next_seq()) {
        if ($is_first) {
            $len = $seq->length;
            $is_first = 0;
        }
        else {
            if ($len != $seq->length) {
                die "Error: Sequences in alignment file do NOT have same length!";
            }
        }
    }
    return $len;
}


#get eviromental path for a given file (if exists)
sub get_PATH {
    my $query_file = shift @_;
    foreach ( split(':', $ENV{'PATH'}) ) {
        return $_ if -e  catfile($_,  $query_file);
    }
}


sub check_aln_tree {
   my $aln_file = shift;
   my $tree_file = shift;
   
   my %aln_ids;
   my $aln_in = Bio::AlignIO->new (-file => $aln_file, -format => 'fasta') or die $!;
   my $aln = $aln_in->next_aln();
   for my $seq ($aln->each_seq) {
      $aln_ids{$seq->id} ++;
      #print "id: ", $seq->id, "\n";
   }

   my $tree_in = Bio::TreeIO->new (-file => $tree_file, -format => 'newick') or die $!;
   my $tree = $tree_in->next_tree();
   my $problems = 0;
   for my $l ($tree->get_leaf_nodes) {
      my $leaf = $l->id;
      #print "leaf: $leaf\n";
      unless ($aln_ids{$leaf}) {
         $problems ++;
         print STDERR "Sequence $leaf present in tree, but not in alignment.\n";
      } else {
         $aln_ids{$leaf} = 0;
      }
   }
   for my $id (keys %aln_ids) {
       #print "id: $id\n"
      if ($aln_ids{$id}){
         $problems ++;
         print STDERR "Sequence $id present in alignment, but not in tree.\n";
      }
   }
   return $problems;
}


