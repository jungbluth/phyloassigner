#!/usr/bin/perl

######################## Copyright and License #########################
#
# PhyloAssigner - Phylogenetic placement and taxonomic label assignment
# Copyright 2012 Fabian Kilpert, Bank Beszteri, Adam Monier
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
use File::Basename;
use File::Spec::Functions qw( catfile catdir );
use File::Path qw( make_path rmtree );
use File::Copy;
use Getopt::Long;


#- Variables -----------------------------------------------------------

my $program_info = "
        -------------------------------------------------------------
        PhyloAssigner - Pipeline for phylogenetic placement and 
          taxonomic label assignment, v0.2.7.17
        Copyright (C) 2012 Fabian Kilpert, Bank Beszteri, Adam Monier
        This program comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions. Use option '--license' for details.
        -------------------------------------------------------------
";

print $program_info , "\n";


# usage
my $usage = "Usage:  phyloassigner.pl [options] database fasta

    database            : .phyloassignerdb database (dir OR .tar.gz)
    fasta               : input FASTA file with query sequences
    
    [options]
    -o dir              : output directory
    -t threads          : number of parallel processing threads (default: 1)
    -p                  : Bayesian posterior probability (PP) mode 
                          (default: maximum likelihood (ML) mode)
    -v                  : verbose mode
    --leafnodes         : include mapping to leaf nodes
    --nopost            : disable post-processing
    --nocleanup         : disable cleanup (default: cleanup enabled)
    --nocheck           : disable input file consistency check 
                          (default: check enabled)
    --jplacefile        : output placefile
    --tabfile           : output placefile
    --lca threshold     : LCA cumulative threshold (default: 0.9)
    --nocleanup         : disable cleanup (default: cleanup enabled)
    --alginer mothur    : set aligner (default: 'hmmer')
    --hmmerdir dir      : dir conatining HMMER3 binaries (DEFAULT: ./)
    --pplacerdir dir    : dir conatining pplacer 1.1 binaries (DEFAULT: ./)
    
    -h                  : show this help
    --help              : show extended help (README)
    --license           : show GNU GPL v3 (COPYING)

";


## main ##
my $phyloassigner_dir = dirname(File::Spec->rel2abs($0));       #main dir of PhyloAssigner
my $tools_dir = catfile($phyloassigner_dir, "tools");           #tools dir of PhyloAssigner
my $cwd = getcwd;

## options ##
my $threads = 1;        # number of parallel output files for threading
my $pp;                 # Bayesian posterior probability (PP) mode
my $verbose;            # verbose mode
my $leafnodes;          # allow mapping to leaf nodes
my $no_post;            # disable postprocessing
my $no_check;           # disable input file consistency check
my $jplaceoutfile;      # output placefile
my $taboutfile;         # output placefile
my $lca = 0.9;          # LCA cumulative threshold
my $no_cleanup;         # disable cleanup
my $hmmer_dir;          # Hmmer directory
my $phyml_dir;          # PhyML directory
my $pplacer_dir;        # Pplacer directory
my $aligner;            # mothur OR hmmalign (default)
my $mothur_dir;         # mothur dir
my $help;               # flag to show help (usage)
my $help_ext;           # flag to show extended help (README)
my $license;            # flag to show license

## arguments ##
my $in_file;            # input file
my $in_db;              # input phyloassignerdb (dir or tar.gz)
my $db;                 # input phyloassignerdb
my $db_base;            # database base name without extension

## other ##
my $output_dir;         # output directory
my $shell_cmd;          # command line for starting external programs
my @args;               # collect command line arguments
my @garbage;            # collect files and folders for eventual clean up
my $basename;           # basename of input file
my $base;               # basename of input file without extension
my @childs;             #collect child processes for forking
my @temp_dirs;          #collect all temporary dirs
my @files_to_process;   #collect files for phylogenetic placement
my @place_file_list;    #collect place-files
my @tab_list;           #collect tab-files
my $placefile;          #place-file
my $tabfile;            #tab-file

my $ref_aln;            # ref alignment file
my $ref_hmm;            # ref profile HMM (HMMER3) file
my $ref_tree;           # ref NEWICK tree to place reads on
my $ref_stats;          # ref stats file for ref tree
my $ref_mapping;        # ref mapping table
my $query_id_list;      # query id list file


#- Get options and arguments -------------------------------------------

#Options
GetOptions ('t=i' => \$threads,
            'output_dir=s' => \$output_dir,
            'pp' => \$pp,
            'verbose' => \$verbose,
            'nopost' => \$no_post,
            'jplacefile=s' => \$jplaceoutfile,
            'tabfile=s' => \$taboutfile,
            'nocleanup' => \$no_cleanup,
            'nocheck' => \$no_check,
            'lca=f' => \$lca,
            'aligner=s' => \$aligner,
            'hmmerdir=s' => \$hmmer_dir,
            'pplacerdir=s' => \$pplacer_dir,
            'mothurdir=s' => \$mothur_dir,
            'leafnodes' => \$leafnodes,
            'h' => \$help,
            'help' => \$help_ext,
            'license' => \$license),
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

## arguments ##
die $usage if not defined @ARGV;
print "WARNING: Too many arguments in command line!\n" if scalar @ARGV > 2;

# mothur dir
if ($aligner and $aligner eq 'mothur') {
    if ( $mothur_dir ) {
        # if --hmmdir
        -e catfile($mothur_dir, "mothur") or die "ERROR! Check --mothurdir option! mothur NOT found in:\n$mothur_dir\n\n", $usage;
    } else {
        $mothur_dir = catfile( $phyloassigner_dir, "binaries") if -e catfile( $phyloassigner_dir, "binaries/mothur");
        $mothur_dir = get_PATH('mothur') if get_PATH('mothur');
        $mothur_dir = $cwd if -e catfile($cwd, "mothur");
    }
    die "ERROR! mothur NOT found! Set path e.g. by using --mothurdir option!!!\n\n" if not $mothur_dir;
}

else {
   # HMMER3 dir (default)
    if ( $hmmer_dir ) {
        # if --hmmdir
        -e catfile($hmmer_dir, "hmmbuild") or die "ERROR! Check --hmmerdir option! hmmbuild NOT found in:\n$hmmer_dir\n\n", $usage;
    } else {
        $hmmer_dir = catfile( $phyloassigner_dir, "binaries") if -e catfile( $phyloassigner_dir, "binaries/hmmbuild");
        $hmmer_dir = get_PATH('hmmbuild') if get_PATH('hmmbuild');
        $hmmer_dir = $cwd if -e catfile($cwd, "hmmbuild");
    }
    die "ERROR! hmmbuild NOT found! Set path e.g. by using --hmmerdir option!!!\n\n" if not $hmmer_dir;
}

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




# options check
$threads >= 1 or die "Error: Threads : $threads -> Specify one or more threads!";

#input database
if ( defined $ARGV[0] ) {
    $in_db = $ARGV[0];
    -e $in_db or die "ERROR: Reference database (.phyloassignerdb) NOT found!!!\n$in_db\n\n"; 
    if ( not -d $in_db ) { 
        if ( $in_db =~ m/\.tar\.gz/) {
            
            ## If packed ##
            
            print "Unpacking phyloassignerdb to current working directory...\n";
            if ( $verbose ) { $shell_cmd = "tar xfvz ".$in_db; }
            else { $shell_cmd = "tar xfz ".$in_db; }
            print $shell_cmd, "\n" if $verbose;
            system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
        
            print "-" x 80, "\n";
            ( $db_base = basename($in_db) ) =~ s/.phyloassignerdb.tar.gz//;
            $db = basename(catdir($cwd, $db_base.".phyloassignerdb"));
            push(@garbage, $db);
        }
        else { die "ERROR: Unknown file type!\n$in_db\n\n" } 
    } else { 
        
        ## If dir ##
        
        ( $db_base = basename($in_db) ) =~ s/.phyloassignerdb//;
        $db = catdir(dirname($in_db), $db_base.".phyloassignerdb");
    }
    
    ## set paths to reference db files ##
    # Reference alignment
    foreach ( $db_base, 'REF' ) { $ref_aln = catfile( $db, $_ . '.aln' ) if -e catfile( $db, $_ . '.aln' ) }
    # Reference HMM profile
    foreach ( $db_base, 'REF' ) { $ref_hmm = catfile( $db, $_ . '.hmm3' ) if -e catfile( $db, $_ . '.hmm3' ) }
    # Reference tree
    foreach ( $db_base, 'REF' ) { $ref_tree = catfile( $db, $_ . '.phy_phyml_tree.txt' ) if -e catfile( $db, $_ . '.phy_phyml_tree.txt' ) }
    # Reference tree stats
    foreach ( $db_base, 'REF' ) { $ref_stats = catfile( $db, $_ . '.phy_phyml_stats.txt' ) if -e catfile( $db, $_ . '.phy_phyml_stats.txt' ) }
    #Reference mapping
    foreach ( $db_base, 'REF' ) { $ref_mapping = catfile( $db, $_ . '.mapping' ) if -e catfile( $db, $_ . '.mapping' ) }
    
} else { die "ERROR: Reference database file (.phyloassignerdb) NOT specified!!!\n\n"; }

#reference alignment
if ( defined $ARGV[1] ) {
    $in_file = $ARGV[1];
    -e $in_file or die "ERROR: Input FASTA file NOT found!!!\n$in_file\n\n"; 
    $basename = basename($in_file);                     #basename of reference alignment
    ( $base = basename($in_file) ) =~ s/\.[^.]+$//;     #basename without extension
} else { die "ERROR: Input FASTA file NOT specified!!!\n\n"; }


#input query file
if ( defined $ARGV[1] ) {
    $in_file = $ARGV[1];
    -e $in_file or die "ERROR: Input query file NOT found!!!\n$in_file\n\n"; 
    #$basename = basename($in_file);                        #basename of reference alignment
    ( $base = basename($in_file) ) =~ s/\.[^.]+$//;  #basename without extension
} else { die "ERROR: Input query file NOT specified!!!\n\n"; }


## other ##

# set output dir if not defined by -o option above 
$output_dir = $in_file.'.output' if not defined $output_dir;    #default


## Path check ##



## set variables ##

#set default alinger
if (not defined $aligner) { $aligner = 'hmmalign'; }
else { $aligner = 'hmmalign' if not (grep (/^$aligner$/, ('hmmalign','mothur'))); }

#basename of input file
$basename = basename($in_file); 

#hmmalign alignment (hmm output)
my $aligner_out = catfile($output_dir, $basename.".aln.hmmer3");

# converted to fasta
my $aln = catfile($output_dir, $basename.".aln");


## some output ##
print "S E T T I N G S :\n\n";
printf "%-15s : %s\n" , "Database", $db;
printf "  %-13s : %s\n" , "Alignment", $ref_aln if $verbose;
printf "  %-13s : %s\n" , "Profile HMM", $ref_hmm if $verbose;
printf "  %-13s : %s\n" , "Tree", $ref_tree if $verbose;
printf "  %-13s : %s\n" , "Stats", $ref_stats if $verbose;
if ( $verbose ) { 
    if ( $ref_mapping ) { printf "  %-13s : %s\n" , "Mapping", $ref_mapping }
    else { printf "  %-13s : %s\n" , "Mapping", '-'}
    }
printf "%-15s : %s\n" , "Query file", $in_file;
printf "%-15s : %s\n" , "Output dir", $output_dir;
if ($aligner eq 'mothur') {
    printf "%-15s : %s\n" , "Aligner", "Mothur";
    printf "%-15s : %s\n" , "mothur dir", $mothur_dir;
} elsif ($aligner eq 'hmmalign'){
    printf "%-15s : %s\n" , "Aligner", "hmmalign (HMMER3)";
    printf "%-15s : %s\n" , "HMMER 3 dir", $hmmer_dir;
}
printf "%-15s : %s\n" , "pplacer1.1 dir", $pplacer_dir;
if ($pp) { printf "%-15s : %s\n" , "pplacer mode", "Bayesian posterior probability (PP)"; }
else { printf "%-15s : %s\n" , "pplacer mode", "Maximum Likelihood (ML)"; }
printf "%-15s : %s\n" , "Threads", $threads;
printf "%-15s : %s\n" , "LCA", $lca;
if ($leafnodes) { printf "%-15s : %s\n" , "Leaf nodes", "Yes"; } else { printf "%-15s : %s\n" , "Leaf nodes", "No"; }
if (not $no_check) { printf "%-15s : %s\n" , "Infile check", "Yes"; } else { printf "%-15s : %s\n" , "Infile check", "No"; }
if (not $no_post) { printf "%-15s : %s\n" , "Postprocessing", "Yes"; } else { printf "%-15s : %s\n" , "Postprocessing", "No"; }
if (not $no_cleanup) { printf "%-15s : %s\n" , "Cleanup", "Yes"; } else { printf "%-15s : %s\n" , "Cleanup", "No"; }
if ($verbose) { printf "%-15s : %s\n" , "Verbose", "Yes"; } else { printf "%-15s : %s\n" , "Verbose", "No"; }
if ($jplaceoutfile) { printf "%-15s : %s\n" , "jplace outfile", $jplaceoutfile; }
if ($taboutfile) { printf "%-15s : %s\n" , "tab outfile", $taboutfile; }


#=======================================================================


#- Remove output directory if already existing -------------------------
if (-d $output_dir) {
    print "-" x 80, "\n";
    print "New output will OVERWRITE data in output directory:\n$output_dir\n";
    rmtree($output_dir) or die "Unable to remove existing output directory!\n";
}
make_path($output_dir, {mode => 0744,}) or die "Unable to create output directory!\n";



#- make query id list --------------------------------------------------
print "-" x 80, "\n";
print "Making query ID list...\n";

#file containing list of query ids
$query_id_list = catfile($output_dir, $basename.".qids");
@args = ('perl', catfile($tools_dir, 'query_ids.pl'),
            "-o", $query_id_list, 
            $in_file);
splice(@args, scalar @args-1, 0, "-v") if $verbose;     # set -v flag
system(@args) == 0 or die "\nQuery sequence file did NOT pass the consistency check!!!\n$in_file\n\n\n";
push(@garbage, $query_id_list);



#~ #- convert DOS input file to Unix -----------------------------------
#~ print "-" x 80, "\n";
#~ print "Converting from DOS to UNIX...\n";
#~ 
#~ @args = ("perl", catfile($tools_dir, "dos_to_unix.pl"),
            #~ '-o', $in_file,
            #~ $in_file);
#~ splice(@args, scalar @args-1, 0, "-v") if $verbose;     # set -v flag
#~ system(@args) == 0 or die "Unable to run:\n@args\n";



#- Check sequence input file (query file) for consistency --------------
if (not $no_check) {
    print "-" x 80, "\n";
    print "Checking query sequences consistency...\n";
    @args = ('perl', catfile($tools_dir, 'check_fasta.pl'),
                $in_file);
    system(@args) == 0 or die "\nQuery sequence file did NOT pass the consistency check!!!\n$in_file\n\n\n";
}



#= Aligner =============================================================

if ($aligner eq 'mothur') {

    #- Run mothur ------------------------------------------------------
    print "-" x 80, "\n";
    print "Align query sequences to reference profile using mothur...\n";

    copy($ref_aln, $output_dir) or die "Can't copy ref alignment to output dir! $!\n";     # copy reference alignment to output dir
    copy($in_file, $output_dir) or die "Can't copy query file to output dir! $!\n";        # copy query file to output dir
    
    chdir $output_dir or die "Can't change working dir to: $!\n";    
    # run mothur
    $shell_cmd = sprintf "%s '#align.seqs(candidate=%s, template=%s, processors=%i)'", 
                            catfile($mothur_dir, "mothur"), basename($in_file), basename($ref_aln), $threads;
    print "$shell_cmd\n" if $verbose;
    system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
    chdir $cwd or die "Can't change back to working dir: $!\n";

    my $mothur_aln = catfile( $output_dir, $basename.'.aln.q.fasta' );
    move( catfile( $output_dir, $base.'.align' ), $mothur_aln );
    
    move( catfile( $output_dir, basename($ref_aln)), catfile( $output_dir, $basename.'.aln.ref.fasta' ) );
    print "$mothur_aln\n" if $verbose;
    process_query_aln($mothur_aln, $mothur_aln);            # change '.' by '-'

    push(@garbage, catfile( $output_dir, basename($in_file) ) );
    push(@garbage, catfile( $output_dir, $base.'.align.report' ) );
    push(@garbage, catfile( $output_dir, $db_base.'.8mer' ) );
    push(@garbage, catfile( $output_dir, $basename.'.aln.q.fasta' ) );
    push(@garbage, catfile( $output_dir, $basename.'.aln.ref.fasta' ) );
    
    #- Combine query and reference to single file (for pplacer) --------
    #~ if ($threads == 1) {
        print "-" x 80, "\n";
        print "Combining query and reference to single file (for pplacer)...\n";
        
        $shell_cmd = sprintf "cat %s %s > %s", 
                        $ref_aln,
                        $mothur_aln,
                        $aln.".fasta";
        print "$shell_cmd\n" if $verbose;
        system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
    #~ }
    copy($aln.".fasta", $aln);
    push(@garbage, $aln.".fasta");
}


#= OR hmmer aligner ====================================================

else {

    #- Run hmmalign ----------------------------------------------------
    print "-" x 80, "\n";
    print "Align query sequences to reference profile using hmmalign (HMMER3)...\n";
    
    print "hmmer_dir: ", $hmmer_dir, "\n";
    print "ref_aln: ", $ref_aln, "\n";
    print "aligner_out: ", $aligner_out, "\n";
    print "ref_hmm: ", $ref_hmm, "\n";
    print "in_file: ", $in_file, "\n";


    #- Convert alignment from FASTA to STOCKHOLM format --------------------
    print "-" x 80, "\n";
    print "Converting alignment format from FASTA to STOCKHOLM format...\n";

    #(my $basename = basename($ref_aln)) =~ s/\.[^.]+$//, "\n"; 
    my $ref_aln_new = $aln.".sto";

    $shell_cmd = sprintf "perl %s %s %s %s",
                    (File::Spec->catfile($tools_dir, 'fasta2stockholm.pl'),
                    $ref_aln,
                    '>', $ref_aln_new);
                    
    print $shell_cmd, "\n" if $verbose;
    system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
    push(@garbage, $ref_aln_new);
    
    -e $ref_aln or die "File does NOT exist!\n$ref_aln\n";
    
    
    # hmmer3 -----------------------------------------------------------
    @args = (catfile($hmmer_dir, "hmmalign"),
                "--trim",
                "--informat", "FASTA",
                "--outformat", "afa",
                "--allcol",
                "--mapali", $ref_aln_new,
                #"--cpu", $threads,
                "-o", $aligner_out,
                $ref_hmm,
                $in_file);
    print "@args\n" if $verbose;
    system(@args) == 0 or die "Unable to run:\n@args\n";
    
    
    #- Parsing hmmer3 outfile (AFA); replace . by -, make UPPERCASE ----
    print "-" x 80, "\n";
    print "Parsing hmmer3 outfile (AFA)...\n";
    
    $shell_cmd = sprintf "perl %s %s > %s", 
                        catfile($tools_dir, "parse_hmmer3_outfile.pl"), 
                        $aligner_out,
                        $aln;
    #print "$shell_cmd\n" if $verbose;
    system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";

    print "$aln\n" if -e $aln and $verbose;
    push(@garbage, $aligner_out);
    


    # Copy and rename to a file with .fasta suffix (recommended by pplacer)
    print "-" x 80, "\n";
    print "Renaming alignment for use in pplacer...\n";
    @args = ("cp", $aln, $aln.".fasta");
    
    print "@args\n" if $verbose;
    system(@args) == 0 or die "Unable to run:\n@args\n";
    
    push(@garbage, $aln.".fasta");


    
    # Split alignment (FASTA) for multithreading 
    if ($threads > 1) {
        print "-" x 80, "\n";
        print "Splitting alignment into reference and query alignment...\n";

        @args = ("perl", catfile($tools_dir, "split_hmm_alignment.pl"),
                    $aln,
                    $query_id_list);
        
        #print "@args\n" if $verbose;
        system(@args) == 0 or die "Unable to run:\n@args\n";

        print $aln.".ref.fasta", "\n" if -e $aln.".ref.fasta" and $verbose;
        print $aln.".q.fasta", "\n" if -e $aln.".q.fasta" and $verbose;
        push(@garbage, $aln.".ref.fasta");
        push(@garbage, $aln.".q.fasta");
        
    }
    
}
#= alignment finished ==================================================



#- Split query fasta file into several files ---------------------------
if ($threads > 1) {
    print "-" x 80, "\n";
    print "Split alignment file for parallel processing...\n";
    
    my $splitfasta_path = catfile($tools_dir, "split_fasta.pl");
    
    @args = ("perl", $splitfasta_path,
                "-n", $threads,
                $aln.".q.fasta");
    splice(@args, scalar @args-1, 0, "-v") if $verbose;     # set -v flag
    print "@args\n" if $verbose;
    system(@args) == 0 or die "Unable to run:\n@args\n";
    
    #make list of partial query files (files_to_process)
    for(0..$threads-1) { 
        my $q_part = sprintf ("${aln}.q.part%.2i.fasta", $_);
        if (-e $q_part) { 
            push(@files_to_process, $q_part);
            push(@garbage, $q_part);
            #print $q_part, "\n";
        }
        else { die "ERROR: Partial q files does NOT exist!\n$q_part\n\n"; }
    }
}



#- Start pplacer -------------------------------------------------------
print "-" x 80, "\n";

# Single threaded
if ($threads == 1) {
    print "Starting pplacer...\n";
    
    #print -e $aligner_out.".q.fasta", "\n";
    my @args = (catfile($pplacer_dir, "pplacer"), 
                    "--out-dir", $output_dir,
                    #"-r", $aln.".ref.fasta",
                    "-t", $ref_tree,
                    "-s", $ref_stats,
                    $aln.".fasta");
    splice(@args, scalar @args-1, 0, "-p") if $pp;     # set -p flag
    print "@args\n" if $verbose;
    system(@args) == 0 or die "Unable to run: @args\n";
    
    #$placefile = $aligner_out.".jplace" if -e $aligner_out.".jplace";
    $placefile = $aln.".jplace" if -e $aln.".jplace";
}
# Multiple threads
elsif ($threads > 1) {
    for(0..$threads-1) {
        my $pid = fork();
        if ($pid) {
            # parent
            push(@childs, $pid);
        } 
        elsif ($pid == 0) {
            # child
            print "Starting pplacer #$_...\n";
            
            my @args = (catfile($pplacer_dir, "pplacer"), 
                        "--out-dir", $output_dir,
                        "-r", $aln.".ref.fasta",
                        "-t", $ref_tree,
                        "-s", $ref_stats,
                        $files_to_process[$_]);
            splice(@args, scalar @args-1, 0, "-p") if $pp;     # set -p flag
            print "@args\n" if $verbose;
            system(@args) == 0 or die "Unable to run: @args\n";
            
            exit(0);
        } 
        else {
            die "Couldn’t fork: $!\n";
        }
    }
    foreach (@childs) { waitpid($_, 0); }
    
    #make a list of output files (jplace files)
    foreach(@files_to_process) { 
        $_ =~ s/\.\w+$/\.jplace/; 
        push(@garbage, $_);
    }
}



#- Merge placement (jplace) files... ---------------------
if ($threads > 1) {
    print "-" x 80, "\n";
    print "Merging placement files (jplace)...\n";
        
    #make file name for combined placement file
    $placefile = $aln.".jplace";
    
    @args = (catfile($pplacer_dir, 'guppy'),
                    'merge',
                    '-o', $placefile, 
                    @files_to_process);
    print "@args\n" if $verbose;
    system(@args) == 0 or die "Unable to run: @args\n";
    print $placefile, "\n" if $verbose;
}


#= if mapping file is available do mapping related processing ==========

if ($ref_mapping) {

#=======================================================================

    #- LCA identification --------------------------------------------------
    print "-" x 80, "\n";

    # Single threaded
    if ($threads == 1) {
        print "Starting LCA identification...\n";
        
        $shell_cmd = sprintf "perl %s %s %s %s %s %s %s > %s", 
                        catfile($tools_dir, "assign_labels.pl"), 
                        '--lca', $lca,
                        ($leafnodes ? '--leafnodes' : ''),
                        ($pp ? '--pp' : ''),
                        $ref_mapping,
                        $aln.".jplace",
                        $aln.".jplace.tab";
        print "$shell_cmd\n" if $verbose;
        system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
        
        $tabfile = $aln.".jplace.tab";
    }
    # Multiple threads
    elsif ($threads > 1) {
        printf "Starting LCA identification (using %s threads)...\n", $threads;
         
        for(0..$threads-1) {
            (my $tabfile = $files_to_process[$_]) .= ".tab";
            
            my $pid = fork();
            if ($pid) {
                # parent
                push(@childs, $pid);
            } 
            elsif ($pid == 0) {
                # child
                print "Starting label assignment.pl #$_...\n";
                
                my $shell_cmd = sprintf "perl %s %s %s %s %s %s %s > %s",
                                            catfile($tools_dir, "assign_labels.pl"), 
                                            '--lca', $lca,
                                            ($leafnodes ? '--leafnodes' : ''),
                                            ($pp ? '--pp' : ''),
                                            $ref_mapping,
                                            $files_to_process[$_],
                                            $tabfile;
                print $shell_cmd, "\n" if $verbose;
                system("$shell_cmd") == 0 or die "Unable to run: $shell_cmd\n";
                
                exit(0);
            } 
            else {
                die "Couldn’t fork: $!\n";
            }
            
            push(@tab_list, $tabfile);
            push(@garbage, $tabfile);
            #foreach ( @files_to_process ) { push(@garbage, $_); }
        }
        
        foreach (@childs) {
            waitpid($_, 0);
        }
    }


    #- Merging tab files -----------------------------------------------
    if ($threads > 1) {
        print "-" x 80, "\n";
        print "Merging tab-files...\n";

        $tabfile = $aln.".jplace.tab";
        my $description_line;

        open(my $fh1, '>', $tabfile) or die "Unable to open file \"$tabfile\": $!";

        foreach my $t (@tab_list) {
            open(my $fh2, '<', $t) or die "Unable to open file \"$t\": $!";
            while (my $line=<$fh2>) {
                
                #take care that the description line is printed only the first time
                if ( $line =~ /^\#/) {
                    print ( {$fh1} $line ) if not defined $description_line;
                    $description_line = $line;
                } else {
                    print {$fh1} $line;
                }
            }
            close $fh2;
        }
        close $fh1;
    }
    
    print $aln.".jplace.tab", "\n" if -e $aln.".jplace.tab" and $verbose;
    


#=======================================================================

}   #end of mapping related processing
else {
    print "-" x 80, "\n";
    print "No mapping data available! Skipping LCA identification!!\n";
}
#=======================================================================

#- Clean up of temporary files and folders -----------------------------
if (not $no_cleanup) {
    print "-" x 80, "\n";
    print "Cleaning up files and folders no longer needed...\n";
    print "Now deleting:\n" if $verbose;
    
    foreach (@garbage) {
        rmtree ($_) or print "Unable to remove: $_\n";
        print $_, "\n" if $verbose;
    }
}



#=======================================================================
#= POST-PROCESSING =====================================================
#=======================================================================

if (not $no_post) {

    print "=" x 80, "\n";
    print "Doing some post-processing...\n";

    #- Count taxon labels from tab-file ------------------------------------
    print "-" x 80, "\n";
    print "Count taxon labels...\n";

    @args = ("perl", catfile($tools_dir, "label_count.pl"), 
                $aln.".jplace.tab");
    splice(@args, scalar @args-1, 0, "-v") if $verbose;
    print "@args\n" if $verbose;
    system(@args) == 0 or die "Unable to run: @args\n";


    #pplacer's guppy tools...

    #~ visualization 
        #~ fat                     makes trees with edges fattened in proportion to the number of reads
        #~ heat                    maps an an arbitrary vector of the correct length to the tree
        #~ ref_tree                writes a taxonomically annotated reference tree and an induced taxonomic tree
        #~ sing                    makes one tree for each query sequence, showing uncertainty
        #~ tog                     makes a tree with each of the reads represented as a pendant edge
    #~ 
      #~ statistical comparison
        #~ bary                    draws the barycenter of a placement collection on the reference tree
        #~ bootviz                 makes a phyloXML tree showing the bootstrap values
        #~ edpl                    calculates the EDPL uncertainty values for a collection of pqueries
        #~ kr                      calculates the Kantorovich-Rubinstein distance and corresponding p-values
        #~ kr_heat                 makes a heat tree
        #~ pca                     performs edge principal components
        #~ splitify                writes out differences of masses for the splits of the tree
        #~ squash                  performs squash clustering
    #~ 
      #~ classification
        #~ classify                outputs classification information in a tabular or SQLite format
    #~ 
      #~ utilities
        #~ check_refpkg            check a reference package
        #~ demulti                 splits apart placements with multiplicity, undoing a round procedure
        #~ distmat                 prints out a pairwise distance matrix between the edges
        #~ filter                  filters one or more placefiles by placement name
        #~ info                    writes the number of leaves of the reference tree and the number of pqueries
        #~ merge                   merges placefiles together
        #~ redup                   restores duplicates to deduped placefiles
        #~ round                   clusters the placements by rounding branch lengths
        #~ taxtable                makes SQL enabling taxonomic querying of placement results
        #~ to_jplace                 converts old-style .place files to .jplace placement files

    
    
    #- Calculate EDPL values (edpl) -----------------------------------------
    print "-" x 80, "\n";
    print "Generating EDPL table...\n";
    
    chdir $output_dir or die "Can't change working dir to: $!\n"; 
    @args = (catfile($pplacer_dir, 'guppy'),
                    'edpl',
                    basename($aln).".jplace");
    splice(@args, 2, 0, "--pp") if $pp;
    print "@args\n" if $verbose;
    system(@args) == 0 or die "Unable to run: @args\n";
    chdir $cwd or die "Can't change back to working dir: $!\n";
    

    #- Generate tree with fattened edges (fat) -----------------------------
    print "-" x 80, "\n";
    print "Generating tree with fattened edges (fat)...\n";

    @args = (catfile($pplacer_dir, 'guppy'),
                    'fat',
                    #'--width-factor', '1',
                    '-o', $aln.".fat.xml", 
                    $aln.".jplace");
    splice(@args, 2, 0, "--pp") if $pp;
    print "@args\n" if $verbose;
    system(@args) == 0 or die "Unable to run: @args\n";



    #~ #- Generate one tree for each query sequence showing uncertainty ----
    #~ print "-" x 80, "\n";
    #~ print "Generating one tree for each query sequence showing uncertainty...\n";
    #~ 
    #~ @args = (catfile($pplacer_dir, 'guppy'),
                    #~ 'sing',
                    #~ $aln.".jplace");
    #~ splice(@args, 2, 0, "--pp") if $pp;
    #~ print "@args\n" if $verbose;
    #~ system(@args) == 0 or die "Unable to run: @args\n";



    #- Generate tree with each read represented as a pendant edge (tog) ---
    print "-" x 80, "\n";
    print "Generating tree with each read represented as a pendant edge (tog)...\n";
    chdir $output_dir or die "Can't change working dir to: $!\n"; 
    @args = (catfile($pplacer_dir, 'guppy'),
                    'tog',
                    '--xml',
                    #'-o', $aln.".tog.xml",
                    basename($aln).".jplace");
    splice(@args, 2, 0, "--pp") if $pp;
    print "@args\n" if $verbose;
    system(@args) == 0 or die "Unable to run: @args\n";
    chdir $cwd or die "Can't change back to working dir: $!\n";


    #~ #- Generate tree file with barycenters (bary) --------------------------
    #~ print "-" x 80, "\n";
    #~ print "Generating tree file with barycenters (bary)...\n";
    #~ 
    #~ @args = (catfile($pplacer_dir, 'guppy'),
                    #~ 'bary',
                    #~ #'-o', $aln.".bary.xml", 
                    #~ $aln.".jplace");
    #~ splice(@args, 2, 0, "--pp") if $pp;
    #~ print "@args\n" if $verbose;
    #~ system(@args) == 0 or die "Unable to run: @args\n";

}
#=======================================================================
#=======================================================================
#=======================================================================

 
#- End of script. ------------------------------------------------------
print "-" x 80, "\n";
print "PhyloAssigner finished successfully.\n\n";
print "All output files stored in:\n$output_dir\n\n";

#~ print "\n";
#~ print $placefile, "\n";
#~ print $tabfile, "\n";

if ( $jplaceoutfile or $taboutfile ) {
    print "Finally moving files...\n";
}

if ( $jplaceoutfile ) {
    copy($placefile, $jplaceoutfile) if -e $placefile or die "ERROR: Unable to generate place file at: $jplaceoutfile\n";
    print "Output place file: $jplaceoutfile\n";
}

if ( $taboutfile ) {
    copy($tabfile, $taboutfile) if -e $tabfile or die "ERROR: Unable to generate tab file at: $taboutfile\n";
    print "Output tab file: $taboutfile\n";
}
print "\n";

=begin comment
=end comment
=cut

#= Subroutines =========================================================

# Remove bad characters from fasta
sub process_query_aln {
    my $infile = shift @_;
    my $outfile = shift @_;
    
    my @new;
    
    #read infile
    open (my $in_fh, "<", $infile) or die "Could not open $infile\n$!\n";
    while (my $line = <$in_fh>) {
        $line =~ s/\./-/g;
        push(@new, $line);
    }
    close $in_fh;
    
    #write outfile
    open (my $out_fh, ">", $outfile) or die "Could not open $outfile for writing\n$!\n";
    print {$out_fh} @new;
    close $out_fh;
    
    return 1;
}


#get eviromental path for a given file (if exists)
sub get_PATH {
    my $query_file = shift @_;
    foreach ( split(':', $ENV{'PATH'}) ) {
        return $_ if -e  catfile($_,  $query_file);
    }
}




