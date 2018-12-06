# phyloassigner

                    P h y l o A s s i g n e r   R E A D M E
                                
                                17th July 2012


                        I. What is PhyloAssigner?

PhyloAssigner is a software pipeline allowing phylogenetic placement and 
taxonomic classification of nucleotide sequences. It was primarily 
developed to facilitate analyses of amplicon pyrotags from environmental 
samples, but may be useful in other studies as well. In a way, it is a 
convenience wrapper around the phylogenetic placement tool pplacer
(http://matsen.fhcrc.org/pplacer/); for more details of this method, see
the pplacer publications and those on the EPA method available in RAxML
(http://www.exelixis-lab.org/).

The program is written in Perl and is released under GNU GENERAL PUBLIC 
LICENSE v3 (see file COPYING). So, feel encouraged to work with the code 
yourself. It is currently developed and tested on Linux systems (we are 
using Ubuntu).

Copyright (C) 2012 Fabian Kilpert, Bank Beszteri, Adam Monier.

When you use the pipeline and / or the marine bacterial reference data set, 
please cite:

Vergin KL, Beszteri B, Monier A, Thrash JC, Temperton B, Treusch A, Kilpert F,
Worden AZ & Giovannoni SJ (2013). High resolution SAR11 ecotype dynamics at the 
Bermuda Atlantic Time-series Study Site by phylogenetic placement of 
pyrosequences. ISME Journal, doi:10.1038/ismej.2013.32.

When using the marine eukaryotic reference data set, please cite:

Wolf, C., Frickenhaus, S., Kilias, E.S., Peeken, I. and Metfies, K. (2013).
Regional variability in eukaryotic protist communities in the Amundsen Sea. Antarctic 
Science, in press.


                            II. Getting started

In case you downloaded the packed "tag.gz" archive, unpack it using your 
Linux shell:

    tar  xfvz  phyloassigner.tar.gz 

The new ./phyloassigner folder contains the main script 
phyloassigner.pl. Use

    perl  phyloassigner.pl  -h

in order to show the quick help and look-up the options available by the 
program. Option "--help" will show you the extended help, which will 
bring you back to this document (README).
It is also possible to download the latest version of PhyloAssigner via 
Apache Subversion (SVN). Go to 
http://aforge.awi.de/gf/project/phyloassigner on further instructions 
how to actually do the check-out.
The PhyloAssigner software should be accessible on your machine by now. 
However, as the program utilizes other packages (for good reasons), we 
are not done yet. 



                    III. Dependencies on other packages

PhyloAssigner is designed of Perl code (e.g. assigning taxonomic 
labels to individual sequences based on their placement on a reference 
tree), but also clamps proven external applications to its pipeline.  

Absolutely required are:

BioPerl 1.6.1       Linux distribution package manager installation
                    recommended (BioPerl 1.4 will NOT work! and problems
		    have been reported with BioPerl 1.6.9 as well; we
		    are working on the latter). The simplest way to install
		    the package is to download and unpack BioPerl-1.6.1.tar.gz
		    from 
		    http://www.bioperl.org/wiki/Getting_BioPerl#Bioperl_1.6.1.2C_Stable_Release
		    and specify its location in the system variable PERL5LIB
		    (e.g., export PERL5LIB=/home/user/BioPerl-1.6.1 when using
		    bash).
JSON 2.53           The simplest way to install the JSON::XS Perl module is to 
                    do:
                    cpan -i JSON
pplacer 1.1         http://matsen.fhcrc.org/pplacer/
                    version alpha 10 or higher required!
		    We now provide the pplacer executables as part of the 
		    PhyloAssigner distribution package (courtesy E. Matsen),
		    so you do not need to worry about this. More current versions 
		    of pplacer are not supported by PhyloAssigner.

And either one or both of the following tools for aligning queries to
reference profiles:

HMMER 3.0           http://hmmer.janelia.org/
                    Download and installation as recommended
mothur              As an alternative to HMMER, the mothur aligner can
                    be used: see http://www.mothur.org/.


Only required for building new reference databases:

PhyML 3             http://www.atgc-montpellier.fr/phyml/
                    Download and installation as recommended
                    (PhyML is NOT always required)

If you just want to perform an analysis with a precompiled database 
(e.g. from download), you won't need PhyML.

PhyloAssigner needs these pieces of software to be installed and their parent folders
be in the path (HMMER / mothur, PhyML) or in the Perl library search path (BioPerl,
JSON::XS), respectively.


                    IV. Convenience application paths

For convenience, PhyloAssigner will try to locate the executable 
binaries (HMMER, pplacer, PhyML, etc.) in your PATH enviroment
variable or - if unsuccessful - in your current working directory.
If it finds the binaries, running an analysis should be as simple as in 
this example:

    perl  phyloassigner.pl  example.phyloassignerdb.tar.gz  
    ./example/queries.fas

If the binaries are located elsewhere - which is the normal case if you 
just downloaded and installed the external applications - you must 
specify the paths to the directories which contain the external 
binaries, e.g.:

    perl  phyloassigner.pl  --hmmerdir /path/to/hmmer/binaries  
    --pplacerdir /path/to/pplacer  example.phyloassignerdb.tar.gz  
    ./example/queries.fas

PhyloAssigner takes use of more than one binary from each of these 
external packages. That is why we (the coders) decided to let 
PhyloAssigner ask (options: --hmmerdir, --pplacerdir, --phymldir) for 
the DIRECTORY containing these binaries and not for each binary itself. 
The fewer paths stated the better. Use option -h for quick help on 
options.
Again, if you use PhyloAssigner a lot, you should definitely think about
setting PATHs enviroment variables (e.g. "export PATH=$PATH:/path/to/
pplacer") or copying all binaries in your current working directory from
which you usually start PhyloAssigner. No need to worry about these 
paths any more. We assume one of these convenience solutions for being 
used in the examples of this document!



                        V. Running PhyloAssigner

The most common case should be that one wants to run an analysis with 
own data (FASTA file with nucleotide sequences) on a PhyloAssigner 
reference database (.phyloassignerdb) either from a publication download 
or compiled by fellow experts beforehand. You will find an example 
database in the installation directory of PhyloAssigner. Give it a try:

    perl  phyloassigner.pl  example.phyloassignerdb.tar.gz  queries.fas 
    
Keep in mind that you have to use paths valid to your local system:

    perl  /path/to/phyloassigner/phyloassigner.pl  
    /path/to/phyloassigner/example.phyloassignerdb.tar.gz  
    /path/to/phyloassigner/queries.fas 

The command 'perl' launches the main script of PhyloAssigner, 
phyloassigner.pl, which is written in Perl. The script itself is invoked 
with two parameters, which tell PhyloAssigner to analyse the nucleotide 
sequences in the FASTA file 'queries.fas' by using the reference data in
the preformatted database 'example.phyloassiger.tar.gz'. 
The file extension 'tar.gz' signals that this database is in fact a gzip 
compressed folder. PhyloAssigner also accepts the database if it is 
unpacked:

    perl  phyloassigner.pl  example.phyloassignerdb  queries.fas

It should usually suffice to use the compressed database, though. Being 
a data product consisting of only one file, it should be much more 
convenient to work with or hand over to colleagues.



                        VI. Pipeline workflow

The script starts the pipeline using the reference data in the 
reference database. If it is compressed it is decompressed in the 
current working directory first. The query sequences in the FASTA file 
are aligned (using hmmalign) to the reference profile HMM (HMMER 3). The 
aligned sequences are placed (pplacer 1.1) on the phylogenetic tree. 
pplacer outputs a number of possible placements in the tree for every 
query sequence - each with likelihood, ML weight ratio, and labelling 
according to the reference taxonomy mapping. Instead of just taking the 
top hit, which might only have a marginally better likelihood than the 
next best hit, the pipeline identifies the last common ancestor (LCA) of 
a range of hits. The sum of ML weight ratios (relative likelihoods) of 
all hits is 1 by definition. A threshold considers all best hits that 
fall within a given LCA threshold sum of 0.85 (user-definable) and 
determines the LCA taxon labelling. In doing so the phylogenetic 
assignment becomes more robust, while uncertainty is reflected by 
placements lower in the reference tree.



                            VII. Output files

PhyloAssigner puts all generated files in an output folder, which by 
default gets the same name as the input query file but with an '.output' 
extension. It contains:

.fat.xml        PHYLOXML Reference tree showing the number of placed     
                query sequences as branch widths (best viewed with     
                Archaeopteryx tree viewer (forester.jar))
.tog.xml        PHYLOXML Reference tree showing the exact location and 
                branch length for every placed sequence
.edpl           Text file with EDPL values for each sequence
.counts         Text file with overview on assigned taxon frequencies
.tab            Text file showing parsed placement results with         
                assignments of sequences to taxa (best and LCA)
.jplace         Main output file (JSON XML) of pplacer containing     
                multiple placements on the tree for every query sequence
.aln            Input reference alignment (FASTA)

PhyloAssigner generates many other files during runtime, which are 
needed temporarily and get deleted when the program finishes. Use option 
'--nocleanup' to keep all files in the end.



                VIII. Compiling a new phyloassignerdb

The PhyloAssigner package comes with an extra script for building new 
reference databases from scratch, setupdb.pl. The script also gives help 
on its usage when called with the '-h' option. Use it like here 
exemplarily demonstrated for the files provided in example subfolder of 
the PhyloAssigner installation (adjust paths to your system):

    perl  setupdb.pl  example/example.nwk  example/example.fna  
    "Uncultured Alveolates"

The Perl script setupdb.pl requires for input:

- NEWICK tree (here example.nwk) containing a pre-calculated tree 
  topology of reference sequences
- reference alignment in FASTA format (here example.fna)
- valid outgroup name for re-rooting the tree if necessary (here 
  "Uncultured Alveolates")

Note: Sequence names of reference tree and alignment must match 
perfectly, not tolerating any discrepancies in naming!! Really!!!
It is therefore recommended to export tree and alignment from a reliable 
database management system, e.g. Arb (http://www.arb-home.de/). It not 
only helps in defining named subclades but also naturally keeps sequence 
names identical in all files. As these files are text "only", you 
can of course create and rework them by hand or using your favourite 
programs (tree builders, aligners, text editors etc.). But you must 
conscientiously take care yourself to not produce files that somehow 
became corrupt!

We recommend using reference data sets with around 1000-1500 sequences, 
this enables high taxonomic resolution at acceptable computational 
costs. Of course the best number for your project will depend on your 
taxa of interest, availability of reference sequences for them and the 
computational capacity at your disposal.

The reference database setup script substitutes Us by Ts, trims the 
alignment to a region, and discards gap-only positions. Alignments 
exported from e.g. Silva (http://www.arb-silva.de) always have 
containing many gap-only positions with no meaningful information for 
that reason. The setup script generates and calibrates a profile HMM 
(HMMER 3), optimizes model parameters and branch lengths on the given 
phylogeny using PhyML 3, and finally generates a table mapping pplacer 
node IDs to taxonomic names from the labelled tree. The output database 
has the .phyloassigner extension or .phyloassigner.tar.gz if packed 
(default). Both can be used right away for the analysis with 
PhyloAssigner (see chapter V.). 



                        IX. phyloassignerdb content

A newly generated phyloassignerdb is folder containing several files:

REF.aln                     Reference alignment
REF.hmm3                    Profile HMM made from the reference alignment
REF.phy_phyml_tree.txt      Reference phylogenetic tree 
REF.phy_phyml_stats.txt     Reference phylogenetic tree statistics file

It is recommended adding a file with additional information on the 
dataset about content, authors, publication, license, etc. manually.
We suggest using a new file for that, ABOUT.txt. You can use

    tar xfvz example.phyloassignerdb.tar.gz

to unpack it, make your additions to the folder, and use

    tar cfvz example.phyloassignerdb.tar.gz example.phyloassignerdb

to pack it again.


                            X. Tips and Issues
                        
- Sequence names must contain word characters only! These 
    are letters of both cases, numbers, and underscores.
    regex: [a-zA-Z0-9_]
    
- Sequence names length is limited to 120 characters. You can yourself
    raise this value yourself by editing tools/fas2phy.pl
    
- PhyloAssigner is designed to continue with the files in its current
    working folder if they exist! This is meant as a convenience 
    solution allowing to skip time consuming processing steps if 
    PhyloAssigner it gets interrupted by a crash or other reasons. 
    IT MAY BE NECESSARY TO DELETE THE LATEST FILES OR THE WHOLE FOLDER
    BY HAND when data has become corrupted! Please, try a complete new
    run - deleting old outputs first - when there are unknown processing
    issues!

                        
                        
                    XI. Available PhyloAssiger databases
                    
                                (more to come)



