# ConSurf
[Link to Elixir](https://bio.tools/tool/mytum.de/ConSurf/1.0.6)  
Taken from [http://consurf.tau.ac.il/overview.html](http://consurf.tau.ac.il/overview.html)

The ConSurf server [1] is a bioinformatics tool for estimating the evolutionary conservation of amino/nucleic acid positions in a protein/DNA/RNA molecule based on the phylogenetic relations between homologous sequences. The degree to which an amino (or nucleic) acid position is evolutionarily conserved is strongly dependent on its structural and functional importance; rapidly evolving positions are variable while slowly evolving positions are conserved. Thus, conservation analysis of positions among members from the same family can often reveal the importance of each position for the protein (or nucleic acid)'s structure or function. In ConSurf, the evolutionary rate is estimated based on the evolutionary relatedness between the protein (DNA/RNA) and its homologues and considering the similarity between amino (nucleic) acids as reﬂected in the substitutions matrix [2,3]. One of the advantages of ConSurf in comparison to other methods is the accurate computation of the evolutionary rate by using either an empirical Bayesian method or a maximum likelihood (ML) method [3].

1. Glaser, F., Pupko, T., Paz, I., Bell, R.E., Bechor-Shental, D., Martz, E. and Ben-Tal, N. (2003) Bioinformatics, 19, 163-164.

2. Mayrose,I., Graur,D., Ben-Tal,N. and Pupko,T. (2004) Mol. Biol. Evol., 21, 1781-1791.

3. Pupko,T., Bell,R.E., Mayrose,I., Glaser,F. and Ben-Tal,N. (2002) Bioinformatics, 18 S71-S77.

## HOWTO Install

Below You will find the guidlines on how to succesfully install ConSurf on your machine. Current manual is designed for Debian OS.

### Installing libraries

In order for ConSurf to work properly please, install the following libraries:

```shell
autoconf
g++
libexpat-dev
libtool
perl
```
or simply run
```shell
sudo apt-get install autoconf g++ libexpat-dev libtool perl -y
```

#### Perl and BioPerl

It is mandatory to have perl installed on your machine.
After Perl is installed please run the following commands in order to download and install perl modules
```shell
cpan
```
in CPAN run the following commands
```shell
sudo cpan force install Bio::Perl Config::IniFiles List::Util
```

### 3rd Party Software

ConSurf uses several programs to calculate the outputs. The Copyrights to these programs belongs to their owner. PLEASE MAKE SURE YOU FOLLOW THE LICENSE BY EACH OF THESE PROGRAMS.

The software can be found in the folder `3rd party software`. Here are the instructions to install the software itself.

#### Rate4Site

Rate4Site - is the main algorithm behind ConSurf.

**Installation**:
The sources are already in the folder `3rd party software`. ConSurf requires the *fast* and *slow* versions. In order to obtain both of them simply run `rate4site_build.sh` script:

```shell
./rate4site_build.sh
```

#### ClustalW

ClustalW - please [download from EBI; `ftp://ftp.ebi.ac.uk/pub/software/clustalw2/`](ftp://ftp.ebi.ac.uk/pub/software/clustalw2/)

In `3rd party software` folder You can already find the executable.

#### blastpgp

blastpgp can be found in `3rd party software/blast-2.2.26/bin`

#### CD-HIT

cd-hit can be found in `3rd party software/cdhit`

To compile simply run:
```shell
./cd-hit_build.sh
```
("Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences", Weizhong Li & Adam Godzik Bioinformatics, (2006) 22:1658-9)

#### MUSCLE MSA program

can be found in `3rd aprty software`

("MUSCLE: a multiple sequence alignment method with reduced time and space complexity". Edgar R.C. (2004), BMC Bioinformatics 5: 113.)

####  Databases

Please download **SwissProt and TrEMBL databases** and format the databases for Blast search (using blast formadb which is part of blast system. [manual](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/doc/blast/formatdb.html) )

  0. update the location of *SwissProt* to point the location of `SWISSPROT_DB`
  0. update the location of *TrEMBL* to point the location of `UNIPROT_DB`.
  
Databases can be downloaded in FASTA format from [here] (http://www.uniprot.org/downloads) or just run the following commands:
```shell
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
gunzip uniprot_trembl.fasta.gz
```

### Configuration

In order to run ConSurf is necessary to update the configuration file - `consurfrc.default`

1. In `[programs]` section please enter the full path to the executables of 3rd party software covered earlier  (Example: `MUSCLE=\home\path_to_ConSurf\3rd party software\muscle3.8.31`)
2. in `[databases]` section please enter the full path to the DBs.   (Example:`SWISSPROT_DB=\home\path_to_DB_Folder\uniprot_sprot.fasta`)

### Installing ConSurf

To install ConSurf run the following commands in ConSurf folder:
```shell
aclocal
autoheader
automake --force-missing --add-missing
autoconf
./configure
sudo make
sudo make install
```

## HOWTO Run

The ConSurf works in several modes

1. Given Protein PDB File.
2. Given Multiple Sequence Alignment (MSA) and Protein PDB File.
3. Given Multiple Sequence Alignment (MSA), Phylogenetic Tree, and PDB File.

The script is using the user provided MSA (and Phylogenetic tree if available) to calculate the conservation score for each position in the MSA based on the Rate4Site algorithm [2].

When running in the first mode, the script builds the MSA for the given protein based on the ConSurf protocol.

Usage:

```shell
consurf -PDB <PDB FILE FULL PATH>  -CHAIN <PDB CHAIN ID> -Out_Dir <Output Directory>
```

Requirements for successful run: 

1. User must have `<PDB_ID>_<CHAIN_ID>_protein_seq.fas` in the directory specified as `-workdir`. The file should contain the protein sequence in FASTA format. More information and the solution is available in [Issue #14](https://github.com/Rostlab/ConSurf/issues/14)  
2. User is strongly encouraged to specify the `<PDB_ID>_<CHAIN_ID>.protein_query.blast` file. It is created in the directory specified as `-workdir` by the ConSurf, but sometimes ConSurf might fail to find it. More information and the solution is available in [Issue #17](https://github.com/Rostlab/ConSurf/issues/17)

### MANDATORY INPUTS

`-PDB <PDB FILE FULL PATH> - PDB File`  
`-CHAIN <PDB CHAIN ID> - Chain ID`  
`-Out_Dir <Output Directory> - Output Path`

### MSA Mode (Not Using -m)

`-MSA <MSA File Name>`	(MANDATORY IF -m NOT USED)  
`-SEQ_NAME <"Query sequence name in MSA file">`  (MANDATORY IF -m NOT USED)  
`-Tree <Phylogenetic Tree (in Newick format)>` (optional, default building tree by Rate4Site).

### Building MSA (Using -m)

```shell
-m Builed MSA mode  
-MSAprogram ["CLUSTALW"] or ["MUSCLE"] # default: MUSCLE
-DB ["SWISS-PROT"] or ["UNIPROT"] # default: UniProt
-MaxHomol <Max Number of Homologs to use for ConSurf Calculation> #deafult: 50
-Iterat <Number of PsiBlast iterataion> # default: 1
-ESCORE <Minimal E-value cutoff for Blast search> # default: 0.001
-BlastFile <pre-calculated blast file provided by the user>
```

### Rate4Site Parameter

(see http://consurf.tau.ac.il/overview.html#methodology and http://consurf.tau.ac.il/overview.html#MODEL for details)

```shell
-Algorithm [LikelihoodML] or [Bayesian] # default: Bayesian
-Matrix [JTT] or [mtREV] or [cpREV] or [WAG] or [Dayhoff] # default JTT
```


### Examples

1. Basic Build MSA mode (using defaults parameters)
	`consurf -PDB  example/1lk2.pdb -CHAIN A -Out_Dir output -m`

Example:
```shell
cd ConSurf
consurf -PDB example/1lk2.pdb -CHAIN A -Out_Dir output -m --workdir workdir -BlastFile example/1lk2_A.protein_query.blast
```
## Method Description

Taken from [http://consurf.tau.ac.il/overview.html](http://consurf.tau.ac.il/overview.html)

Given the amino or nucleic acid sequence (can be extracted from the 3D structure), ConSurf carries out a search for close homologous sequences using BLAST (or PSI-BLAST) [4,5]. The user may select one of several databases and specify criteria for defining homologues. The user may also select the desired sequences from the BLAST results. The sequences are clustered and highly similar sequences are removed using CD-HIT [6]. A multiple sequence alignment (MSA) of the homologous sequences is constructed using MUSCLE(default) or CLUSTALW. The MSA is then used to build a phylogenetic tree using the neighbor-joining algorithm as implemented in the Rate4Site program [7]. Position-specific conservation scores are computed using the empirical Bayesian or ML algorithms [2,3]. The continuous conservation scores are divided into a discrete scale of nine grades for visualization, from the most variable positions (grade 1) colored turquoise, through intermediately conserved positions (grade 5) colored white, to the most conserved positions (grade 9) colored maroon. The conservation scores are projected onto the protein/nucleotide sequence and on the MSA.

1. Glaser, F., Pupko, T., Paz, I., Bell, R.E., Bechor-Shental, D., Martz, E. and Ben-Tal, N. (2003) Bioinformatics, 19, 163-164.

2. Mayrose,I., Graur,D., Ben-Tal,N. and Pupko,T. (2004) Mol. Biol. Evol., 21, 1781-1791.

3. Pupko,T., Bell,R.E., Mayrose,I., Glaser,F. and Ben-Tal,N. (2002) Bioinformatics, 18 S71-S77.

4. Altschul,S.F., Wootton,J.C., Gertz,E.M., Agarwala,R., Morgulis,A., Schaffer,A.A. and Yu,Y.K. (2005) FEBS J., 272, 5101-5109.

5. Altschul,S.F., Madden,T.L., Schaffer,A.A., Zhang,J., Zhang,Z., Miller,W. and Lipman,D.J. (1997) Nucleic Acids Res., 25, 3389-3402.

6. Li,W. and Godzik,A. (2006) Bioinformatics, 22, 1658-1659.

7. Pupko, T., Bell, R.E., Mayrose, I., Glaser, F. and Ben-Tal, N. (2002) Bioinformatics, 18, S71-77.


## HOWTO Using Vagrant

In the folder `vagrant_ConSurf` You can find the Vagrantfile, file necessary to run a Vagrant Box with ConSurf on it.
Intallation of Vagrant will not be explained in here.

In order to start the Box, execute the following:
```shell
cd /<path_to_dir>/vagrant_ConSurf
vagrant up --no-destroy-on-error
vagrant ssh
```

After the box has been started and You are connected to the VM through SSH, run `export LC_ALL=C`.
All dependency libraires from section *Installing Libraries* should have been installed during execution of `vagrant up`, rest can be set up according to *Installation* guidelines presented earlier.
