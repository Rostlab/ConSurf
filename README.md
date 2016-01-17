# ConSurf

## HOWTO Install

Below You will find the guidlines on how to succesfully install ConSirf on Your machine. Current manual is designed for Debian OS.

### Installing libraries

In order for ConSurf to work properly please install the following libraries:

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

It is mandatory to have perl installed on Your machine.
After Perl is installed please run the following commands in order to download and install perl modules
```shell
cpan
```
in CPAN run the following commands
```shell
force install Bio::Perl
force install Config::IniFiles
force install List::Util
```

### 3rd Party Software

ConSurf uses several programs to calculate the outputs. The Copyrights to these programs belongs to their owner. PLEASE MAKE SURE YOU FOLLOW THE LICENSE BY EACH OF THESE PROGRAMS.

The software can be found in the folder `3rd party software`. Here are the instructions to install the software itself.

#### Rate4Site

Rate4Site - is the main algorithm behind ConSurf.

**Installation**:
The sources are already in the folder `3rd party software`. ConSurf requires the *fast* and *slow* versions. In order to obtain both of them simply run `rate4site_build.sh` script:

```shell
(optional) chmod +x rate4site_build.sh
./rate4site_build.sh
```

#### ClustalW

ClustalW - please [download from EBI; `ftp://ftp.ebi.ac.uk/pub/software/clustalw2/`](ftp://ftp.ebi.ac.uk/pub/software/clustalw2/)

In `3rd paty software` folder You can already find the executable.

#### blastpgp

blastpgp can be found in `3rd aprty software/blast-2.2.26/bin`

#### CD-HIT

cd-hit can be found in `3rd aprty software/cdhit`

To compile simply run:
```shell
(optional) chmod +x cd-hit_build.sh
./cd-hit_build.sh
```
("Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences", Weizhong Li & Adam Godzik Bioinformatics, (2006) 22:1658-9)

#### MUSCLE MSA program

can be found in `3rd aprty software`

("MUSCLE: a multiple sequence alignment method with reduced time and space complexity". Edgar R.C. (2004), BMC Bioinformatics 5: 113.) 

####  Databases

Please download **SwissProt and TrEMBL databases** and format the databases for Blast search (using blast formadb which is part of blast system [manual](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/doc/blast/formatdb.html) )

  0. update the location of *SwissProt* to point the location of `SWISSPROT_DB`
  0. update the location of *TrEMBL* to point the location of `UNIPROT_DB`.

### Configuration

In order to run ConSurf is is necessary to update the configuration file - `consurfrc.default`

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

The script is using the user provided MSA (and Phylogenetic tree if available) to calculate the conservation score for each position in the MSA based on the Rate4Site algorithm
(Mayrose, I., Graur, D., Ben-Tal, N., and Pupko, T. 2004. Comparison of site-specific rate-inference methods: Bayesian methods are superior. Mol Biol Evol 21: 1781-1791).

When running in the first mode the scripts the MSA is automatically build the MSA for the given protein based on ConSurf protocol.

Usage:

```shell
ConSurf -PDB <PDB FILE FULL PATH>  -CHAIN <PDB CHAIN ID> -Out_Dir <Output Directory>
```

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
	`perl ConSurf.pl -PDB  MY_PDB_FILE.pdb -CHAIN MY_CHAIN_ID -Out_Dir /MY_DIR/ -m`
2. using build MSA mode and advanced options:
	`perl ConSurf.pl -PDB MY_PDB.pdb -CHAIN A -Out_Dir /MY_DIR/ -m -MSAprogram CLUSTALW -DB "SWISS-PROT" -MaxHomol 100 -Iterat 2 -ESCORE 0.00001 -Algorithm LikelihoodML -Matrix Dayhoff`  
	- This will run Consurf in building MSA mode (-m) for Chain A (-CHAIN) of PDB file: `/groups/bioseq.home/1ENV.pdb` (-PDB).  
	- The sequences for the MSA will be chosen according to iterations of PSI-Blast (-Iterat) against SwissProt Database (-DB "SWISS-PROT"), with E value cutoff of 0.00001 (-ESCORE), considering maximum 100 homologues (MaxHomol).  
	- The Rate4Site will use Maximum Liklihood algorithm (-Algorithm) and Dayhoff model (-Matrix)
3. Simple Run With prepared MSA.
	`perl ConSurf.pl -PDB MY_PDB_FILE.pdb -CHAIN A -Out_Dir /MY_DIR/ -MSA MY_MSA_FILE -SEQ_NAME MY_SEQ_NAME`

For any questions or suggestions please contact: bioSequence@tauex.tau.ac.il

## Method Description

## Evaluation
