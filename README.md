# ConSurf

## Installation

ConSurf uses several programs to calculate the outputs. The Copyrights to these programs belongs to their owner. PLEASE MAKE SURE YOU FOLLOW THE LICENSE BY EACH OF THESE PROGRAMS.

After installing these programs please update the relevant lines on the `CONSURF_CONSTANTS.pm` file accordingly:

0. **Rate4Site** - is the main algorithm behind ConSurf. Please download the [Rate4Site program](http://www.tau.ac.il/~itaymay/cp/rate4site.html) and follow the installation instructions relevant to your operating system. Please install both _Rate4Site_slow_ (using `Makefile_slow`) and _Rate4Site_ (Makefile).
  0. Update the line that points to `RATE4SITE` in `CONSURF_CONSTANTS.pm` file; replace `/db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4site.exe` with the location path to your  `Rate4Site` installation.
  0. Update the line points to _RATE4SITE_SLOW_; `/db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4siteSlow.exe` with the location of your `rate4site_Slow` installation
0. **ClustalW** - please [download from EBI](ftp://ftp.ebi.ac.uk/pub/software/clustalw2/):
  0. update the line `/usr/local/bin/clustalw` to point your installation.
0. In `ConSurf.pl` please replace the `/usr/bin/perl` with the relevant perl command in your system on the line #!/usr/bin/perl -w. (please noticed that ConSurf requiered that perl version above 5 is installed)
0. In `ConSurf.pl` please replace `/groups/bioseq.home/HAIM/ConSurf_Exec` with the location of ConSurf Script in your system
0. The ConSurf system requiered BioPerl to be installed on your system

### (OPTONAL) The programs below are required by the mode that automatically created the MSA

0. **blastpgp** - please download local version of the NIH blast system and install it locally.
  0. please update the line `/usr/local/bin/blastpgp` to point your installation
0. **cd-hit** - please download CD-Hit program and install it locally  
("Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences", Weizhong Li & Adam Godzik Bioinformatics, (2006) 22:1658-9)
  0. please update the line `/db1/Local/src/cd-hit_redundency/` to point your local installation directory
0. Please download **SwissProt and TrEMBL databases**
  0. format the databases for Blast search (using blast formadb which is part of blast system)
  0. update the location of `/biodb/BLAST/Proteins/swissprot` to point the locaation of swissprot
  0. update the location of `/biodb/BLAST/Proteins/uniprot` to point the location of TrEMBL.
0. Download **MUSCLE MSA program**
("MUSCLE: a multiple sequence alignment method with reduced time and space complexity". Edgar R.C. (2004), BMC Bioinformatics 5: 113.) and install it on your system.
  0. update `/usr/local/bin/muscle` to point your installation


##Usage

The ConSurf work in several modes

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

## MANDATORY INPUTS

`-PDB <PDB FILE FULL PATH> - PDB File`  
`-CHAIN <PDB CHAIN ID> - Chain ID`  
`-Out_Dir <Output Directory> - Output Path`

## MSA Mode (Not Using -m)

`-MSA <MSA File Name>`	(MANDATORY IF -m NOT USED)  
`-SEQ_NAME <"Query sequence name in MSA file">`  (MANDATORY IF -m NOT USED)  
`-Tree <Phylogenetic Tree (in Newick format)>` (optional, default building tree by Rate4Site).

## Building MSA (Using -m)

```shell
-m Builed MSA mode  
-MSAprogram ["CLUSTALW"] or ["MUSCLE"] # default: MUSCLE
-DB ["SWISS-PROT"] or ["UNIPROT"] # default: UniProt
-MaxHomol <Max Number of Homologs to use for ConSurf Calculation> #deafult: 50
-Iterat <Number of PsiBlast iterataion> # default: 1
-ESCORE <Minimal E-value cutoff for Blast search> # default: 0.001
```

## Rate4Site Parameter

(see http://consurf.tau.ac.il/overview.html#methodology and http://consurf.tau.ac.il/overview.html#MODEL for details)

```shell
-Algorithm [LikelihoodML] or [Bayesian] # default: Bayesian
-Matrix [JTT] or [mtREV] or [cpREV] or [WAG] or [Dayhoff] # default JTT
```


## Examples

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
