#!/usr/bin/perl

package CONSURF_CONSTANTS; #don't forget: a package must end with a return value (1; in the end)!!!!!

use constant MAXIMUM_MODIFIED_PERCENT => 0.15;
use constant FRAGMENT_REDUNDANCY_RATE => 95;
use constant FRAGMENT_OVERLAP => 0.10;
use constant FRAGMENT_MINIMUM_LENGTH => 0.60;
use constant MINIMUM_FRAGMENTS_FOR_MSA => 5;
use constant LOW_NUM_FRAGMENTS_FOR_MSA => 10;

use constant BAYES_INTERVAL => 3;


#external databases
	# replace /biodb/PDB/data/structures/divided/pdb/ with your location of the PDB
use constant PDB_DIVIDED =>         "/mnt/project/rost_db/data/pdb/entries/";
	# replace /biodb/BLAST/Proteins/swissprot with your location of SwissProt database
use constant SWISSPROT_DB =>        "/mnt/project/rost_db/data/swissprot/uniprot_sprot";
	# replace /biodb/BLAST/Proteins/uniprot with your location of UniProt Database
use constant UNIPROT_DB =>       "/mnt/project/consurf/data/uniprot";

use constant UNIREF90_DB =>	    "/mnt/project/rost_db/data/big/big_80";

use constant CLEAN_UNIPROT_DB =>   "/mnt/project/consurf/data//clean_uniprot";

use constant RCSB_WGET =>           "wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/";

#external programs
	#replace /usr/local/bin/muscle with the location of your muscle instalation
use constant MUSCLE =>              "muscle";
	# replace /usr/local/bin/clustalw with your locationion of your clustalw istalation
#use constant CLUSTALW =>            "/mnt/home/gcelniker/ConSurf_TEST/DB/clustalw/";
use constant CLUSTALW =>            "clustalw";
	# replace /db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4site.exe with the location of your rate4site instalation
#use constant RATE4SITE =>           "/mnt/home/gcelniker/ConSurf_TEST/DB/r4s/rate4site";
use constant RATE4SITE =>           "/mnt/project/consurf/bin/rate4site";
	#replace /db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4siteSlow.exe with the location of your rate4site_Slow instalation
use constant RATE4SITE_SLOW =>      "/mnt/project/consurf/bin/rate4site";
	#replace /db1/Local/src/cd-hit_redundency/ with the directory containing your instalation of CD-HIT
use constant CD_HIT_DIR =>          "cd-hit";
	#replace /usr/local/bin/blastpgp with the location of your blastpgp instalation
use constant BLASTPGP => 	    "blastpgp";
# constant values
use constant BLAST_MAX_HOMOLOGUES_TO_DISPLAY => 50;

# external links
use constant RCSB_WEB =>                "http://www.rcsb.org/";


1;
