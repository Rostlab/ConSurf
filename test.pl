use Getopt::Long
 
my %FORM=();
my %VARS=();

#GetOptions(
#           "PDB=s"=>\$VARS{pdb_file_name},
#	    "CHAIN=s"=>\$FORM{chain},
#	    
#	    "MSAprogram:s"=>\$FORM{MSAprogram},
#	    "DB:s"=>\$FORM{database}, #"SWISS-PROT"
#	    "MaxHomol:s"=>\$FORM{MAX_NUM_HOMOL},
#	    "Iterat:s"=>\$FORM{iterations},
#	    "ESCORE:s"=>\$FORM{ESCORE},
#	    "SEQ_NAME:s"=>\$FORM{msa_SEQNAME},
#	    "Algorithm:s"=>\$FORM{algorithm},
	    
	#    "Out_Dir=s"=>\$VARS{working_dir},
	    
	#    "MSA:s"=>\$FORM{MSA},
	#    "Matrix:s"=>\$FORM{matrix},
	    
	#    "Tree:s"=>\$VARS{user_tree_file_name},
	    
#	    "m"=>\$FORM{buildMSA});
	    
my $command='cd /groups/bioseq.home/HAIM/ConSurf_Exec/TESTS/; perl test1.pl';
print $command;
`$command`;
# GetOptions(
#             "PDB=s"=>\$FORM{PDB},
# 	    "CHAIN=s"=>\$FORM{chain},
# 	    "MSA:s"=>\$FORM{MSA},
#             "MSAprogram:s"=>\$FORM{MSAprogram},
# 	    "Matrix:s"=>\$FORM{matrix},
# 	    "m"=>\$FORM{buildMSA});
# # 	    "int=i"=> \$moo{mandatoryinteger},
# #             "optint:i"=> \$moo{optionalinteger},
# #             "float=f"=> \$moo{mandatoryfloat},
# #             "optfloat:f"=> \$moo{optionalfloat});

foreach (keys %FORM) {
 print "$_ = *$FORM{$_}*\n";
}
foreach (keys %VARS) {
 print "$_ = *$VARS{$_}*\n";
}
print "Unprocessed by Getopt::Long\n" if $ARGV[0];
foreach (@ARGV) {
  print "$_\n";
}

# MUSCLE
# if($FORM{MSAprogram} eq 'CLUSTALW'){
# 
#  if ($FORM{algorithm} eq "LikelihoodML"){
#  $FORM{matrix}
#  $FORM{chain}
#  
