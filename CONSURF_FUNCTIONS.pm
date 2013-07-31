#!/usr/bin/perl

package CONSURF_FUNCTIONS; #don't forget: a package must end with a return value (1; in the end)!!!!!

#use strict;
use Fcntl ':flock'; # import LOCK_* constants
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::AlignI;

# global vars
my %aa;
my %modified_residues;
my $nucleic_acid_flag = "!";

my $bayesInterval=3;
my %ColorScale = (0 => 9, 1 => 8, 2 => 7, 3 => 6, 4 => 5, 5 => 4, 6 => 3, 7 => 2, 8 => 1);
my %tr_aa;
    $tr_aa{LYS}="K";$tr_aa{ARG}="R";$tr_aa{HIS}="H";$tr_aa{ASP}="D";$tr_aa{GLU}="E";
    $tr_aa{TYR}="Y";$tr_aa{TRP}="W";$tr_aa{SER}="S";$tr_aa{THR}="T";$tr_aa{PHE}="F";
    $tr_aa{LEU}="L";$tr_aa{ILE}="I";$tr_aa{MET}="M";$tr_aa{CYS}="C";$tr_aa{ASN}="N";
    $tr_aa{GLN}="Q";$tr_aa{ALA}="A";$tr_aa{VAL}="V";$tr_aa{PRO}="P";$tr_aa{GLY}="G";

#---------------------------------------------		
# gives the number in minimum 2 digits
sub new_num{
    my $num = shift;
    ($num < 10) ? return "0".$num : return $num;
}
#---------------------------------------------		
sub printTime {
   my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
   my $year = 1900 + $yearOffset;

   $second = &new_num($second);
   $minute = &new_num($minute);
   $hour = &new_num($hour);
   $month = &new_num($month+1);
   $dayOfMonth = &new_num($dayOfMonth);

   return "$hour:$minute:$second $dayOfMonth-".$month."-$year";
}
#---------------------------------------------
# input:  path to pdb file
# output: 3 options:
# 1. --PDB_NOT_OPEN if couldn't open the pdb file
# 2. --NO_CHAINS if no chain was founded in column 22
# 3. string with all the chains founded in this pdb.

sub which_chain_in_pdb_and_seqres{
    my $input_pdb = shift;
    my $chain_founded;
    my %all_chains;
    my @ret;
    my $seqres_found = "--SEQRES_no";
    
    unless (open PDB, $input_pdb){
        @ret = ("--PDB_NOT_OPEN $input_pdb $!");
        return \@ret;}
    while (<PDB>){
        if (/^ATOM/){
            $chain_founded = substr $_, 21, 1;
            if (!(exists $all_chains{$chain_founded})){
                $all_chains{$chain_founded} = 1;
            }
        }
        if ($seqres_found eq "--SEQRES_no" && /^SEQRES/){
            $seqres_found = "--SEQRES_yes";
        }
    }
    close PDB;    
    $chain_founded = "";
    foreach my $key (keys %all_chains){
        $chain_founded.=$key;
    }
    if($chain_founded !~ /\S/){
        @ret = ("--NO_CHAINS", $seqres_found);}
    else{
        @ret = ($chain_founded, $seqres_found);}
    return \@ret;
}
#---------------------------------------------
# input : 1. path to a pdb file, where there is no chain identifier in the 22 column of ATOM and 12 column of SEQRES
#         2. one letter denotes a chain identifier to add
# output : the same file, in the same path, where the letter given as input is added to the previously empty 22 column.
sub add_chain_to_pdb{
    my $input_pdb = shift;
    my $chain_id_to_add = shift;
    
    my ($beg_line, $end_line, $line);
    
    open PDB_IN, "+>>".$input_pdb;
    seek PDB_IN, 0, 0;
    my @all_lines_in_pdb = <PDB_IN>;
    truncate PDB_IN, 0;
    foreach(@all_lines_in_pdb){
        if (/^ATOM/){
            $line = $_;
            $beg_line = substr $line, 0, 21;
            $end_line = substr $line, 22, length($line);
            $_ = $beg_line.$chain_id_to_add.$end_line;
        }
        elsif (/^SEQRES/){
            $line = $_;
            $beg_line = substr $line, 0, 11;
            $end_line = substr $line, 12, length($line);
            $_ = $beg_line.$chain_id_to_add.$end_line;
        }
        print PDB_IN $_;
    }    
    close PDB_IN;    
}
#---------------------------------------------
sub convertNewline{
    # runs dos2unix, the program that converts plain text files in DOS/MAC format to UNIX format.
    my $inputFilePath = shift;
    my $WorkingDir = shift;
    my $dos2unix="cd $WorkingDir;dos2unix -q $inputFilePath";        
    system "$dos2unix";
    # if the input file was in mac format, the simple dos2unix will not work.
    # read the file - if it is only one line, it might mean that the new line characters
    # are not read well (for example: ^M). Trying to run dos2unix again, saying the format is mac
    $WorkingDir.='/' unless $WorkingDir =~ /\/$/;
    if (open FILE, $WorkingDir.$inputFilePath){
        my $num_of_lines = 0;
        while (<FILE>){
            $num_of_lines++;
        }
        close FILE;
        if ($num_of_lines==1){
            $dos2unix="cd $WorkingDir;dos2unix -c mac $inputFilePath -q ";
            system "$dos2unix";
        }
    }
    
}
#---------------------------------------------
# FROM rate4site_routines.pm
####################################################################################
# There are some tests to see if rate4site failed.
# Since we can't trust only one of them, we do all of them. If one of them is tested to be true - than a flag will get TRUE value
# 1. the .res file might be empty.
# 2. if the run failed, it might be written to the log file of r4s.
# 3. in a normal the r4s.log file there will lines that describe the grades. if it fail - we won't see them
# In one of these cases we try to run the slower version of rate4site.
# We output this as a message to the user.
sub check_if_rate4site_failed
{
    my $res_flag = shift; #path of the file, including working dir
    my $r4s_log = shift;    
    my $return_ans = "no";
    my $print_to_html = "";
    my $print_to_log = "";
    my $r4s_process_id = "";
    my $error_found = "no";
    
    if(!-e $res_flag){
        $print_to_log = "CONSURF_FUNCTIONS::check_if_rate4site_failed : the file $res_flag does not exsits. \n";
        $error_found = "yes";
    }
    elsif (-e $res_flag && -z $res_flag) #1
    {
        $print_to_log = "CONSURF_FUNCTIONS::check_if_rate4site_failed : the file $res_flag was found to be of size 0. \n";        
        $error_found = "yes";
    }
    if(-e $res_flag){
        unless (open R4S_RES, $res_flag){
            $print_to_log = "CONSURF_FUNCTIONS::check_if_rate4site_failed : can not open file: $res_flag. aborting.\n";
            $error_found = "yes";
        }
        while(<R4S_RES>){
            if(/In the tree file there is the name: (.+) that is not found in the sequence file/){
                $error_found = "yes";
                $print_to_log = "CONSURF_FUNCTIONS::check_if_rate4site_failed : sequence name $1 was found in the tree file, was not found in the MSA\n";
                last;
            }
        }
        close R4S_RES;
    }
    if (-e $r4s_log && !(-z $r4s_log)) #2,3
    {
        unless (open R4SLOG, $r4s_log) {
            $print_to_log = "CONSURF_FUNCTIONS::check_if_rate4site_failed : can not open file: $r4s_log. aborting.\n";
            $error_found = "yes";
        }
        while (<R4SLOG>)
        {
            if (/^.Process_id= (\d+)/){
                $r4s_process_id = $1;
            }
            if ($_ =~ m/likelihood of pos was zero/){
                $print_to_log = "CONSURF_FUNCTIONS::check_if_rate4site_failed : the line: \"likelihood of pos was zero\" was found in $r4s_log.\n";
                $error_found = "yes";
                last;
            }
            if ($_ =~ m/rate of pos\:\s\d\s=/){ #if we see this line, we change the flag{                
                $return_ans = "no";
                last;
            }
            if($_ =~ m/The amino-acid sequences contained the character: (.+)/){
                $print_to_log = "CONSURF_FUNCTIONS::check_if_rate4site_failed : illegal character $1 was found in the MSA\n";
                $error_found = "yes";
                last;
            }
            if($_ =~ m/Could not find a sequence that matches the sequence name/){
                my $seq_name = <R4SLOG>;
                $print_to_log = "CONSURF_FUNCTIONS::check_if_rate4site_failed : the submitted query sequence name $seq_name was not found in MSA";
                $error_found = "yes";
                last;
            }
            if($_ =~ m/The sequence name: (.+)was found in the tree file but not found in the sequence file/){
                my $seq_name = $1;
                $print_to_log = "CONSURF_FUNCTIONS::check_if_rate4site_failed : the sequence name $1 was found in the tree file, but not in the MSA";
                $error_found = "yes";
                last;
            }
            if ($_ =~ m/Bad format in tree file/){
                $print_to_log = "CONSURF_FUNCTIONS::check_if_rate4site_failed : Bad Tree Format";
                $error_found = "yes";
                last;
            }
            if ($_ =~ m/not all sequences are of the same lengths/){
                $print_to_log = "CONSURF_FUNCTIONS::check_if_rate4site_failed : problem with the MSA : not all sequences are of the same lengths";
                $error_found = "yes";
                last;
            }
        }
        close R4SLOG;
    }
    $return_ans = "yes" if ($error_found eq "yes");
    return ($return_ans, $print_to_log, $r4s_process_id, $print_to_html);
}

#---------------------------------------------

#############################################
#subroutines of rasmol_gradesPE_and_pipe.pm #
#############################################

#************************************************************************************************
# read the rate4site output into an array.
# calculates layers according to the max and min grades from the output
# $Output : a pointer to array ; assigns each position in the MSA (array index) with a grade (arracy value)

sub assign_colors_according_to_r4s_layers{	
    
	my ($rate4site_filename, $Output) = @_;
    
    unless (open RATE4SITE, $rate4site_filename){
        return ("assign_colors_according_to_r4s_layers : can't open $rate4site_filename","PANIC");}
    my $line;
    my $i = 0;
    while (<RATE4SITE>) {
        $line = $_;
        chomp $line;
		# baysean
        if ($line =~ /^\s+(\d+)\s+(\w)\s+(\S+)\s+\[\s*(\S+),\s*(\S+)\]\s+\S+\s+(\d+)\/(\d+)/){
                $Output->[$i]{POS} = $1;
                $Output->[$i]{SEQ} = $2;
                $Output->[$i]{GRADE} = $3;
                $Output->[$i]{INTERVALLOW} = $4;
                $Output->[$i]{INTERVALHIGH} = $5;
                $Output->[$i]{MSA_NUM}=$6;
                $Output->[$i]{MSA_DENUM}=$7;
                $i++;
        }
		# Maximum likelihood
        elsif($line=~m/(\d+)\s+(\w)\s+(\S+)\s+(\d+)\/(\d+)/){
            $Output->[$i]{POS} = $1;
            $Output->[$i]{SEQ} = $2;
            $Output->[$i]{GRADE} = $3;
            $Output->[$i]{INTERVALLOW} = $3;
            $Output->[$i]{INTERVALHIGH} = $3;
            $Output->[$i]{MSA_NUM}=$4;
            $Output->[$i]{MSA_DENUM}=$5;
            $i++;
        }
    }
    close RATE4SITE;
	
    my $element;
    my $max_cons = $Output->[0]{GRADE};
    my $ConsColorUnity; #unity of conservation to be colored
    foreach $element (@$Output){
        if ($$element{GRADE} < $max_cons) {$max_cons = $$element{GRADE};}
    }
    $ConsColorUnity = $max_cons / 4.5 * -1; 
    if ($max_cons !~ /^\-/){$ConsColorUnity = $max_cons;}        
        
#calculates the grades for each color

    my $NoLayers = 9;
    my @ColorLayers;
    for (my $i = 0; $i <= $NoLayers; $i++) {
        $ColorLayers[$i] = $max_cons + ($i * $ConsColorUnity);
    }
#gives the color to the interval

    my $Count = 0;

    foreach my $element (@$Output) {
        
        for (my $i = 0; $i <= $#ColorLayers; $i++){
            if( ($i==$#ColorLayers) and !exists $$element{INTERVALLOWCOLOR}){
                $$element{INTERVALLOWCOLOR} = 8;
            }
            elsif ($$element{INTERVALLOW} >= $ColorLayers[$i] and $$element{INTERVALLOW} < $ColorLayers[$i + 1]) {
                $$element{INTERVALLOWCOLOR} =$i;
            } 
            elsif ( ($$element{INTERVALLOW} < $ColorLayers[$i]) and !exists $$element{INTERVALLOWCOLOR}){
                $$element{INTERVALLOWCOLOR} = 0;
            } 
            if (($i == $#ColorLayers)  and !exists $$element{INTERVALHIGHCOLOR}){
                $$element{INTERVALHIGHCOLOR} = 8;
            }
            elsif ($$element{INTERVALHIGH} >= $ColorLayers[$i] and $$element{INTERVALHIGH} < $ColorLayers[$i + 1]) {
                $$element{INTERVALHIGHCOLOR} =$i;
            } 
            elsif ( ($$element{INTERVALHIGH} < $ColorLayers[$i]) and !exists $$element{INTERVALHIGHCOLOR}){
                $$element{INTERVALHIGHCOLOR} = 0;
            }
        } # END FOR
    } # END FOREACH
	
#give the color for each position based on the grades	
	# match the colors to the grades
    foreach my $element (@$Output){ 
        for (my $i = 0; $i <= $#ColorLayers; $i++) {
            if ($i == $#ColorLayers) {
                $$element{COLOR} = $ColorScale{$i-1};
            }
            elsif ($$element{GRADE} >= $ColorLayers[$i] && $$element{GRADE} < $ColorLayers[$i + 1]) {
                $$element{COLOR} = $ColorScale{$i};         
                last;
            }            
        }
        if ((($$element{INTERVALHIGHCOLOR}-$$element{INTERVALLOWCOLOR})>$bayesInterval ) or ($$element{MSA_NUM} <= 5)){
            $$element{ISD}=1;			
        }
        else{$$element{ISD}=0;}		
    }
	return ("OK");
} 
#************************************************************************************************
#printing the the gradesPE file
sub create_gradesPE{
	my ($Output,$ref_match, $ref_residue_freq, $no_isd_residue_color, $isd_residue_color, $gradesPE_file) = @_;
    my ($seq3d_grades_isd, $seq3d_grades);
    # open file
    unless (open PE, ">$gradesPE_file" ){
        return ("create_gradesPE : can't open '$gradesPE_file'","PANIC");}
    print PE "\t Amino Acid Conservation Scores\n";
    print PE "\t===============================\n\n";
    print PE "- POS: The position of the AA in the SEQRES derived sequence.\n";
    print PE "- SEQ: The SEQRES derived sequence in one letter code.\n";
    print PE "- 3LATOM: The ATOM derived sequence in three letter code, including the AA's positions as they appear in the PDB file and the chain identifier.\n";
    print PE "- SCORE: The normalized conservation scores.\n";
    print PE "- COLOR: The color scale representing the conservation scores (9 - conserved, 1 - variable).\n";
    print PE "- CONFIDENCE INTERVAL: When using the bayesian method for calculating rates, a confidence interval is assigned to each of the inferred evolutionary conservation scores.\n"; 
    print PE "- CONFIDENCE INTERVAL COLORS: When using the bayesian method for calculating rates. The color scale representing the lower and upper bounds of the confidence interval.\n"; 
    print PE "- MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.\n";
    print PE "- RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.\n\n";
    print PE " POS\t SEQ\t    3LATOM\tSCORE\t\tCOLOR\tCONFIDENCE INTERVAL\tCONFIDENCE INTERVAL COLORS\tMSA DATA\tRESIDUE VARIETY\n";
    print PE "    \t    \t        \t(normalized)\t        \t               \n";
    foreach my $elem (@$Output){
		my $pos = $$elem{POS};
		my $var = "";
		my $atom_3L=$$ref_match{$pos};
		my $score = $$elem{COLOR};
        
		printf (PE "%4d", $pos);
        printf (PE "\t%4s", "$$elem{SEQ}");	        
        printf (PE "\t%10s", $atom_3L);
        printf (PE "\t%6.3f", "$$elem{GRADE}");
        if($$elem{ISD}==1){
            printf (PE "\t\t%3d", "$$elem{COLOR}");
            printf (PE "%1s", "*");
        }
        else{printf (PE "\t\t%3d", "$$elem{COLOR}");}
        printf (PE "\t%6.3f", "$$elem{INTERVALLOW}");
        printf (PE "%1s", ",");
        printf (PE "%6.3f", "$$elem{INTERVALHIGH}");
        printf (PE "\t\t\t%5d", "$ColorScale{$$elem{INTERVALLOWCOLOR}}");
        printf (PE "%1s", ",");
        printf (PE "%1d\t\t", "$ColorScale{$$elem{INTERVALHIGHCOLOR}}");
        printf (PE "\t%8s", "$$elem{MSA_NUM}\/$$elem{MSA_DENUM}");
        for my $_aa (keys %{$ref_residue_freq->{($pos)}}){
            $var.= "$_aa,";
        }
        chop($var) if ($var =~ /,$/);
        print PE "\t$var\n";		
		# the amino-acid in that position, must be part of the residue variety in this column
		if ($var !~ /$$elem{SEQ}/){
			close PE;
			return ("create_gradesPE : in position $pos, the amino-acid ".$$elem{SEQ}." does not match the residue variety: $var.","PANIC");}
		#printing the residue to the rasmol script
        #assigning grades to $seq3d strings
        if($atom_3L !~ m/\-/){
            $atom_3L =~ m/(.+):/;
            $atom_3L = $1;
            if($score=~ m/(\d)/){
                my $color = $1;
                push @{ $no_isd_residue_color->[$color] }, $atom_3L;
                #if($score=~ m/\*/){
				if ($$elem{ISD}==1){
                    push @{ $isd_residue_color->[10] }, $atom_3L;
                    $seq3d_grades_isd.="0";
                }
                else{
                    push @{ $isd_residue_color->[$color] }, $atom_3L;                        
                    $seq3d_grades_isd.="$color";
                }
                $seq3d_grades.="$color";
            }
        }
        else{
            $seq3d_grades_isd.=".";
            $seq3d_grades.=".";
        }
    }
    print PE "\n\n*Below the confidence cut-off - The calculations for this site were performed on less than 6 non-gaped homologue sequences,\nor the confidence interval for the estimated score is equal to- or larger than- 4 color grades.\n";
    close(PE);
	return ("OK",$seq3d_grades_isd, $seq3d_grades);
}

#************************************************************************************************
#matches the position in the seqres/msa sequence to the position in the pdb
sub match_seqres_pdb{
	
	my ($seqres_atom_aln, $atom_position, $chain, $ref_fas2pdb) = @_;
	my (@seqres, @atoms, $length_of_seqres, $length_of_atom);
	my $pdbseq="";
	my $query_seq="";
	
	#creating arrays containing each sequences
	unless(open ALN, $seqres_atom_aln){return ("match_seqres_pdb : Could not open the file $seqres_atom_aln for reading $!", "PANIC");}
	flock ALN, LOCK_EX;
	while(<ALN>){
		if($_=~ m/^ATOM_\S+\s+(\S+)/){
			$pdbseq=$pdbseq.$1;
		}
		elsif($_=~ m/^SEQRES_\S+\s+(\S+)/){
			$query_seq=$query_seq.$1;
		}
		elsif($_=~ m/^MSA_\S+\s+(\S+)/){
			$query_seq=$query_seq.$1;
		}
	}
	flock ALN,LOCK_UN;
	close(ALN);
	#arrays that conatin the sequences from the pairwise for each sequneces (including gaps).
	#the sequnces start from position 1 in the array!
	$pdbseq="X".$pdbseq;
	$query_seq="X".$query_seq;
	@atoms=split("",$pdbseq);
	@seqres=split("",$query_seq);
	
	#creating the hash that matches the position in the ATOM fasta to its position
	#in the pdb file and also the fasta ATOM position to the correct residue
	unless(open MATCH, $atom_position){return ("match_seqres_pdb : Could not open the file $atom_position for reading $!", "PANIC");}
	
	flock MATCH, LOCK_EX;
	my %match_ATOM;
	my %res_ATOM;
	while(<MATCH>){
		$_=~ m/(\w{3})\s+(\d+)\s+(\S+)/;
		my $res=$1;
		my $fas_atom=$2;
		my $pdb_atom=$3;
		$match_ATOM{$fas_atom}=$res.$pdb_atom.":".$chain;
	}
	flock MATCH,LOCK_UN;
	close MATCH;
    #creating a hash in which the key is the position in the aln (i.e the
    #position in the @seqres) and the value is the correct position in the fas.
    my $pos_count=0;
    my %aln_fas_seqres;
    for(my $n=1;$n<@seqres+0;$n++){
#        if($seqres[$n] !~ m/\-/){ #HAIM
        if(($seqres[$n] !~ m/\-/) and ($seqres[$n] !~ m/X/)){ 
            $pos_count++;
            $aln_fas_seqres{$n}=$pos_count;
        }
    }
    $length_of_seqres = $pos_count;
    #creating a hash in which the key is the position in the aln (i.e the
    #position in the @atoms) and the value is the correct position in the fas.
    $pos_count=0;
    my %aln_fas_atoms;
    for(my $n=1;$n<@atoms+0;$n++){
#        if($atoms[$n] !~ m/-/){ #HAIM
	  if(($atoms[$n] !~ m/-/) and ($atoms[$n] !~ m/X/)){
            $pos_count++;
            $aln_fas_atoms{$n}=$pos_count;
        }
    }	
    $length_of_atom = $pos_count;    
    for(my $i=1;$i<@seqres+0;$i++){
		my $fas_pos=$aln_fas_seqres{$i};
        if($seqres[$i]!~m/-/){           
#            if($atoms[$i] eq "-"){ #HAIM
            if(($atoms[$i] eq "-") or ($atoms[$i] eq "X")){
                $ref_fas2pdb->{$fas_pos}="-";
            }
            else{
                 my $match=$match_ATOM{$aln_fas_atoms{$i}};
                $ref_fas2pdb->{$fas_pos}=$match;
            }
        }
    }
    return ("OK",$length_of_seqres, $length_of_atom);
}
#************************************************************************************************

#adds the 3LATOM colomn to the new gradesPE file. also creates a rasmol script.
# its input is a reference to a hash where key is the position in the SEQRES sequence, and the value is that 3LATOM according to the PDB sequence
sub add_pdb_gradesPE{
    my ($pre_gradesPE, $outPE, $ref_match, $no_isd_residue_color, $isd_residue_color, $ref_residue_freq)=@_;
	my ($seq3d_grades_isd, $seq3d_grades, $aa_position);
    my ($grade, $score, $high, $low, $gradel, $gradeh, $msa_m, $msa_d, $var);
    my $seq = "";
    unless (open GRADESPE, $pre_gradesPE){return ("add_pdb_gradesPE : Could not open the file $pre_gradesPE for reading $!", "PANIC");}
    unless(open PE, ">$outPE"){return ("add_pdb_gradesPE : Could not open the file $outPE for writing $!", "PANIC");}
    while(<GRADESPE>){
        if($_!~ m/\s+\d+\s+\w\s+atom_res/){
            print PE "$_";    
            next;
        }
        else{
            $_=~m/^\s+(\d+)\s+/;
            my $pos=$1;
            my $atom_3L=$$ref_match{$pos};
			$atom_3L=~ m/(\w{3})/;
			my $aa=$1;			
            $_=~ m/\s+(\d+)\s+(\w)\s+atom_res\s+(\S+)\s+(\S+)\s+(\S+),\s*(\S+)\s+(\S+),\s*(\S+)\s+(\d+)\/(\d+)/;
            $pos=$1; $seq = $2;
            $grade=$3; $score=$4; $high=$5;$low=$6;$gradel=$7;
            $gradeh=$8; $msa_m=$9; $msa_d=$10; $var="";
                
            printf (PE "%4d", "$pos");
            printf (PE "\t%4s", "$seq");
            printf (PE "\t%10s", "$atom_3L");
            printf (PE "\t%6.3f", "$grade");
            if($score=~m/\*/){
                printf (PE "\t\t%3d", "$score");
                printf (PE "%1s", "\*");
            }
            else{printf (PE "\t\t%3d", "$score");}
            printf (PE "\t%6.3f", "$high");
            printf (PE "%1s", ",");
            printf (PE "%6.3f", "$low");
            printf (PE "\t\t\t%5d", "$gradel");
            printf (PE "%1s", ",");
            printf (PE "%1d\t\t", "$gradeh");
            printf (PE "\t%8s", "$msa_m\/$msa_d");
            for my $_aa (keys %{$ref_residue_freq->{$pos}}){
                $var.= "$_aa,";
            }
            
            chop($var) if ($var =~ /,$/);
            print PE "\t$var\n";
            #printing the residue to the rasmol script
            #assigning grades to $seq3d strings
            if($atom_3L !~ m/\-/){
                $atom_3L =~ m/(.+):/;
                $atom_3L = $1;
                if($score=~ m/(\d)/){
                    my $color = $1;
                    push @{ $no_isd_residue_color->[$color] }, $atom_3L;
                    if($score=~ m/\*/){
                        push @{ $isd_residue_color->[10] }, $atom_3L;
                        $seq3d_grades_isd.="0";
                    }
                    else{
                        push @{ $isd_residue_color->[$color] }, $atom_3L;                        
                        $seq3d_grades_isd.="$color";
                    }
                    $seq3d_grades.="$color";
                }
            }
            else{
                $seq3d_grades_isd.=".";
                $seq3d_grades.=".";
            }
        }
    }
    close PE;
	close GRADESPE;
	return ("OK",$seq3d_grades_isd, $seq3d_grades);
}
#************************************************************************************************
#************************************************************************************************
#input: $input_pdb_file, $ref_return_arr
#output: ("OK"/"ERR",$header_line,$title_lines,$compnd_lines)

sub extract_data_from_pdb{
	my $input_pdb_file = shift;
	my @return_arr = ();
	unless (open PDB, $input_pdb_file){
		$return_arr[0] = "extract_data_from_pdb : Could not open the file $input_pdb_file $!";
		$return_arr[1] = "PANIC";
	}	
	else{
		flock PDB, 2;
		while (<PDB>){
			if(/^HEADER/)   {
				$_ =~ s/\s+$//;
				$return_arr[1] = $_;
			}
			elsif (/^TITLE\s+/){
				$_ =~ s/\s+$//;
				$_ =~ /^TITLE\s+\d*\s(.*)\s*/;
				$return_arr[2].=$1." ";
			}
			 elsif (/^COMPND\s+/) {
				$_ =~ s/\s+$//;
				$_ =~ /^COMPND\s+\d*\s(.*)/;
				$return_arr[3].=$1." ";
			}
			elsif (/^SOURCE/ || /^KEYWDS/ || /^AUTHOR/ || /^SEQRES/ || /^ATOM/) {# no nead to go over all the pdb
				last;
			}
		}
		flock PDB, 8;
		close PDB;
		$return_arr[0] = "OK";
	}
	return @return_arr;
}
#************************************************************************************************

# for each position, calculate the % for each residue which is found in the MSA in that position
# the values in the line might not sum to exactly 100%, as we round the value to be written with
# 2 digits after the '.'
sub print_precentage_text{
    my ($ref_residue_freq,$ref_position_totalAA, $out_file) = @_;    
    my @aa_arr = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","OTHER");
    my ($val, $aa_found, $total, @aa_in_position, $other, $other_val);
    unless (open OUT, ">".$out_file) {return ("Could not open the file $out_file for writing $!", 'PANIC');}
	print OUT "The table details the residue variety in % for each position in the query sequence.\nEach column shows the % for that amino-acid, found in position (\"pos\") in the MSA.\nIn case there are residues which are not a standard amino-acid, they are represented under column \"OTHER\"\n\n";
    print OUT "pos  |  ";
    print OUT "   -$_-  " foreach (@aa_arr);
    print OUT "\n" ;
    print OUT "---------" foreach (@aa_arr);
    print OUT "\n";
    # order the lines according to AA position in the sequence
    for my $position (sort {$a <=> $b} keys %$ref_residue_freq){
        $total=0; $aa_found = "no"; $other="";$other_val=0;
        printf OUT '%4s', $position;
        print OUT " | ";
        # the total number of animo acids found in the MSA for that position
        $total += $ref_position_totalAA->{$position};
        #For each position , sort the amino acids variety
        @aa_in_position = sort keys %{$ref_residue_freq->{$position}};
        my $j=0;
        my $aa = $aa_in_position[$j];
        $val = ($ref_residue_freq->{$position}->{$aa})/$total*100;
        # in order to print a table, we go by the sorted array
        my $i=0; 
		while ($i<@aa_arr){
			# if there are non standart aa, we calculate its total value seperately and add it under column OTHER
			while ($aa !~ /^[KRHDEYWSTFLIMCNQAVPG]$/ and $j < @aa_in_position){
				$other_val+=$val;
				$j+=1;
				if ($j < @aa_in_position){
					$aa = $aa_in_position[$j];
					$val = ($ref_residue_freq->{$position}->{$aa})/$total*100;
				}
			}
			if ($aa_arr[$i] eq 'OTHER' and $other_val!=0){
				$other = sprintf '%*2$.2f', $other_val;
				print OUT "     ".$other;
			}			
            elsif ($aa_arr[$i] eq $aa){				
                printf OUT '%*2$.2f', $val, 8;
                $j+=1;
                if ($j < @aa_in_position){
                    $aa = $aa_in_position[$j];
                    $val = ($ref_residue_freq->{$position}->{$aa})/$total*100;
					
                }
            }			
            else{
                print OUT "        ";
            }
			$i++;
        }
        print OUT "\n";
    }
    close OUT;
	return ("OK");
}
#************************************************************************************************
sub less_than{
    my $_val = shift;
    my $in = shift;
    if ($in>=0 and $in<$_val) {return 1;}
    else{return 0;}
}
#************************************************************************************************
# design the frequencies array
sub freq_array{
	my ($isd_residue_color, $no_isd_residue_color) = @_;
	my ($consurf_grade_freqs_isd, $consurf_grade_freqs);
# the insufficient data should be the first in this array
    $consurf_grade_freqs_isd = "Array(";
    $consurf_grade_freqs_isd .=$#{$isd_residue_color->[10]}+1;
    for (my $i=1; $i<10; $i++){
        $consurf_grade_freqs_isd .=",";
        $consurf_grade_freqs_isd .= $#{$isd_residue_color->[$i]}+1;    
    }
    $consurf_grade_freqs_isd .=")";
    $consurf_grade_freqs = "Array(0";
    print "\n";
    for (my $i=1; $i<10; $i++){
        $consurf_grade_freqs .=",";
        $consurf_grade_freqs .= $#{$no_isd_residue_color->[$i]}+1;    
    }
    $consurf_grade_freqs .=")";
	return ($consurf_grade_freqs_isd, $consurf_grade_freqs);
}

#************************************************************************************************
# Go over the PDB file and find all the chains that are identical to the given chain
sub find_identical_chains_on_PDB_File{
	my $PDBfile=shift;
	my $chain=shift;
	unless (open PDBFILE , "$PDBfile") {return "find_identical_chains_on_PDB_File : cannot open $PDBfile for reading $!\n";}
        my $IdenticalChainsLine = $chain;
        my %chainsHash=();
	while (<PDBFILE>)
	{
	      chomp($_);
    	      if (/^SEQRES\s+\d+\s+([A-Z0-9])\s+\d+\s+([A-Z\s]+)\s+/i)
		{
              	my $key = $1;
              	my  $line = $2;
              	$line =~ s/\s+//g;
              	$chainsHash{$key} .= $line;
		}
	}
	close PDBFILE;
	foreach my $key (keys  %chainsHash)
	{
    		if ($chainsHash{$key} eq $chainsHash{$chain} && ($key ne $chain))
		{
        	$IdenticalChainsLine .= $key;
    		}
	}
	return ("OK",$IdenticalChainsLine);
}

#************************************************************************************************
# create a file with all the ATOM records from a pdb.
# a list of residues, serial number and their position according to the pdb
# all the outpus are written to a hash (which is an input to the routine), so in order to extract output information:
# $ref_to_return_hash->{ERROR} : error
# $ref_to_return_hash->{NMR} : will be "yes" if it is NMR
# $ref_to_return_hash->{INFO} : information regarding nmr model
# $ref_to_return_hash->{WARNING} : in case a non-standard residue was found, will give its description
# $ref_to_return_hash->{AA_SEQ} : will contain the ATOM sequence in 1 letter amino-acids
sub create_atom_position_file{
	my ($pdb_filename ,$atom_position_filename, $chain_id, $ref_to_return_hash) = @_;
	my $nmr_model = "no";
	my $info_ret = "";
	my $warn_ret = "";
	my ($residue,$pos,$last_pos);
	my $first=1;
	my $fas="";
	my $pdb_line =0;
	
	unless (open PDB, $pdb_filename){
		$ref_to_return_hash->{ERROR} = "rasmol_gradesPE_and_pipe::create_atom_position_file : could not open the file $pdb_filename for reading $!";
		return;
	}
	unless (open CORR, ">$atom_position_filename"){
		$ref_to_return_hash->{ERROR} = "rasmol_gradesPE_and_pipe::create_atom_position_file : could not open the file $atom_position_filename for writing $!";
		return;
	}


    #going over the PDB file and extracting the chain sequence.
    while(<PDB>){
		$pdb_line++;
        # in case of NMR, we only take the first model's info
        if (/^MODEL\s+1\s*$/){
            $nmr_model = "yes";
            $ref_to_return_hash->{INFO} = "found NMR model 1.";
        }
        elsif (/^MODEL\s+(\d+)\s*$/ and $nmr_model eq "yes"){
            $ref_to_return_hash->{INFO} .= " stop reading PDB at NMR model $1";
            last;
        }
        #reaches an ATOM line, extracts the residue and prints out its position.
        #elsif($_=~ m/^ATOM.................$chain_id/){
		elsif($_=~ m/^ATOM/ and (substr ($_,21,1) eq $chain_id)){
			$residue=substr ($_,17,3);
            $pos=substr ($_,22,5);;
            if(exists($tr_aa{$residue})){
				if($first ==1){
					$fas=$fas.$tr_aa{$residue};   
                    $last_pos=$pos;
                    $fas_pos=1;
					print CORR "$residue\t$fas_pos\t$pos\n";
				}
				$first=0;
                if($pos ne $last_pos){
                    $fas=$fas.$tr_aa{$residue};   
                    $fas_pos++;
					print CORR "$residue\t$fas_pos\t$pos\n";
                }
                $last_pos=$pos;
            }
			else{
				$ref_to_return_hash->{WARNING}.= "line $pdb_line : residue $residue is not legal\n";
			}
		}
	}
	close PDB;
	close CORR;
	$ref_to_return_hash->{AA_SEQ} = $fas;
	$ref_to_return_hash->{NMR} = "yes" if ($nmr_model eq "yes");
}
#************************************************************************************************
#************************************************************************************************
sub read_gradesPE{
# the routine matches each position in the gradesPE file its grade. In case there was a grade mark with *, we put it in a seperate hash with the grade 0.
# the routine returns "yes" if a * was found and "no" otherwise
    my $gradesPE_file = shift;
    my $gradesPE_hash_ref = shift;
    my $gradesPE_0_hash_ref = shift;
	my $insufficient = "no";
    
    open GRADES, $gradesPE_file;
    while (<GRADES>){
        if (/^\s*\d+\s+\w/ ){
		    my @grades=split;            
		    $grades[2] =~ s/[a-z\:]//gi;
            if ($grades[4] =~/\d\*?/){
				# if it is insufficient color - we change its grade to 0, which will be read as light yellow
				if ($grades[4] =~/(\d)\*/){
					$gradesPE_hash_ref->{$grades[2]} = $grades[4];
                    $gradesPE_0_hash_ref->{$grades[2]} = 0;
                    $insufficient = "yes";
				}
				else{
					$gradesPE_hash_ref->{$grades[2]} = $grades[4];
                    $gradesPE_0_hash_ref->{$grades[2]} = $grades[4];
				}
            }
        }
    }
    close GRADES; 
    return $insufficient;
}


1;
