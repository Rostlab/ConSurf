#!/usr/bin/perl -w

package prepareMSA;

use Bio::SearchIO;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::AlignI;

use strict;

#----------------------------------
sub create_cd_hit_output{
#----------------------------------    
    my $working_dir = shift;
    my $input_file = shift;   # fasta file with sequences to cluster
    my $output_file = shift;  # final fasta file with chosen homolougs
    my $cutoff = shift;       # the decimal fraction of identity we allow between the homolougs
    my $cd_hit_dir = shift;   # directory from where we run cd-hit program
    my $ref_cd_hit_hash = shift; # hash that holds all the chosen homolougs sequences
    my $which_server = shift;
    
    my ($seq, $seq_name);
    
    ##################
    # running cd-hit #
    ##################
    my $cmd="./cd-hit -i $working_dir"."$input_file -o $working_dir"."$output_file -c $cutoff";
    chdir($cd_hit_dir);
    
    my $ans = `$cmd`;    
    
    unless ((-e $working_dir.$output_file)||(-z $working_dir.$output_file)){
        return ("sys", "parseFiles::create_cd_hit_output : CD-HIT produced no output!\n");    
    }
    #####################################
    # insert chosen homolougs to a hash #
    #####################################    
    my $cd_hit_result  = Bio::SeqIO->new('-file' => "$working_dir"."$output_file" , '-format' => 'Fasta'); 
    while ( my $seqObj = $cd_hit_result->next_seq() ) {
        $seq = $seqObj->seq();        
        $seq_name = $seqObj->primary_id();
	my $description = $seqObj->desc();
   #     $ref_cd_hit_hash->{$seq_name} = $seq;
   	$ref_cd_hit_hash->{$seq_name}{SEQ} = $seq;
	$ref_cd_hit_hash->{$seq_name}{DESCRIPTION} = $description;
    }
    $cd_hit_result->close();
    
    return "ok";
    
}

#----------------------------------
sub choose_homologoues_from_blast{
#----------------------------------    
    my $working_dir = shift;
    my $query = shift;  # query sequence, given in Fasta Format
        
    my $redundancyRate = shift;  #a number between 0-100 that will decide under which % we take homolougs    
    my $frag_overlap = shift;   # the maximum % of overlapping between 2 fragments of the same hit (should be given as decimal fraction
    my $min_length_percent = shift;  # the minimum percent a homoloug should have ( given as decimal fraction)
    my $min_num_of_homolougs = shift; # the minimum number of homolougs we must collect
    
    my $blast_output = shift;  # output of blast run, to read from, full path
    my $fasta_output = shift;  # the final homolouges file, to write to
    my $rejected_seqs = shift; # file that will have the rejected sequences and the reason for it
    
    my $ref_blast_hash = shift; # a hash to hold all the sequences names and their e-values
            
    my ($query_seq_name, $query_seq_length, $query_AAseq, $min_length, $seq_details, $s_name, $AAseq, $s_eval, $s_beg, $s_end, $s_ident,$s_description,$seq_description);
    my $seq_name_exists = "no";
    my $ans = "no";
    
    my %sequences = ();
    # hash of sequences names. each sequence name (unique) points to array of hashes.
    # {seq_name} => [{e_val => <e_val1>, AAseq => <AAseq1>, beg => <S1_beg>, end => <S1_end>},...,{e_val => <e_valn>, AAseq => <AAseqn>, beg => <Sn_beg>, end => <Sn_end>}]
    # seq_name : the name,
    # S1_beg : the index for the begining of the match for S1,
    # S1_end : the index to the end of the match for S1
    # e_val1 : the e-value. 
    # AAseq1 : the sequence itself, as read from the "Subjct" line of the HSP    
    
    unless (open QUERY, $query){
        return ("sys", "parseFiles::choose_homologoues_from_blast : can't open file $query for reading\n");
    }
    ################################
    # Extracting Query information #
    ################################    

    while (<QUERY>){
        if ($_ !~ m/>/){
            $_=~m/(\S+)/;
            $query_seq_length += length($1);
        }
    }
    close QUERY;    
    
    ######################################################
    # defining the minimum length a homoloug should have #
    ######################################################
    $min_length = $query_seq_length*$min_length_percent; #60% the query's length
    
    ##################################################
    # Reading blast output and collect the homolougs #
    ##################################################
    
    unless (open OUT_REJECT, ">".$working_dir.$rejected_seqs){
        {return ("sys", "parseFiles::choose_homologoues_from_blast : can't open file $rejected_seqs for writing\n");}
    }    
    my $searchio = new Bio::SearchIO(-format => 'blast',
                                     -file   => $blast_output);
    while( my $result = $searchio->next_result ) {  # $result is a Bio::Search::Result::ResultI object
        while( my $hit = $result->next_hit ) { # $hit is a Bio::Search::Hit::HitI object or undef if there are no more
            $s_description=$hit->description();
	    $s_name= $hit->name();
            if ($s_name =~ m/.+\|(\S+\|\S+_\S+)/) {$s_name = $1;}
            if ($s_name =~ m/.+\|(\S+_\S+)\|(\S+)/) {$s_name = $2.'|'.$1;}
            elsif ($s_name =~ m/(\S+\|\S+_\S+)/) {$s_name = $1;}
            elsif ($s_name =~ m/(\S+_\S+)\|(\S+)/) {$s_name = $2.'|'.$1;}
            $seq_name_exists = "no";
            while( my $hsp = $hit->next_hsp ) { #hsp is the next available High Scoring Pair, Bio::Search::HSP::HSPI object or null if finished
                # extracting relevant details from the fragment
                ($s_beg, $s_end) = $hsp->range("sbjct");
                $AAseq = $hsp->hit_string();
                $AAseq =~ s/-//g;
                $s_eval = $hsp->evalue();
                $s_eval =~ s/,//g;
                if ($s_eval =~ m/^e/) {$s_eval = "1".$s_eval;}
                $s_ident = $hsp->percent_identity();                
                # deciding if we take the fragment
                # in case there is already a fragemnt with the same name, we do another validity test for overlapping
                if (exists $sequences{$s_name}){
                    $seq_name_exists = "yes";
                    $seq_details = $sequences{$s_name};
                    $ans = &check_if_seq_valid($redundancyRate, $min_length, $s_ident, $AAseq, $s_name);
                    if ($ans eq "yes"){
                        {$ans = &check_if_no_overlap($frag_overlap, $seq_details, $s_beg, $s_end);}
                    }                    
                }
                else {
                    $ans = &check_if_seq_valid($redundancyRate, $min_length, $s_ident, $AAseq, $s_name);
                }
                # after taking the info, check if the currecnt sequence is valid. If so - insert it to the hash
                if ($ans eq "yes"){
                    # in case there is more than one fragment for this seq_name: add another hash to the details array
                    if ($seq_name_exists eq "yes"){
#                        push @$seq_details, {e_val => $s_eval, AAseq => $AAseq, beg => $s_beg, end => $s_end}; #HAIM
			push @$seq_details, {e_val => $s_eval, AAseq => $AAseq, beg => $s_beg, end => $s_end, description => $s_description};
                    }
                    # in case it is the first fragment for this seq_name: insert a details array as a value for this seq_name key
                    else{
#                        $sequences{$s_name} = [{e_val => $s_eval, AAseq => $AAseq, beg => $s_beg, end => $s_end}]; #HAIM
			$sequences{$s_name} = [{e_val => $s_eval, AAseq => $AAseq, beg => $s_beg, end => $s_end, description => $s_description}];
                    }
                }
                else {print OUT_REJECT "Fragment $s_name"."_$s_beg"."_$s_end rejected: $ans\n";}
            }
        }
    }
    $searchio->close();
    close OUT_REJECT;
    
    ##############################################################################
    # Print the selected homolougs to a file and insert the e-value info to hash #
    ##############################################################################
    
    my $i;
    my $seq_frag_name;
    
    unless (open OUT ,">".$working_dir.$fasta_output)
        {return ("sys", "parseFiles::choose_homologoues_from_blast : can't open file $fasta_output for writing\n");}
    while (($s_name, $seq_details) = each (%sequences)){        
        for ($i=0; $i<=$#$seq_details; $i++){
            $seq_frag_name = "$s_name"."_".$seq_details->[$i]{beg}."_".$seq_details->[$i]{end};
	    $seq_description= $seq_details->[$i]{description};
            #print OUT ">$seq_frag_name\n"; #HAIM
            print OUT ">$seq_frag_name | $seq_description\n";
	    print OUT $seq_details->[$i]{AAseq}."\n";
            $ref_blast_hash->{$seq_frag_name} = $seq_details->[$i]{e_val};
        }
    }
    close OUT;
    
    #####################################################
    # Check that the number of homolougs found is legal #
    #####################################################
    
    my $final_num_homolougs = (keys %$ref_blast_hash);
    my $ret = "";
    if ($final_num_homolougs == 1){
        $ret .= "only <a href=$fasta_output>one unique sequence</a> ";
    }
    elsif ($final_num_homolougs <= $min_num_of_homolougs){
        $ret .= "only <a href=$fasta_output>$final_num_homolougs unique sequences</a> ";
    }
    else{
        $ret = "ok";
    }
    if ($ret !~ /^ok/){
        my $blast_link;
        $blast_link = $1 if ($blast_output =~ /$working_dir(.+)/);
        $ret .= "were chosen from <a href=$blast_link>PSI-BLAST output</a>. (<a href =$rejected_seqs>Click here</a> for the list of the sequences which were not chosen from PSI-BLAST output).<br />The minimal number of sequences required for the calculation is $min_num_of_homolougs.<br>";
        return ("user", $ret);
    }
    return $ret;
}

#----------------------------------
sub check_if_seq_valid{
# checks if a sequence read from blast is valid according to some parameters.
# Found valid - returns "yes"
#----------------------------------
#   GLOBAL INPUT
    my $redundancyRate = shift; # The maximum % of identity
    my $min_length = shift;     # The minimum length a sequence should have
    
#   HOMOLOUG SPECIFIC INPUT
    my $ident_percent = shift;  # percentage identity
    my $aaSeq = shift;          # the Amino Acid sequence
    my $seqName = shift;        # the sequence identifier
    
    my $seq_length = length($aaSeq);
    my $ans;
    
    #parameter1: the sequence identity is less than the idnentity percent that was defined
    if ($ident_percent >= $redundancyRate){ $ans = "identity percent $ident_percent is too big";}
    
    #parameter2: the sequnece length is greater than the minimum sequence length
    elsif ($seq_length<$min_length){$ans = "the sequence length $seq_length is too short" ;}
    
    #parameter2: the sequnece letters should be legal to rate4site
    elsif ($aaSeq!~ m/^[ACDEFGHIKLMNPQRSTVWYBZX]+$/){$ans = "illegal character was found in sequence: $seqName";}
    
    else {$ans="yes"};
    
    return $ans;      
}

#----------------------------------
sub check_if_no_overlap{

# tests:
# 1. seq2 is inside seq1 or seq1 is inside seq2 : exit (since we assume, according to the order HSP are written,  that the one which was already chosen has a better e-val)
# 2. seq1 starts before seq2
# 3. seq2 starts before seq1
# -> calculate the overlapping fragment
# since seq1 has better e-val, we check if the overlapping fregment is more than 10% ($max_overlap) the length of seq 1
# yes: exit (return "no"). no: check that there is no overlap for all other fragments
# 
#----------------------------------
#   GLOBAL INPUT    
    my $max_overlap = shift;   # the maximum % of overlapping between 2 fragments of the same hit (should be given as decimal fraction) 
    my $ref_seq_details = shift;  # reference to the sequence details array of hash which holds the info for that specific fragment

#   FRAGMENT (SUBJCT) IN QUESTION SPECIFIC INPUT
    my $s_bgn = shift;          # first amino_acid index
    my $s_end = shift;          # last amino acid sequence

    my ($fragment_beg, $fragment_end, $fragment_length, $overlap_length);
    my $ans = "check_if_no_overlap : no ans was picked";
    
    my $i=0;
    while ($i<=$#$ref_seq_details){ # read data from each array
        $fragment_beg = $ref_seq_details->[$i]{beg};
        $fragment_end = $ref_seq_details->[$i]{end};
        $fragment_length = $fragment_end-$fragment_beg+1;
        
        # fragment is inside subjct or subjct is inside fragment
        if ($s_bgn <= $fragment_beg && $s_end >= $fragment_end) {
            return "previous fragment found $fragment_beg"."_$fragment_end is fully inside new fragment";
        }
        elsif ($s_bgn >= $fragment_beg && $s_end <= $fragment_end){
            return "new fragment is fully inside previous fragment found $fragment_beg"."_$fragment_end";
        }
                        
        # fragment begins before subjct
        elsif($fragment_end<$s_end && $fragment_end>$s_bgn){
            $overlap_length = $fragment_end - $s_bgn + 1;
            if ($overlap_length > $fragment_length*$max_overlap) {
                return "overlap length of fragment is ". $overlap_length." which is greater than maximum overlap: ".($fragment_length*$max_overlap);
            }
            # when the fragment might be a good match, we can only insert it if it did not match to all the fragments
            elsif($i==$#$ref_seq_details){
                $ans = "yes";
            }
        }
        # fragment begins after subjct
        elsif($fragment_beg>$s_bgn &&  $fragment_beg < $s_end){
            $overlap_length = $s_end - $fragment_beg + 1;
            if ($overlap_length > $fragment_length*$max_overlap) {
                return "overlap length of fragment is ". $overlap_length." which is greater than maximum overlap: ".($fragment_length*$max_overlap);
                }
            # when the fragment might be a good match, we can only insert it if it did not match to all the fragments
            elsif($i==$#$ref_seq_details){
                $ans = "yes";
            }
        }
        #no overlap
        elsif($fragment_beg>=$s_end || $fragment_end<=$s_bgn){
            if($i==$#$ref_seq_details){
                $ans = "yes";
            }
        }
        $i++;
    }
    return $ans;
}
1;

#----------------------------------
sub get_blast_round{
#   reads the blast file, look for the last round found.
#   also check if converged
#
# Params:
#   blast_output_filename - name of the blast file to process
#
# Returns:
#  If error:("err", error_description)
#  If no hits found : ("no_hits")
#  else: (last_round_number, 0 or 1 if found_converged)
#----------------------------------
    my $blast_output_filename = shift;    
    my $last_round_number = 0;
    my $found_converged = 0;
    my $ret = "";
# search file for last round

    unless (open(BLAST, "$blast_output_filename")){
		return("err", "cannot open file '$blast_output_filename' : $!");
	}

  	while (<BLAST>){
        if (/\*+ No hits found/){
            return("no_hits");
        }
		if(/Results from round (\d+)/){
        	$last_round_number = $1;
        }
        if(/CONVERGED!/){
        	# found last round
        	$found_converged = 1;
        	last;
        }
    }
	close BLAST;
    return ($last_round_number, $found_converged);
}
#----------------------------------
sub print_blast_according_to_round{
# prints the content of the file from blast, starting from $round_number
# Returns:
#  If error : ("err", error_description)
#  Else:  ("no_err", 0 if round not found; 1 otherwise)
#----------------------------------
    my $blast_output_filename = shift;
    my $round_number = shift; # since blast can return less rounds than requested, the blast file will be searched for last round
    my $output_filename = shift;
    my $round_found = 0;
    
    # open files
    unless (open(BLAST, "$blast_output_filename")){
    	return("err", "Can't open file '$blast_output_filename': $!");
    }
    unless (open(NEW_FILE, ">$output_filename")){
    	return("err", "Can't open file '$output_filename': $!");
    }
    # read blast file
    my $output_data = 1;
    my $line;
    while($line = <BLAST>){
        if($line =~ /Results from round (\d+)/){
            if($round_number == $1){
                $output_data = 1;
                $round_found = 1;
            }
            else{
                $output_data = 0;
            }            
        }
        if ($output_data == 1){
            print NEW_FILE $line;
        }
    }
    
    # close files
    close NEW_FILE;
    close BLAST;
    
    return("no_err",$round_found);
    
}
#----------------------------------
sub sort_sequences_from_eval{
# input: reference to hash which holds sequences and their evalue
# optional: query sequence input file (to extract the query sequence and add it)
# output: a file with the top-evalued sequences out of the hash
#----------------------------------
    my $ref_to_seqs_hash = shift;
    my $ref_to_cd_hash = shift;
    my $max_num_homologs = shift;
    my $output_file = shift;
    my $input_file = shift;
    my ($s_name, $s_aa_sq, $query_name, $counter);
    my $query_AAseq = "";
    
    # open original seq fasta file
    if (defined($input_file) and -e $input_file){
        unless (open QUERY, $input_file) {return ("err","Can't open '$input_file' for reading $!");}
        
        # get query name & amino acid seq from file
        while (<QUERY>)    {
            chomp;
            if (/^>(.+)/){
                $query_name = $1;
            }
            else{
                $query_AAseq.= $_;
            }
        }
        close QUERY;
    }

    # choose best homologs and create final file final_homolougs_filename
    $counter = 1;
    unless(open FINAL, ">>".$output_file) {return ("err","Can't open '$output_file' for writing $!");}

    # write query details
    print FINAL ">$query_name\n"."$query_AAseq\n" if ($query_AAseq ne "");
    
    my %blast_hash = %$ref_to_seqs_hash;
    my %cd_hit_hash = %$ref_to_cd_hash;
    # write homologs
    foreach $s_name (sort { $blast_hash{$a} <=> $blast_hash{$b} } keys %blast_hash ){
        # quit if reached max number of homologs    
        if ($counter > $max_num_homologs){
            last;
        }
        
        # write next homolog
        #if (defined($cd_hit_hash{$s_name})){ #HAIM
	if (defined($cd_hit_hash{$s_name}{SEQ})){
            #$s_aa_sq = $cd_hit_hash{$s_name};
	    $s_aa_sq = $cd_hit_hash{$s_name}{SEQ};
	    my $s_description = $cd_hit_hash{$s_name}{DESCRIPTION};
            print FINAL ">$s_name $s_description\n".$s_aa_sq."\n";
            $counter++;
        }    
    }
    close FINAL;
    return ("ok");
}

#----------------------------------
#----------------------------------