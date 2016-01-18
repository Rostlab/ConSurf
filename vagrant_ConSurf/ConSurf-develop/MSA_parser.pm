#!/usr/bin/perl -w

package MSA_parser;

use Bio::SearchIO;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::AlignI;
use Fcntl ':flock'; # import LOCK_* constants

use strict;

#----------------------------------
sub msa_format_extract {
#   input: MSA, query sequence name
#   output: query sequence
# The MSA is first converted to clustalw format. than it is parsed and searched for the query sequence name
#----------------------------------
    my ($working_dir,$original_msa_file,$clustalw_msa_file,$clustalw_out,$query_seq_name,$from_which_server,$report) = @_;
    $query_seq_name =~ s/\(|\)/_/g; # replace () with _, this is a conversion done by clustalw
    open OUT, ">$working_dir/$report";
    my @ret = ();
    my $cmd = "clustalw -convert -infile=$original_msa_file -outfile=$clustalw_msa_file -output=msf > $clustalw_out";    
    
    if ($from_which_server ne "biocluster"){
        $cmd = 'ssh bioseq@biocluster \'cd '.$working_dir."; $cmd'";
    }
    else{
        chdir $working_dir;
        print OUT "run from $working_dir\n";
    }
    print OUT "running $cmd\n";
    `$cmd`;
    # ERROR if no output
    if (!-e "$working_dir/$clustalw_msa_file" or (-e "$working_dir/$clustalw_msa_file" and -z "$working_dir/$clustalw_msa_file")){
        print OUT "file $working_dir/$clustalw_msa_file not found or zero";
        if (!-e "$working_dir/$clustalw_out" or -z "$working_dir/$clustalw_out"){
            @ret = ("err", "no_msa");
        }
        else{
            @ret = ("err", "read $clustalw_out");
        }        
    }
    # parse the new clustalw file, extract query sequence, eliminate gaps
    else{
        my $aln_file_obj = Bio::AlignIO->new(-file => "< $working_dir/$clustalw_msa_file", -format => "clustalw");
        my $aln_obj = $aln_file_obj->next_aln; # $aln_obj is a Bio::Align::AlignI compliant object
        my $num_of_seq_in_msa = 0;
        my $query_sequence = "";
        foreach my $seq_obj($aln_obj->each_seq){  # $seq_obj is an array of Seq objects 
            $num_of_seq_in_msa++;
            if($seq_obj->id eq $query_seq_name){
                $query_sequence = $seq_obj->seq();
            }
        }
        if ($num_of_seq_in_msa==0){
            @ret = ("err", "no_sequences_in_msa");
        }
        elsif($query_sequence eq ""){
            @ret = ("err", "no_sequence_found");
        }
        else{
            $query_sequence =~ s/\-//g;
            @ret = ($num_of_seq_in_msa, $query_sequence);
        }
    }
    close OUT;
   return @ret;
}

#----------------------------------
sub determine_msa_format{
# reads the first character in the MSA and returns the msa format
#----------------------------------    
    my $msa_input = shift;   
    my @ret = ();
    my $first_char = "";
    my $rest = "";
    my $found_pir = 0;
    open MSA, $msa_input or return ("err", "cannot open the file $msa_input $!");
    LINE: while (<MSA>){
        my $line = $_;
        if ($line =~ /^\s*(\S)(.+)?/){
            $first_char = $1;            
            $rest = $2 if (defined $2);
            if ($first_char eq 'M'){
                if ($rest =~ /^SF:/){
                    @ret = ("format", "gcg");
                }
            }
            elsif ($first_char eq '>'){
                if ((defined $rest) and $rest =~ /^P1\;/){
                    while (<MSA>){
                        unless (/^\s*>\s*$/){                        
                            if (/\*\s*$/){
                                $found_pir = 1;
                                @ret = ("format", "pir");
                            }
                        }
                        else{
                            @ret = ("format", "fasta");
                        }
                    }
                }
                else{
                    @ret = ("format", "fasta");
                }
            }
            elsif($line =~ /\s*MSF:\s+/){
                @ret = ("format", "gcg");
            }
            last LINE;
        }
        # MSF is not nessecrialy the first phrase, but it hould be in the first line
        
    }
    close MSA;
    unless (exists $ret[0]){
        if ($first_char eq '#'){
            @ret = ("format", "nexus");
        }    
        elsif($first_char eq 'C'){
            @ret = ("format", "clustalw");
        }
        elsif($first_char eq 'P'){  # MSF or GCG format: http://www.bioperl.org/wiki/MSF_multiple_alignment_format
            @ret = ("format", "gcg");
        }
        #elsif($first_char eq '%'){
        #    @ret = ("format", "gde");
        #}
        elsif (!defined($ret[0])){
            @ret = ("err", "unknown_format");
        }
    }
    return @ret;
    
    ## RSF : http://www.hku.edu/bruhk/gcgdoc/figure/using_sequences_8.gif
    ## Rich Sequence Format
    #elsif($first_char eq '!'){
    #    @ret = ("format", "rsf");
    #}
    # GDE : http://ubik.microbiol.washington.edu/HMA/html/gde_format.html
    #elsif($first_char eq ';'){
    #    @ret = ("format", "mase");
    #}    
}
#----------------------------------
sub check_msa_licit{    
# read MSA according to mode, look for ilegal characers
# supported format by bio:seqIO are: fasta, gcg, pir
# supported format by Bio::AlignIO are: clustalw, fasta, msf (gcg), nexus,
#
# If illegal chars were found, returns an array with 3 cells:
# ('user_error',"SEQ_NAME: <the sequence name where the error was found>", "IRR_CHAR: <a string with the iregular chars>")
# If an exception found - the SeqIO and AlighIO may through an exception, returns an array with 2 cells:
# ('user_error',"exception")
# Any other falut:
# ('user_error',"could not read msa")
#----------------------------------
    my ($msa_file,$msa_format) = @_;
    my $ans = "";
    my ($msa_fh, $seq, $seq_name, @ret);
    my $read_seq = 0;
    
    if ($msa_format eq "fasta" or $msa_format  eq "gcg" or $msa_format eq "pir"){
        $msa_fh = Bio::SeqIO->new(-file => $msa_file ,
                                -format => $msa_format);
        unless (defined $msa_fh) {return ('user_error',"could not read msa");}
        # this instruction may throw an exception if parameters are incorrect, for example - if the bioPerl identified too many unvalid characters. In order to catch the exception, we put it in a "catch" block
        eval{
            while ( $seq = $msa_fh->next_seq() ) { # returns a Bio::Seq sequence object
                $read_seq = 1;
                $ans = &check_sequence_licit($seq->seq());
                last if $ans ne "";
                
            }
        };
        if( $@ ) { #Caught exception
            return ('user_error',"exception");
        }
    }
    elsif($msa_format eq "clustalw" or $msa_format eq "nexus"){
        my $msa_fh = Bio::AlignIO->new(-file => $msa_file ,  
                                -format => $msa_format);
        unless (defined $msa_fh) {return ('user_error',"could not read msa");}
        # this instruction may throw an exception if parameters are incorrect, for example - if the bioPerl identified too many unvalid characters. In order to catch the exception, we put it in a "catch" block
        eval{
            while ( my $aln = $msa_fh->next_aln() ) { # Returns : a Bio::Align::AlignI
                foreach $seq ($aln->each_seq() ){ # Returns : an array of Seq objects
                    $read_seq = 1;
                    $ans = &check_sequence_licit($seq->seq());
                    last if $ans ne "";
                }
            }
        };
        if( $@ ) {#Caught exception
            return ('user_error',"exception");
        }
    }
    else{
        return ('user_error',"could not read msa");
    }
    if ($ans ne ""){
        @ret = ('user_error',"SEQ_NAME: $seq_name", "IRR_CHAR: $ans");
    }
    elsif($read_seq ==0){
        @ret = ('user_error',"could not read msa");
    }
    else{
        @ret = ("OK");
    }
    return @ret;
}
#----------------------------------
sub check_sequence_licit{
# look for iregular chars in AA sequence
#----------------------------------
    my $sequence = shift;
    my $no_regular_format_char = '';
    my $no_regular_note = '';
    
    $sequence =~ s/\s*$//;
    my @iregularchars   = split('',$sequence);
    foreach my $char(@iregularchars){
        if ($char !~ /[ACDEFGHIKLMNPQRSTVWXY\-]/i){
            if ($char =~ /([BJOUZ])/){
                $no_regular_note .=  " \"$1\", " if $no_regular_note !~ /$1/;
            }
            elsif ($char =~ /\*/){
                $no_regular_format_char .= qq( \"\*\", ) if $no_regular_format_char !~ /\*/;
            }
            elsif ($char =~ /\$/){
               $no_regular_format_char .= qq( \"\$\", ) if $no_regular_format_char !~ /\$/;
            }
            elsif ($char =~ /\./){
              $no_regular_format_char .= qq( \"\.\", ) if $no_regular_format_char !~ /\./;
            }
            elsif ($char =~ /\?/){
              $no_regular_format_char .= qq( \"\?\", ) if $no_regular_format_char !~ /\?/;                
            }
            elsif ($char =~ /\|/){
                $no_regular_format_char .= qq( \"\|\", ) if $no_regular_format_char !~ /\|/;
            }
            elsif ($char eq '\\'){ # it is imposible to have a RE trying to match the char \
                $no_regular_format_char = "\"\\\", ";
            }
            elsif ($char =~ /(.)/){
                $no_regular_format_char .= "\"$1\", " if $no_regular_format_char !~ /$1/;
                
            }            
        }
    }
    $no_regular_format_char .= $no_regular_note;
  
    if ($no_regular_format_char){
        $no_regular_format_char =~ s/,\s*$//;
    }
    return $no_regular_format_char;
}

#----------------------------------
sub get_info_from_msa{    
# read the MSA, returns a hash with all the sequences id and their corresponding sequences
# input : 1. full path for MSA, 2. MSA format, 3. reference to hash
# OPTIONAL INPUT:
# if you wish to also create a new MSA file in fasta format:
# 4. full path to an empty file 
# outputs: ERROR:# 
# ('user_error',"exception") - If an exception found - the SeqIO and AlighIO may through an exception
# ('user_error',"could not read msa") - the MSA was not opened for some reason
# ('user_error',"duplicity <sequence id>" - in case the same sequence name was found more than once
# ('user_error',"no seq id") - a sequences id is missing
# read successfully :
# ('OK')
#----------------------------------
    my ($msa_file,$msa_format, $ref_hash_seq_names, $msa_out_file) = @_;
    my $ans = "";
    my ($msa_fh, $seq, $seq_name, @ret);
        
    if ($msa_format eq "fasta" or $msa_format  eq "gcg" or $msa_format eq "pir"){
        $msa_fh = Bio::SeqIO->new(-file => $msa_file ,
                                -format => $msa_format);        
        unless (defined $msa_fh) {return ('user_error',"could not read msa");}
        # this instruction may throw an exception if parameters are incorrect, for example - if the bioPerl identified too many unvalid characters. In order to catch the exception, we put it in a "catch" block
        eval{
            while ( $seq = $msa_fh->next_seq() ) { # returns a Bio::Seq sequence object
                if ($seq->id() =~ /^\s*$/){
                    @ret = ('user_error',"no seq id");
                    last;
                }
                elsif (exists $ref_hash_seq_names->{$seq->id()}){
                    @ret = ('user_error',"duplicity ".$seq->id());
                    last;
                }
                else{
                    $ref_hash_seq_names->{$seq->id()} = $seq->seq();
                }
            }
        };
        if( $@ ) { #Caught exception
            @ret =  ('user_error',"exception");
        }
        $msa_fh->close();
    }
    elsif($msa_format eq "clustalw" or $msa_format eq "nexus"){
        $msa_fh = Bio::AlignIO->new(-file => $msa_file ,  
                                -format => $msa_format);
        $msa_fh->verbose(2); # in order to catch warning which are sent from method "warn()", we need to change the verbose value. otherwise, the warnings are silent
        unless (defined $msa_fh) {return ('user_error',"could not read msa");}
        # this instruction may throw an exception if parameters are incorrect, for example - if the bioPerl identified too many unvalid characters. In order to catch the exception, we put it in a "catch" block
        eval{
            while ( my $aln = $msa_fh->next_aln() ) { # Returns : a Bio::Align::AlignI
                foreach $seq ($aln->each_seq() ){ # Returns : an array of Seq objects                    
                    if ($seq->id() =~ /^\s*$/){
                        @ret = ('user_error',"no seq id");
                        last;
                    }
                    elsif (exists $ref_hash_seq_names->{$seq->id()}){
                        @ret = ('user_error',"duplicity ".$seq->id());
                        last;
                    }
                    else{
                        $ref_hash_seq_names->{$seq->id()} = $seq->seq();
                    }
                }
            }
        };
        if( $@ ) {#Caught exception
            # since the msa is not read "by blocks" but by sequences, if the same sequence id appears more than once - the AlignIO object throughs a warning with the duplicity and continues, but concating the 2 sequences. In this was, we catch this warning and report it to the user (instead of letting the file handler continue reading the file)
            if ($@ =~ /MSG: Duplicate sequence : (.+)/){
                @ret = ('user_error',"duplicity ".$1);
            }
            else{
                @ret = ('user_error',"exception");
            }
        }
        $msa_fh->verbose(0); # when we are done, setting the value of verbosity back to 0
        $msa_fh->close();
    }
    else{
        return ('user_error',"could not read msa");
    }
    if (defined $msa_out_file){
        open MSA_OUT, ">$msa_out_file";
        foreach my $key (keys %$ref_hash_seq_names){
            print MSA_OUT ">$key\n".$ref_hash_seq_names->{$key}."\n";
        }
        close MSA_OUT;
    }
    close OUT;
    
    unless (exists $ret[0]){@ret = ("OK")};
    return @ret;
}
#----------------------------------
sub check_validity_tree_file {
# check the validity of the newick format of the uploaded tree
# input: path to a tree file, reference to error hash
#
# if find error, mark it through the error hash:
# the marked items are: 'bootstrap' - if found bootstrap value,
#   'internal_nodes' - if found internal nodes
#   'left_right' - if the brackets are not equal in number
#   'noRegularFormatChar' - a list of all the non standard characters
#----------------------------------
    my ($treeFile, $ref_to_err) = @_;
    my $lineCounter=0;
    my $rightBrackets=0;
    my $leftBrackets=0;
    my (@lineArr, $line, $noRegularFormatChar, $treeFileOneLine);
    my $errorBool = 0;
    my $internal_nodes = 0;
    my $bootstrap = 0;

    my $read_right_bracket = "no";
    
    # in case the uploaded tree file contains more than one line -
    # read the tree and rewrite it
    unless (open TREEFILE, $treeFile) {$ref_to_err->{error} = "could not read $treeFile $!";}
    while (<TREEFILE>) {
       $line = $_;
       chomp($line);
       $treeFileOneLine .= $line;
       $lineCounter++;
    }
    close TREEFILE;   
    $line =  $treeFileOneLine;
    # add a semi-colon if missing
    if ($line !~ m/;\s*$/){
        $line =~ s/\s*$//;
        $line.=';';
    }   
    if ( $lineCounter>1) {
        unless (open TREEFILE, ">$treeFile") {$ref_to_err->{error} = "could not write to file $treeFile $!";}
        print TREEFILE  $line; 
        close TREEFILE;
    }
    # legal tree: same number of left and right brackets, no irregular chars
    @lineArr=split(//,$line); 
    foreach my $chain(@lineArr) {
        if ($chain eq '(')  {
           $leftBrackets++;
           $read_right_bracket = "no";
        }
        elsif ($chain eq ')') {
           $rightBrackets++;
           $read_right_bracket = "yes";
        }	
        elsif ($chain =~ /([\!|\@|\#|\$|\^|\&|\*|\~|\`|\{|\}|\'|\?|\\|\/|\<|\>])/){
           $noRegularFormatChar .= " '$1', " if $noRegularFormatChar !~ /\Q$1\E/;
           $read_right_bracket = "no";
        }
        # if right after a right Bracket we read a character which is not legal (ie: , : ;) we output a message to the user, since we don't handle bootstrap values or internal node names
        else{
            if($read_right_bracket eq "yes"){
                if($chain =~ /\d/){
                    $ref_to_err->{bootstrap} = 1;
                }
                elsif($chain !~ /[,|:|;]/){
                    $ref_to_err->{internal_nodes} = 1;
                }
            }
           $read_right_bracket = "no";
        }
    }
    if ($leftBrackets ne $rightBrackets) {
        $ref_to_err->{left_right}=1;
    } 
    if ($noRegularFormatChar =~ /.+/)  {
       $noRegularFormatChar =~ s/\,\s$//;
       $ref_to_err->{noRegularFormatChar} = $noRegularFormatChar;
    } 
}
#----------------------------------
sub read_residue_variety{
# the routine creates 2 hashes, according to information it reads from a MSA
# for each position in the query sequence :
# 1. it collects all the residues that are aligned to this positions. it also counts the number of time each residue appeared in that position in the MSA
# 2. it counts the total number of residues which aligned to this position
#----------------------------------
	my ($msa_file, $msa_ref_sequence, $msa_format, $ref_residue_frequency, $ref_position_totalAA) = @_;
	my $line;
	my $position = 1;
	my $num_of_seqs = 0;
	my $flag=0;
    my ($msa_fh, $seq, $seq_name, $data_length);
	my @elements_to_remove = ();
	my $ret_seq = "";

	# open msa file
    if ($msa_format eq "fasta" or $msa_format  eq "gcg" or $msa_format eq "pir"){
        $msa_fh = Bio::SeqIO->new(-file => $msa_file ,
                                -format => $msa_format);        
        unless (defined $msa_fh) {return ('MSA_parser::read_residue_variety : could not read msa $msa_format');}
        # this instruction may throw an exception if parameters are incorrect, for example - if the bioPerl identified too many unvalid characters. In order to catch the exception, we put it in a "catch" block
        eval{
            while ( $seq = $msa_fh->next_seq() ) { # returns a Bio::Seq sequence object
                $num_of_seqs++;
                $ret_seq = set_position($seq->seq(), $seq->id(), $msa_ref_sequence, \@elements_to_remove, $ref_residue_frequency, $ref_position_totalAA);
            }
        };
        if( $@ ) { #Caught exception
            return  ('MSA_parser::read_residue_variety : exception');
        }
        $msa_fh->close();
    }
    elsif($msa_format eq "clustalw" or $msa_format eq "nexus"){
        my $msa_fh = Bio::AlignIO->new(-file => $msa_file ,  
                                -format => $msa_format);
        unless (defined $msa_fh) {return ('MSA_parser::read_residue_variety : could not read msa $msa_format');}
        eval{
            while ( my $aln = $msa_fh->next_aln() ) { # Returns : a Bio::Align::AlignI
                foreach $seq ($aln->each_seq() ){ # Returns : an array of Seq objects
                    $num_of_seqs++;
                    $ret_seq = set_position($seq->seq(), $seq->id(), $msa_ref_sequence, \@elements_to_remove, $ref_residue_frequency, $ref_position_totalAA);
                }
            }
        };
        if( $@ ) {#Caught exception
            return ('MSA_parser::read_residue_variety : exception');
        }
        $msa_fh->close();
    }
    # remove all positions where the query seq was skipped
	my $index = pop(@elements_to_remove);
	while (defined($index))
	{
        residues_shift_left($index, $ref_residue_frequency);
        residues_shift_left($index, $ref_position_totalAA);
		$index = pop(@elements_to_remove);
	}    

	# If there was an unexpected character in the query sequence in the MSA, report it and exit
	if ($ret_seq ne "")
		{return ($ret_seq);}
	else {return ("OK", $num_of_seqs);}
}
#----------------------------------
sub residues_shift_left{
#----------------------------------
    my ($start_position, $ref_residue_frequency) = @_;
    my $index = $start_position;
    
    while (exists($ref_residue_frequency->{$index+1}))
    {
        $ref_residue_frequency->{$index} = $ref_residue_frequency->{$index+1};
        $index++;
    }
    delete($ref_residue_frequency->{$index});    
}
#----------------------------------
sub set_position{
# reads the given sequence, $seq_data. For each position in the sequence, which is not a gap (-) put in the hash $ref_residue_frequency the other residues from the MSA which aligns to this position
#----------------------------------
    my ($seq_data, $seq_name, $msa_ref_sequence, $ref_elements_to_remove, $ref_residue_frequency, $ref_position_totalAA)  = @_;
    my $data_length = length $seq_data;
    my $ret_seq = ""; my $position_in_MSA = 0;
    for (my $position_in_seq = 0;$position_in_seq < $data_length; $position_in_seq++){
        $position_in_MSA++;
        # get residue at this position
        my $current_aa = substr($seq_data,$position_in_seq,1);
        
        # if this is the query seq and no aa is found we need to save
        # this position in order to erase it in the end
        if ($seq_name eq $msa_ref_sequence and $current_aa !~ /^[ABCDEFGHIKLMNPQRSTVWY]$/i){
                push @{$ref_elements_to_remove},$position_in_MSA;
                $ret_seq .= "in seq $seq_name ignored $current_aa in position ".($position_in_seq)." " if ($current_aa !~ /^[-Xx]$/);
        }				
        elsif ($current_aa !~ /^-$/){
            # save residue value
            if (!exists $ref_residue_frequency->{$position_in_MSA}){
                $ref_residue_frequency->{$position_in_MSA}->{$current_aa} = 1;
            }
            elsif (!exists $ref_residue_frequency->{$position_in_MSA}->{$current_aa}){
                $ref_residue_frequency->{$position_in_MSA}->{$current_aa} = 1;
            }
            else{                
                $ref_residue_frequency->{$position_in_MSA}->{$current_aa}++;
            }
            $ref_position_totalAA->{$position_in_MSA}++;
        }
    }
    return $ret_seq;
}
#----------------------------------
1;