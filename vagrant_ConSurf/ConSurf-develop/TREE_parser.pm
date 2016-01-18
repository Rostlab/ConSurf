#!/usr/bin/perl -w

package TREE_parser;
use strict;

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
    my (@lineArr, $tree, $noRegularFormatChar, $treeFileOneLine);
    my $errorBool = 0;
    my $internal_nodes = 0;
    my $bootstrap = 0;

    my $read_right_bracket = "no";
    
    # in case the uploaded tree file contains more than one line -
    # read the tree and rewrite it
    unless (open TREEFILE, $treeFile) {$ref_to_err->{error} = "could not read $treeFile $!";die "could not read $treeFile $!";}
    while (<TREEFILE>) {
       $tree = $_;
       chomp($tree);
       $treeFileOneLine .= $tree;
       $lineCounter++;
    }
    close TREEFILE;   
    $tree =  $treeFileOneLine;
    # add a semi-colon if missing
    if ($tree !~ m/;\s*$/){
        $tree =~ s/\s*$//;
        $tree.=';';
    }   
    if ( $lineCounter>1) {
        unless (open TREEFILE, ">$treeFile") {$ref_to_err->{error} = "could not write to file $treeFile $!";}
        print TREEFILE  $tree; 
        close TREEFILE;
    }
    # legal tree: same number of left and right brackets, no irregular chars
    @lineArr=split(//,$tree); 
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
    if ($leftBrackets != $rightBrackets) {
        $ref_to_err->{left_right}=1;
    } 
    if ($noRegularFormatChar =~ /.+/)  {
       $noRegularFormatChar =~ s/\,\s$//;
       $ref_to_err->{noRegularFormatChar} = $noRegularFormatChar;
    } 
}
#----------------------------------
sub extract_nodes_from_tree{
    
    my ($tree, $ref_tree_nodes) = @_;
    my @tree_arr = split(/\(/, $tree);
    my @sub_tree = ();
    my @temp_arr;
    my $sub_counter = 0;

    # building the array @sub_tree, so that each cell will hold maximum one sequence name
    for(my $i=0; $i<@tree_arr; $i++){
        if ($tree_arr[$i] ne ""){
            $tree_arr[$i] = "(".$tree_arr[$i];
        }
        if ($tree_arr[$i] =~ m/.*,.+/){
            @temp_arr = split(/,/, $tree_arr[$i]);
            foreach (@temp_arr){
                $sub_tree[$sub_counter] = $_.",";
                $sub_counter++;
            }
        }
        else{
            $sub_tree[$sub_counter] = $tree_arr[$i];
            $sub_counter++;        
        }
    }

   # extract the nodes
    my ($exp, $new_rest_exp);
    my $seq_found = "no";
    for (my $k=1; $k<@sub_tree; $k++){
        #in this part we wish to split the expression to 2 parts; left part : (?seq_name ; right part: all the rest       
        if ($sub_tree[$k] ne ""){
            if ($sub_tree[$k] =~ m/(.+)(:.+)/){
                $exp = $sub_tree[$k];
                while ($exp =~ m/(.+)(:.+)/){
                    $exp = $1;
                }
            }
           # in case the expression is of format:  seq_name:distance,
            elsif($sub_tree[$k] =~ m/(.+)(\);.+)/){
                $exp = $1;
                while ($exp =~ m/(.+)(\))/){
                    $exp = $1;
                }
            }
            #  in case the expression is of format:  seq_name)*,
            elsif($sub_tree[$k] =~ m/(.+)(\)?.+)/){
                $exp = $1;
                while ($exp =~ m/(.+)(\))/){
                    $exp = $1;
                }            
            }
            $exp =~ m/(\(?)(.+)/;
            if (exists $ref_tree_nodes->{$2}){
                return ('user_error', "duplicity: $2");
            }
            $ref_tree_nodes->{$2} = 1;
        }
    }
    return ("OK");
}

#----------------------------------
sub removeBPvalues {
   my $IN_treeFile=shift;
   my $OLD_treeFile=shift;
   my $treeFileOneLine;
   open(TREEFILE,"$IN_treeFile");
   while (<TREEFILE>) {
      my $line = $_;
      chomp($line);
      $treeFileOneLine .= $line;
   }
   close TREEFILE;
   my $changed = "no";
   if ($treeFileOneLine =~ m/\)\d*\.?\d+\:/) {
      $treeFileOneLine =~ s/\)\d*\.?\d+\:/\)\:/g; #replace bootstrap values  which look like this: ((A:0.02,B:0.03)40:0.3);
      $changed = "yes";
   }
   if ($treeFileOneLine =~ m/\d*\.?\d+\[\d*\.?\d+\]/) {
      $treeFileOneLine =~ s/(\d*\.?\d+)\[\d*\.?\d+\]/$1/g;#replace bootstrap values  which look like this:(A:0.4,(B:0.1,C:0.1):0.3[40]);
      $changed = "yes";
   }
   if ($changed eq "yes") {
      rename $IN_treeFile, $OLD_treeFile;
      open (TREE_REMOVED,">$IN_treeFile");
      print TREE_REMOVED $treeFileOneLine."\n";
      close TREE_REMOVED;
   }
   return $changed;
}
#----------------------------------
#----------------------------------
1;

