#!/usr/bin/perl -w

use strict;
use Storable;
use lib "/bioseq/ConSurf";
use lib "/bioseq/bioSequence_scripts_and_constants";
use lib "/bioseq/pupkoSVN/trunk/www/bioSequence_scripts_and_constants";
use GENERAL_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;
use CONSURF_CONSTANTS;
use rate4site_routines;
use rasmol_gradesPE_and_pipe;
use cp_rasmol_gradesPE_and_pipe;
use MSA_parser;

my $stored_data_file = $ARGV[0];
my $stored_form_data = $ARGV[1];
my $vars_ref = retrieve($stored_data_file);
my %VARS = %$vars_ref;
my $form_ref = retrieve($stored_form_data);
my %FORM = %$form_ref;

#
my @gradesPE_Output = ();  # an array to hold all the information that should be printed to gradesPE
# in each array's cell there is a hash for each line from r4s.res.
# POS: position of that aa in the sequence ; SEQ : aa in one letter ;
# GRADE : the given grade from r4s output ; COLOR : grade according to consurf's scale
my %residue_freq = (); # for each position in the MSA, detail the residues 
my %position_totalAA = (); # for each position in the MSA, details the total number of residues

# these arrays will hold for each grade, the residues which corresponds to it.
# there are 2 arrays: in the @isd_residue_color, a grade with insufficient data, *, will classify to grade 10
# in the @no_isd_residue_color, the grade will be given regardless of the * mark
# PLEASE NOTE : the [0] position in those arrays is empty, because each position corresponds a color on a 1-10 scale
my @no_isd_residue_color = ();
my @isd_residue_color = ();
# these variables will be used in the pipe block, for view with FGiJ.
# $seq3d_grades_isd : a string. each position in the string corresponds to that ATOM (from the PDB) ConSurf grade. For Atoms with insufficient data - the grade will be 0
# $seq3d_grades - same, only regardeless of insufficient data
my ($seq3d_grades_isd, $seq3d_grades);

#These variables will hold the length of pdb ATOMS and the lenght of SEQRES/MSA_REFERENCE seq
# The data is filled by cp_rasmol_gradesPE_and_pipe::match_seqres_pdb
my ($length_of_seqres,$length_of_atom);

# programs
my $rate4s = GENERAL_CONSTANTS::RATE4SITE;
my $rate4s_slow = GENERAL_CONSTANTS::RATE4SITE_SLOW;
my $FGiJ_path = "/fgij/"; 

# files
$VARS{run_log_Q} = CONSURF_CONSTANTS::CONSURF_LOGS_DIR."/$VARS{run_number}_Q.log";

if ($FORM{uploaded_TREE} ne "") {$VARS{tree_file}=$VARS{user_tree_file_name};}
else {$VARS{tree_file}="TheTree.txt";}
#print "pdb_FILE:*$FORM{pdb_FILE}*";<STDIN>;
if ($FORM{pdb_FILE} ne "") # User PDB
{
	if ($FORM{pdb_FILE} =~ /([^\/]+)$/){
		$VARS{Used_PDB_Name}=$1;
		}
} 

elsif ($FORM{pdb_ID} ne "") {$VARS{Used_PDB_Name}=$FORM{pdb_ID};} # Given PDB_ID
else {$VARS{Used_PDB_Name}="";} #ConSeq Mode
$VARS{r4s_log} = "r4s.log";
$VARS{r4s_out} = "r4s.res";
$VARS{r4s_slow_log} = "r4s_slow.log";
$VARS{atom_positionFILE} = "atom_pos.txt";
$VARS{gradesPE} = "consurf.grades";
$VARS{rasmolFILE} = "rasmol.scr";
$VARS{rasmol_isdFILE} = "isd_rasmol.scr";
$VARS{pipeFile} ="$VARS{working_dir}/$VARS{Used_PDB_Name}" .  "_consurf" . $VARS{run_number} . "_pipe.pdb";
$VARS{Server_Results_Path}="results/$VARS{run_number}";
$VARS{FGiJ_link} = $FGiJ_path . "fg.htm?mol=/$VARS{Server_Results_Path}/". $VARS{Used_PDB_Name} .  "_consurf" . $VARS{run_number} . "_pipe.pdb";
$VARS{Confidence_link}="http://consurf.tau.ac.il/overview.html#CONFIDENCE";

#chimera files
$VARS{chimerax_script_for_figure} = $VARS{Used_PDB_Name}.'_consurf_'.$VARS{run_number}.'_Figure.chimerax';
$VARS{chimerax_script_for_figure_isd} = $VARS{Used_PDB_Name}.'consurf_'.$VARS{run_number}.'_Figure_isd.chimerax';

$VARS{chimera_color_script} = "/chimera/chimera_consurf.cmd";
$VARS{chimera_instructions} = "chimera_instructions.html";

$VARS{scf_for_chimera} = $VARS{Used_PDB_Name} .  "_consurf_" . $VARS{run_number}. ".scf"; 
$VARS{header_for_chimera} = $VARS{Used_PDB_Name} .  "_consurf_" . $VARS{run_number}. ".hdr";
$VARS{chimerax_file} = $VARS{Used_PDB_Name} .  "_consurf_" . $VARS{run_number}. ".chimerax";

$VARS{isd_scf_for_chimera} = $VARS{Used_PDB_Name} .  "_consurf_" . $VARS{run_number}. "_isd.scf";
$VARS{isd_header_for_chimera} = $VARS{Used_PDB_Name} .  "_consurf_" . $VARS{run_number}."_isd.hdr";
$VARS{isd_chimerax_file} =   $VARS{Used_PDB_Name}.  "_consurf_" . $VARS{run_number}. "_isd.chimerax";

$VARS{insufficient_data_pdb} = "";

# Atoms Section with consurf grades instead TempFactor Field
$VARS{ATOMS_with_ConSurf_Scores} = $VARS{Used_PDB_Name} . "_ATOMS_section_With_ConSurf.pdb";
$VARS{ATOMS_with_ConSurf_Scores_isd} =  $VARS{Used_PDB_Name}. "_ATOMS_section_With_ConSurf_isd.pdb";

&open_log_file;
open OUTPUT, ">>$VARS{working_dir}/$VARS{output_page}" or &exit_on_error('sys_error',"create_output_php : could not open the file $VARS{working_dir}/$VARS{output_page} for writing $!");
$VARS{insufficient_data}="no";

#---------------------------------------------
# mode : no msa - with PDB or without PDB
#---------------------------------------------
if ($VARS{running_mode} eq "_mode_pdb_no_msa" or $VARS{running_mode} eq "_mode_no_pdb_no_msa"){        
    $VARS{protein_MSA} = "query_msa.aln";
    &create_MSA;
    if ($VARS{running_mode} eq "_mode_pdb_no_msa")
    {
    	$VARS{msa_SEQNAME}="Input_pdb_SEQRES_A"; #CONSIDER TO CHANGE    
    }
}
#---------------------------------------------
# mode : include msa
#---------------------------------------------
elsif ($VARS{running_mode} eq "_mode_pdb_msa" or $VARS{running_mode} eq "_mode_msa" or $VARS{running_mode} eq "_mode_pdb_msa_tree" or $VARS{running_mode} eq "_mode_msa_tree"){
    #$VARS{protein_MSA} = $VARS{user_msa_file_name};
    $VARS{protein_MSA} = $VARS{user_msa_fasta};
    $VARS{msa_SEQNAME}=$FORM{msa_SEQNAME};
}    
&run_rate4site;
&assign_colors_according_to_r4s_layers(\@gradesPE_Output);
&read_residue_variety(\%residue_freq, \%position_totalAA); # put value in $VARS{num_of_seqs_in_MSA}

#------------------------------
# TEST: to print the content of the residue frequency hash
#foreach my $key (sort { $a <=> $b } (keys %residue_freq)){
#    print OUTPUT "$key ";
#    foreach my $keyy (keys %{$residue_freq{$key}}){
#        print OUTPUT "$keyy $residue_freq{$keyy}";
#    }
#    print OUTPUT "<br />\n";
#------------------------------
#}
#---------------------------------------------
# mode : include pdb
#---------------------------------------------
# in order to create 3D outputs, we need to compare the ATOM to the sequence from rate4site 
if ($VARS{running_mode} eq "_mode_pdb_no_msa" or $VARS{running_mode} eq "_mode_pdb_msa" or $VARS{running_mode} eq "_mode_pdb_msa_tree"){
    &create_atom_position_file; # this file will be used later to create the output which aligns rate4site sequence with the ATOM records
    my %r4s2pdb = (); # key: poistion in SEQRES/MSA, value: residue name with position in atom (i.e: ALA22:A)
    &match_pdb_to_seq(\%r4s2pdb);            
    &create_gradesPE(\%r4s2pdb); # this routine, $VARS{working_dir}apart from creating the file "consurf.grades" also collects the information in order to create the rasmol scripts and a variable which will be used in the "pipe" file, that holds a string with the grades. in the pipe file this var is called: seq3d_grades_isd and seq3d_grades
    &create_rasmol; #This will create the 2 rasmol scripts (one with isd and one without)
    &create_pipe_file; # This will create the pipe file for FGiJ
    &replace_TmpFactor_Consurf; #Will replace the TempFactor Column with the ConSurf Grades (will create also isd file if relevant)
    &create_chimera; #This Will create the script for chimera coloring
}

#foreach my $key (sort { $a <=> $b } (keys %residue_freq)){
#    print OUTPUT "$key ";
#    foreach my $keyy (keys %{$residue_freq{$key}}){
#        print OUTPUT "$keyy $residue_freq{$keyy}";$VARS{unique_seqs};
#    }
#    print OUTPUT "<br />\n";
#}
#print "<br />\nposition_totalAA<br />\n";
#
#foreach my $key (sort { $a <=> $b }(keys %position_totalAA)){
#    print OUTPUT "$key $position_totalAA{$key}<br />\n";
#}

#---------------------------------------------
# Arrange The HTML Output File
#---------------------------------------------
print OUTPUT "\n<H1><center><a name=finish>ConSurf calculation is finished:</a></center></H1>\n";
open OUTPUT, ">>$VARS{working_dir}/$VARS{output_page}" or &exit_on_error('sys_error',"create_output_php : could not open the file $VARS{working_dir}/$VARS{output_page} for writing $!");
print OUTPUT "<font size=+2><i>Final Results</i></font><br><br>\n"; 
if($FORM{pdb_FILE} eq "" and $FORM{pdb_ID} eq "" and $FORM{chain} eq "") { # No PDB is available - ConSeq Outputs
	}
	
else {	# UNIQ ConSurf Results
	print OUTPUT "<font size=+1><A HREF='".$VARS{FGiJ_link}."' TARGET=_blank><b>View ConSurf Results</b></A> with FirstGlance in Jmol</font><br>\n";
	print OUTPUT "<A HREF= \"javascript: top.start_pe('".$VARS{pipeFile}."');\"<b>View ConSurf Results</b></A> with Protein Explorer (Windows only)<br>\n";
	print OUTPUT "<small>(<A HREF='" .$FGiJ_path . "consurf/fg_vs_pe.htm'>Which is best?</A>)</small><br>\n";
	print OUTPUT "<font size=+1><A HREF='consurf.grades' TARGET=Conservation_window>Amino Acid Conservation Scores, Confidence Intervals and Conservation Colors</A></font>\n";
	print OUTPUT "<br><br><font size=+2><i>RasMol Coloring Scripts</b><br></font></i>\n";
   	print OUTPUT "<font size=+1><A HREF='$VARS{rasmol_isdFILE}' TARGET=spt_window> RasMol Coloring Script Showing Insufficient Data</A></font><br>\n";
   	print OUTPUT "<font size=+1><A HREF='$VARS{rasmolFILE}' TARGET=spt_window> RasMol Coloring Script Hiding Insufficient Data</A></font><br>\n";
	print OUTPUT "<br><br><font size=+2><i>PDB Files</b><br></font></i>\n";
  	print OUTPUT "<font size=+1><A HREF='' TARGET=PDB_window> PDB File with Conservation Scores in the tempFactor field</A></font><br>\n";
   	print OUTPUT "<font size=+1><A HREF='$VARS{Used_PDB_Name}_consurf$VARS{run_number}_pipe.pdb' TARGET=PDB_window> PDB File with ConSurf Results in its Header, for FirstGlance in Jmol or Protein Explorer</A></font><br>\n";
	print OUTPUT "<font color=\"green\">Highly recommended for Chimera users:</font><br>\n<b><font size=+1> <a href = \"";
      # if there was insufficient data, we add also the chimerax file showing insufficient data
	if ($VARS{insufficient_data} eq "yes" and -e "$VARS{working_dir}/$VARS{isd_chimerax_file}" and !-z "$VARS{working_dir}/$VARS{isd_chimerax_file}" and -e "$VARS{working_dir}/$VARS{isd_header_for_chimera}" and !-z "$VARS{working_dir}/$VARS{isd_header_for_chimera}"){
         print OUTPUT "$VARS{isd_chimerax_file}";
     }
     else{
         print OUTPUT "$VARS{chimerax_file}";
     }
     print OUTPUT "\"  type=\"application/x-chimerax\">View ConSurf results</a></b> with Chimera</font> (<a href=\"".GENERAL_CONSTANTS::CHIMERA_DOWNLOAD."\">Download Chimera</a>)<br>\n";
     print OUTPUT "<small>Chimera will open the molecule, tree and alignment, coloured by conservation.";
     
     if ($VARS{insufficient_data} eq "yes" and -e "$VARS{working_dir}/$VARS{isd_chimerax_file}" and !-z "$VARS{working_dir}/$VARS{isd_chimerax_file}" and -e "$VARS{working_dir}/$VARS{isd_header_for_chimera}" and !-z "$VARS{working_dir}/$VARS{isd_header_for_chimera}"){
         print OUTPUT " If you wish to view the molecule and avoid the <a href=/overview.html#CONFIDENCE>insufficient data</a>, use <a href = \"$VARS{chimerax_file}\" type=\"application/x-chimerax\">this link</a> please."
     }
# }
# elsif(-e $WorkingDir.$chimerax_script_for_figure and !-z $WorkingDir.$chimerax_script_for_figure){
#     print OUTPUT "$new_icon <font color=\"green\">Highly recommended for Chimera users:</font> $new_icon <br>\n<b><font size=+1> <a href = \"".$chimerax_script_for_figure."\"  type=\"application/x-chimerax\">View ConSurf results</a></b> with Chimera</font> (<a href=\"".GENERAL_CONSTANTS::CHIMERA_DOWNLOAD."\">Download Chimera</a>)<br>\n";
# }
# 
# print OUTPUT "</small><br>\n";
# print OUTPUT "<h4><u>Output Files:</u></h4>\n";
# print OUTPUT "<b><font color=green>Download all ConSurf outputs in a <a href=\"$zip_output\">click!</a></font></b><br><br>\n" if (-e $WorkingDir.$zip_output and !-z $WorkingDir.$zip_output);
# print OUTPUT "&nbsp;&nbsp;&nbsp;<A HREF= $final_out TARGET=Conservation_window>Amino Acid Conservation Scores, Confidence Intervals and Conservation Colors</A><br><br>\n";






# Create a high resolution figure
#    Produce a figure of your protein, colored by ConSurf's colors:
# 
#     * Follow the instructions to produce a PyMOL figure (For users of PyMOL)
#     * Follow the instructions to produce a Chimera figure (For users of Chimera)
	}

print OUTPUT "<br><br><font size=+2><i>Sequences Data</b><br></font></i>\n";
print OUTPUT "&nbsp;&nbsp;&nbsp;<font size=+1><A HREF='$VARS{BLAST_out_file}' TARGET=Blast_window>PSI-BLAST output</A> (PSI-BLAST hits with E-values and pairwise alignments)<br></font>\n";
print OUTPUT "&nbsp;&nbsp;&nbsp;<font size=+1><A HREF= '$VARS{FINAL_sequences_html}' TARGET=Homologues_window>Unique Sequences Used</A> (displayed in FASTA format, linked to UniProt)<br><br></font>;\n";


#<b>Alignment</b><br>
#&nbsp;&nbsp;&nbsp;<A HREF= 1env.aln TARGET=MSA_window>Multiple Sequence Alignment</A> (in Clustal format)<br><br>
#<small><b>&nbsp;&nbsp;&nbsp;Alignment details</b><br>
#&nbsp;&nbsp;&nbsp;The average number of replacements between any two sequences in the alignment; a distance of 0.01 means that on average, the #expected replacement for every 100 positions is 1.<br>
#&nbsp;&nbsp;&nbsp;<i>Average pairwise distance</i> : 0.22515<br>
#&nbsp;&nbsp;&nbsp;<i>Lower bound</i> : -0<br>

#&nbsp;&nbsp;&nbsp;<i>Upper bound</i> : 0.684321</small><br><br>
#<b>Phylogenetic Tree</b><br>
#&nbsp;&nbsp;&nbsp;<A HREF="$VARS{tree_file}" TARGET=Tree_window>Phylogenetic Tree in Newick format</A><br>
#&nbsp;&nbsp;&nbsp;<A HREF=treeView.html TARGET=TreeView_window>View Phylogenetic Tree </A><br><br>
#<b>RasMol Coloring Scripts</b><br>
#&nbsp;&nbsp;&nbsp;<A HREF= rasmolisd.txt TARGET=spt_window>RasMol Coloring Script Showing Insufficient Data</A><br>
#&nbsp;&nbsp;&nbsp;<A HREF= rasmol.txt TARGET=spt_window>RasMol Coloring Script Hiding Insufficient Data</A><br><br>

#print OUTPUT "<font size=+1>Sequences";

close OUTPUT;

stop_reload("$VARS{working_dir}/$VARS{output_page}");
close OUTPUT;
close LOG;


#---------------------------------------------
sub create_MSA{
#---------------------------------------------
    my $cmd;
    chdir $VARS{working_dir};
    if($FORM{MSAprogram} eq 'CLUSTALW'){
        $cmd = GENERAL_CONSTANTS::CLUSTALW." -infile=$VARS{FINAL_sequences} -outfile=$VARS{protein_MSA}";
    }
    else{
        $cmd = GENERAL_CONSTANTS::MUSCLE." -in $VARS{FINAL_sequences} -out $VARS{protein_MSA} -clwstrict -quiet";
    }
    print LOG "create_MSA : run $cmd\n";
    chdir $VARS{working_dir};
    `$cmd`;
    if (!-e "$VARS{working_dir}/$VARS{protein_MSA}" or -z "$VARS{working_dir}/$VARS{protein_MSA}"){
        exit_on_error('sys_error',"create_MSA : the file $VARS{working_dir}/$VARS{protein_MSA} was not created or of size zero");
    }
    $VARS{msa_format} = "clustalw";
    # remove the dnd file    
    $cmd = "rm *.dnd";
    chdir $VARS{working_dir};
    `$cmd`;
}
#---------------------------------------------
sub open_log_file{
#---------------------------------------------
    if (!-e $VARS{run_log_Q} or -z $VARS{run_log_Q}){
        open LOG, ">".$VARS{run_log_Q} or exit_on_error('sys_error', "Cannot open the log file $VARS{run_log_Q} for writing $!");
        print LOG "--------- ConSurf Log $VARS{run_number}_Q.log -------------\n";
        print LOG "Begin Time: ".(BIOSEQUENCE_FUNCTIONS::printTime)."\n";
    }
    else{
        open LOG, ">>".$VARS{run_log_Q} or exit_on_error('sys_error', "Cannot open the log file $VARS{run_log_Q} for writing $!");
    }
}
#---------------------------------------------
sub exit_on_error{
#---------------------------------------------
    my $which_error = shift;
    my $error_msg = shift;
    my $error_definition = "<font size=+1 color='red'>ERROR! ConSurf session has been terminated:</font><br />\n";
    my $syserror = "<font size=+1 color='red'>A SYSTEM ERROR OCCOURED!</font><br />Plesae try to run ConSurf again in a few minutes.<br />We apologize for the inconvenience.<br />\n";
    
    if ($which_error eq 'user_error'){
        print LOG "\nEXIT on error:\n$error_msg\n";
        print OUTPUT  $error_definition."$error_msg";
        
        # print $error_msg to the screen
    }
    elsif($which_error eq 'sys_error'){
        print LOG "\n$error_msg\n";
        print OUTPUT $syserror;
        #print $error_msg to the log file
    }    
    print OUTPUT "<br>
<?    
    include(\"/var/www/html/ConSurf/php/templates/footer.tpl\");
?>";
    close OUTPUT;
    # finish the output page
    sleep 10;
    open OUTPUT, "$VARS{working_dir}/$VARS{output_page}";
    my @output = <OUTPUT>;
    close OUTPUT;
    # remove the refresh commands from the output page
    open OUTPUT, ">$VARS{working_dir}/$VARS{output_page}";
    foreach my $line (@output){
        print OUTPUT $line unless ($line =~ /REFRESH/ or $line =~ /NO-CACHE/);        
    }
    close OUTPUT;
    
    print LOG "\nExit Time: ".(BIOSEQUENCE_FUNCTIONS::printTime)."\n";
    close LOG;
    chmod 0755, $VARS{working_dir};
    exit;
}
#---------------------------------------------
sub run_rate4site{    
#---------------------------------------------  
    my ($cmd, $algorithm, $tree_file_r4s, $query_name, $msa, $did_r4s_fail);    
    my %MatrixHash = (JTT => '-Mj', mtREV => '-Mr', cpREV => '-Mc', WAG => '-Mw', Dayhoff => '-Md');
    # choose the algorithm
    if ($FORM{algorithm} eq "LikelihoodML"){
        $algorithm = "-im";
    }
    else{
        $algorithm = "-ib";
    }
    $tree_file_r4s = '';
    if ($VARS{running_mode} eq "_mode_pdb_msa_tree" or $VARS{running_mode} eq "_mode_msa_tree"){
        $tree_file_r4s = "-t $VARS{user_tree_file_name}";
    }
    if ($VARS{running_mode} eq "_mode_pdb_no_msa" or $VARS{running_mode} eq "_mode_no_pdb_no_msa"){
        $query_name = $VARS{query_string};
        $msa = $VARS{protein_MSA};
        $VARS{rate4site_msa_format} = "clustalw";
    }
    else{
        $query_name = $FORM{msa_SEQNAME};
        print LOG "run_rate4site : Please note: MSA for rate4site run is '$VARS{user_msa_fasta}' (and not original user file : '$VARS{user_msa_file_name}')\n";
        $msa = $VARS{user_msa_fasta};
        $VARS{rate4site_msa_format} = "fasta";
    }
    my $r4s_comm = "$rate4s $algorithm -a \'$query_name\' -s $msa -zn $MatrixHash{$FORM{matrix}} $tree_file_r4s -bn -l $VARS{r4s_log} -o $VARS{r4s_out}";
    print LOG "run_rate4site : running command: $r4s_comm\n";
    chdir $VARS{working_dir};
    `$r4s_comm`;
    $did_r4s_fail = &check_if_rate4site_failed("$VARS{working_dir}/$VARS{r4s_out}", "$VARS{working_dir}/$VARS{r4s_log}");
    # if the run failed - we rerun using the slow verion
    if ($did_r4s_fail eq "yes"){
        print LOG "run_rate4site : The run of rate4site failed. Sending warning message to output.\nThe same run will be done using the SLOW version of rate4site.\n";
        print_message_to_output("<font color='red'><b>Warning:</b></font> The given MSA is very large, therefore it will take longer for ConSurf calculation to finish. The results will be sent to the e-mail address provided.<br>The calculation continues nevertheless.");
        $FORM{send_user_mail} = "yes";
        $r4s_comm = "$rate4s $algorithm -a \'$query_name\' -s $msa -zn $MatrixHash{$FORM{matrix}} $tree_file_r4s -bn -l $VARS{r4s_slow_log} -o $VARS{r4s_out}";
        print LOG "run_rate4site : running command: $r4s_comm\n";
        chdir $VARS{working_dir};
        `$r4s_comm`;
        $did_r4s_fail = &check_if_rate4site_failed("$VARS{working_dir}/$VARS{r4s_out}", "$VARS{working_dir}/$VARS{r4s_slow_log}");
        if ($did_r4s_fail eq "yes"){            
            my $err = "The run $r4s_comm for $VARS{run_number} failed.\nThe MSA $msa was too large.\n";
            exit_on_error('user_error', "The calculation could not be completed due to memory problem, since the MSA is too large. Please run ConSurf again with fewer sequences.")
        }
    }    
}
#---------------------------------------------    
sub check_if_rate4site_failed{
# There are some tests to see if rate4site failed.
# Since I can't trust only one of them, I do all of them. If onw of them is tested to be true - than a flag will get TRUE value
# 1. the .res file might be empty.
# 2. if the run failed, it might be written to the log file of r4s.
# 3. in a normal the r4s.log file there will lines that describe the grades. if it fail - we won't see them
# In one of these cases we try to run the slower version of rate4site.
# We output this as a message to the user.
#---------------------------------------------
    my ($res_flag, $r4s_log) = @_;
    my $ret = "no";
    my @did_r4s_failed = &rate4site_routines::check_if_rate4site_failed($res_flag, $r4s_log);    
    if ($did_r4s_failed[0] eq "yes"){
        $ret = "yes";
        print LOG "check_if_rate4site_failed : ".$did_r4s_failed[1];
        $VARS{r4s_process_id} = $did_r4s_failed[2];        
        &remove_core();
        # if there was an error in user input which rate4site reported: we output a message to the user and exit
        if (exists $did_r4s_failed[3] and $did_r4s_failed[3] ne ""){
            exit_on_error('user_error', $did_r4s_failed[3])
        }
    }
    return $ret;
}
#---------------------------------------------
sub print_message_to_output{
#---------------------------------------------
    my $msg = shift;
    print OUTPUT "\n<ul><li>$msg</li></ul>\n";
}
#---------------------------------------------
sub remove_core{
#---------------------------------------------
    if (-e "$VARS{working_dir}/core.$VARS{r4s_process_id}"){
        print LOG "remove core file : core.".$VARS{r4s_process_id}."\n";
        unlink "$VARS{working_dir}/core.$VARS{r4s_process_id}";
    }
}
#---------------------------------------------
sub assign_colors_according_to_r4s_layers{
#---------------------------------------------
    my $ref_to_gradesPE = shift;
    print LOG "assign_colors_according_to_r4s_layers : $VARS{working_dir}/$VARS{r4s_out}\n";
    my @ans = rasmol_gradesPE_and_pipe::assign_colors_according_to_r4s_layers("$VARS{working_dir}/$VARS{r4s_out}", $ref_to_gradesPE);
    if ($ans[0] ne "OK") {
        exit_on_error('sys_error',$ans[0]);}
}
#---------------------------------------------
sub read_residue_variety{
#---------------------------------------------
    my ($ref_to_res_freq, $ref_to_positionAA) = @_;
    print LOG "read_residue_variety : Calling: MSA_parser::read_residue_variety($VARS{working_dir}/$VARS{protein_MSA}, $VARS{query_string}, $VARS{rate4site_msa_format}, $ref_to_res_freq, $ref_to_positionAA)\n";
    my @ans = MSA_parser::read_residue_variety("$VARS{working_dir}/$VARS{protein_MSA}", $VARS{query_string}, $VARS{rate4site_msa_format}, $ref_to_res_freq, $ref_to_positionAA);
    if ($ans[0] ne "OK") {
        exit_on_error('sys_error',$ans[0]);}
    if (keys %$ref_to_res_freq<1 or (keys %$ref_to_positionAA <1)){
        exit_on_error('sys_error',"could not extract information from MSA $VARS{protein_MSA} in routine MSA_parser::read_residue_variety");}
    
    $VARS{num_of_seqs_in_MSA} = $ans[1];    
}
#---------------------------------------------
sub create_atom_position_file{
#---------------------------------------------
    my $chain;
    unless ($FORM{chain} =~ /none/i){
        $chain = $FORM{chain};}
    else{
        $chain = " ";}
    my %output;
    print LOG "create_atom_position_file : calling rasmol_gradesPE_and_pipe::create_atom_position_file($VARS{working_dir}/$VARS{pdb_file_name},$VARS{working_dir}/$VARS{atom_positionFILE},$chain,\%output)\n";
    rasmol_gradesPE_and_pipe::create_atom_position_file("$VARS{working_dir}/$VARS{pdb_file_name}","$VARS{working_dir}/$VARS{atom_positionFILE}", $chain, \%output);
    if (exists $output{ERROR}) {exit_on_error('sys_error',$output{ERROR});}	
	if (exists $output{INFO}) {print LOG "create_atom_position_file : $output{INFO}";}
	if (exists $output{WARNING}) {print LOG "create_atom_position_file : $output{WARNING}";}
    if (!-e "$VARS{working_dir}/$VARS{atom_positionFILE}" or -z "$VARS{working_dir}/$VARS{atom_positionFILE}"){
        exit_on_error('sys_error',"create_atom_position_file : The file $VARS{working_dir}/$VARS{atom_positionFILE} does not exist or of size 0");
    }
}
#---------------------------------------------
sub match_pdb_to_seq{
#---------------------------------------------    
    my $ref_r4s2pdb = shift;
    my $chain;
    unless ($FORM{chain} =~ /none/i){
        $chain = $FORM{chain};}
    else{
        $chain = " ";}
    print LOG "match_pdb_to_seq : calling cp_rasmol_gradesPE_and_pipe::match_seqres_pdb($VARS{working_dir}/$VARS{pairwise_aln}, $VARS{working_dir}/$VARS{atom_positionFILE}, $chain, $ref_r4s2pdb)\n";
    my @ans = cp_rasmol_gradesPE_and_pipe::match_seqres_pdb("$VARS{working_dir}/$VARS{pairwise_aln}", "$VARS{working_dir}/$VARS{atom_positionFILE}", $chain, $ref_r4s2pdb);
    unless ($ans[0] eq "OK") {exit_on_error('sys_error', "match_pdb_to_seq : rasmol_gradesPE_and_pipe::".$ans[0]);}
    elsif (keys %$ref_r4s2pdb <1 ){exit_on_error('sys_error',"match_pdb_to_seq : Did not create hash to hold r4s and ATOM sequences");}
    print LOG "match_pdb_to_seq : Total residues in the msa sequence: $ans[1]. Total residues in the ATOM : $ans[2]\n";
    $length_of_seqres = $ans[1];
    $length_of_atom = $ans[2];
    #foreach my $resi (sort (keys %$ref_r4s2pdb)){
    #    print OUTPUT "$resi ".$ref_r4s2pdb->{$resi}."<br />";
    #}
}
#---------------------------------------------
sub create_gradesPE(){
#---------------------------------------------
    my $ref_r4s2pdb = shift;
    my @ans=cp_rasmol_gradesPE_and_pipe::create_gradesPE(\@gradesPE_Output, $ref_r4s2pdb, \%residue_freq, \@no_isd_residue_color, \@isd_residue_color, "$VARS{working_dir}/$VARS{gradesPE}");
    unless ($ans[0] eq "OK"){exit_on_error('sys_error', "create_gradesPE : rasmol_gradesPE_and_pipe::".$ans[0]);}
    elsif(!-e "$VARS{working_dir}/$VARS{gradesPE}" or -z "$VARS{working_dir}/$VARS{gradesPE}") {
        exit_on_error('sys_error', "create_gradesPE : the file $VARS{working_dir}/$VARS{gradesPE} was not found or empty");}
    if ($ans[1] eq "" or $ans[2] eq ""){
        exit_on_error('sys_error', "create_gradesPE : there is no data in the returned values seq3d_grades_isd or seq3d_grades from the routine");
    }
    $seq3d_grades_isd = $ans[1];
    $seq3d_grades = $ans[2];
}

#---------------------------------------------
sub create_rasmol(){
# print 2 rasmol files, one showing insufficient data, one hiding it.
#---------------------------------------------
    my $chain;
    unless ($FORM{chain} =~ /none/i){
        $chain = $FORM{chain};}
    else{
        $chain = " ";}
    my @ans;
    print LOG "Calling cp_rasmol_gradesPE_and_pipe::print_rasmol for files $VARS{working_dir}/$VARS{rasmolFILE} and $VARS{working_dir}/$VARS{rasmol_isdFILE}\n";
    @ans = cp_rasmol_gradesPE_and_pipe::print_rasmol("$VARS{working_dir}/$VARS{rasmolFILE}", "no",\@no_isd_residue_color, $chain, "no"); #Without isd residue Color
    unless ($ans[0] eq "OK") {
        exit_on_error('sys_error', "create_rasmol : cp_rasmol_gradesPE_and_pipe::$ans[0]");
        }
    @ans = cp_rasmol_gradesPE_and_pipe::print_rasmol("$VARS{working_dir}/$VARS{rasmol_isdFILE}", "yes",\@isd_residue_color, $chain, "no");#With isd Residue Color
    unless ($ans[0] eq "OK") {
        exit_on_error('sys_error', "create_rasmol : cp_rasmol_gradesPE_and_pipe::$ans[0]");
        }
    if (!-e "$VARS{working_dir}/$VARS{rasmolFILE}" or -z "$VARS{working_dir}/$VARS{rasmolFILE}" or !-e "$VARS{working_dir}/$VARS{rasmol_isdFILE}"  or -z "$VARS{working_dir}/$VARS{rasmol_isdFILE}") {
    exit_on_error('sys_error', "create_rasmol : Did not create one of rasmol outputs");
    }
}
#---------------------------------------------
 sub create_chimera (){
##---------------------------------------------
# Chimera Output includes several files:
# a. header file *.hdr
# b. scf file *.scf
# c. script to show the colored MSA, Tree, and colored 3D structure. (.chimerax)
# d. the Script for Chimera Image.
# e. The Html with the istructions how to create Chimera High resolution Image
#
# A.+B. creating the *.hdr and *.scf files	
# creating view page for the chimera alingment requires the query name for the input sequence in the MSA. In case MSA was uploaded - this name might not be exact, so the option of viewing the alignment only applies for cases where the user did not supply MSA
    print LOG "Calling: cp_rasmol_gradesPE_and_pipe::color_with_chimera($VARS{working_dir}, $VARS{msa_SEQNAME}, $VARS{protein_MSA}, $VARS{r4s_out}, $VARS{scf_for_chimera}, $VARS{header_for_chimera},$VARS{isd_scf_for_chimera},$VARS{isd_header_for_chimera})\n";
    my @ans=cp_rasmol_gradesPE_and_pipe::color_with_chimera($VARS{working_dir}, $VARS{msa_SEQNAME}, $VARS{protein_MSA}, $VARS{r4s_out}, $VARS{scf_for_chimera}, $VARS{header_for_chimera},$VARS{isd_scf_for_chimera},$VARS{isd_header_for_chimera});
    
    if ($ans[0] ne "OK") {exit_on_error('sys_error', "create_chimera: @ans\n");}
    else {$VARS{insufficient_data}=$ans[1];}
    
# C. creating the script that shows the MSA, Tree and colored 3D structure 
	print LOG "Calling: cp_rasmol_gradesPE_and_pipe::create_chimera_script ($VARS{ATOMS_with_ConSurf_Scores},$VARS{run_url}, $VARS{working_dir}, $VARS{chimerax_file},$VARS{protein_MSA},$VARS{tree_file},$VARS{scf_for_chimera}, $VARS{header_for_chimera})\n";
	cp_rasmol_gradesPE_and_pipe::create_chimera_script ($VARS{ATOMS_with_ConSurf_Scores},$VARS{run_url}, $VARS{working_dir}, $VARS{chimerax_file},$VARS{protein_MSA},$VARS{tree_file},$VARS{scf_for_chimera}, $VARS{header_for_chimera});

	# create also ignoring insufficient data file, only in case the selected algorithm was bayes
	if ($VARS{insufficient_data} eq "yes" and $FORM{algorithm} eq "Bayes")
	{
		print LOG "Calling: cp_rasmol_gradesPE_and_pipe::create_chimera_script ($VARS{ATOMS_with_ConSurf_Scores_isd}, $VARS{run_url}, $VARS{working_dir}, $VARS{isd_chimerax_file},$VARS{protein_MSA},$VARS{tree_file},$VARS{isd_scf_for_chimera}, $VARS{isd_header_for_chimera})\n";
		cp_rasmol_gradesPE_and_pipe::create_chimera_script ($VARS{ATOMS_with_ConSurf_Scores_isd},$VARS{run_url}, $VARS{working_dir}, $VARS{isd_chimerax_file},$VARS{protein_MSA},$VARS{tree_file},$VARS{isd_scf_for_chimera}, $VARS{isd_header_for_chimera});
	}
# D. The Script For Chimera Image
	
	print LOG "Calling: cp_rasmol_gradesPE_and_pipe::create_chimera_image_script($VARS{chimerax_script_for_figure},$VARS{ATOMS_with_ConSurf_Scores},$VARS{run_url})\n";
	cp_rasmol_gradesPE_and_pipe::create_chimera_image_script($VARS{chimerax_script_for_figure},$VARS{ATOMS_with_ConSurf_Scores},$VARS{run_url});
	if ($VARS{insufficient_data} eq "yes")
	{
		print LOG "Calling: cp_rasmol_gradesPE_and_pipe::create_chimera_image_script($VARS{chimerax_script_for_figure_isd},$VARS{ATOMS_with_ConSurf_Scores_isd},$VARS{run_url})\n";
		cp_rasmol_gradesPE_and_pipe::create_chimera_image_script($VARS{chimerax_script_for_figure_isd},$VARS{ATOMS_with_ConSurf_Scores_isd},$VARS{run_url});		
	}
 # E. Create Chimera HTML Page	
 	if ($VARS{insufficient_data} eq "yes")
 	{
 		my @ans=cp_rasmol_gradesPE_and_pipe::create_chimera_page ("$VARS{working_dir}/$VARS{chimera_instructions}",$VARS{Confidence_link}, $VARS{chimera_color_script}, $VARS{chimerax_script_for_figure},$VARS{ATOMS_with_ConSurf_Scores}, $VARS{chimerax_script_for_figure_isd}, $VARS{ATOMS_with_ConSurf_Scores_isd});
		if ($ans[0] ne "OK") {&exit_on_error('sys_error',"cp_rasmol_gradesPE_and_pipe::create_chimera_page FAILED: @ans");}
 	}
 	else
	{
		my @ans=cp_rasmol_gradesPE_and_pipe::create_chimera_page ("$VARS{working_dir}/$VARS{chimera_instructions}",$VARS{Confidence_link}, $VARS{chimera_color_script}, $VARS{chimerax_script_for_figure},$VARS{ATOMS_with_ConSurf_Scores});
		if ($ans[0] ne "OK") {&exit_on_error('sys_error',"cp_rasmol_gradesPE_and_pipe::create_chimera_page FAILED: @ans");}
	}
# 		my $chimera_instructions_file = shift;
# 	my $chimera_consurf_commands = shift;
# 	my $conf_link= shift;
# 	
# 	my $chimerax_script = shift;
# 	my $PDB_atoms_ConSurf = shift;
# 	
# 	my $isd_PDB_atoms_ConSurf = shift;
# 	my $isd_chimerax_script =shift;
# 	
	
		
}


#---------------------------------------------
sub create_pipe_file (){
##---------------------------------------------
    # CREATE PART of PIPE	
    print LOG "Calling cp_rasmol_gradesPE_and_pipe::create_part_of_pipe_new $VARS{pipeFile},$VARS{unique_seqs},$FORM{database},seq3d_grades_isd, seq3d_grades, $length_of_seqres, $length_of_atom, $FORM{ESCORE},$FORM{iterations},$FORM{MAX_NUM_HOMOL},$FORM{MSAprogram},$FORM{algorithm},$FORM{matrix}\n";
    my @ans = cp_rasmol_gradesPE_and_pipe::create_part_of_pipe_new("partOfPipe",$VARS{unique_seqs},$FORM{database},$seq3d_grades_isd, $seq3d_grades, $length_of_seqres, $length_of_atom, \@isd_residue_color, \@no_isd_residue_color, $FORM{ESCORE},$FORM{iterations},$FORM{MAX_NUM_HOMOL},$FORM{MSAprogram},$FORM{algorithm},$FORM{matrix});
    unless ($ans[0] eq "OK") {&exit_on_error('sys_error',"cp_rasmol_gradesPE_and_pipe::create_part_of_pipe_new FAILED: @ans");}	
    elsif(!-e "partOfPipe" or -z $VARS{pipeFile}){&exit_on_error('sys_error',"create_pipe_file: The file partOfPipe was not found or empty");}

    print LOG "going to extract data from the pdb, calling: rasmol_gradesPE_and_pipe::extract_data_from_pdb($VARS{working_dir}"."/"."$VARS{pdb_file_name}";
    my @header_pipe = cp_rasmol_gradesPE_and_pipe::extract_data_from_pdb($VARS{working_dir}."/".$VARS{pdb_file_name});
    if ($header_pipe[0] ne "OK"){print LOG @header_pipe;}
	
    my ($msa_filename,$tree_filename,$run_date,$completion_time,$ref_header_title,$IN_pdb_id_capital,$msa_query_seq_name);	

    #GET THE FILE NAMES
    if ($FORM{uploaded_MSA} ne ""){$msa_filename=$FORM{uploaded_MSA};}
    else {$msa_filename="";}
    if ($FORM{uploaded_TREE} ne ""){$tree_filename=$FORM{uploaded_TREE};}
    else {$tree_filename="";}
    
    #Get msa_query_seq_name
    if ($FORM{msa_SEQNAME} ne ""){$msa_query_seq_name=$FORM{msa_SEQNAME};}
    else {$msa_query_seq_name="";}
    
   #GET THE CURRENT TIME
   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
   $year += 1900;
   $mon += 1;
   $run_date = $year . '-' . $mon . '-' . $mday;
   $completion_time = $hour . ':' . $min . ':' . $sec;
   
   $IN_pdb_id_capital=uc($VARS{Used_PDB_Name});
   
   # FIND IDENTICAL CHAINS
   print LOG "Calling cp_rasmol_gradesPE_and_pipe::find_identical_chains_on_PDB_File($VARS{working_dir}/$VARS{pdb_file_name},$FORM{chain})\n";	
   @ans=cp_rasmol_gradesPE_and_pipe::find_identical_chains_on_PDB_File("$VARS{working_dir}/$VARS{pdb_file_name}",$FORM{chain});	
   unless ($ans[0] eq "OK") {&exit_on_error('sys_error',"cp_rasmol_gradesPE_and_pipe::find_identical_chains_on_PDB_File FAILED: @ans");}
   my $identical_chains=$ans[1];
    print LOG "identical_chains:$identical_chains\n";
   # USE THE CREATED PART of PIPE to CREATE ALL THE PIPE TILL THE PDB ATOMS (DELETE THE PART PIPE)
    print LOG "Calling cp_rasmol_gradesPE_and_pipe::create_consurf_pipe_new $VARS{working_dir},$IN_pdb_id_capital,$FORM{chain},\@header_pipe,$VARS{pipeFile},$identical_chains,partOfPipe,$VARS{working_dir},$VARS{run_number},$msa_filename,$msa_query_seq_name,$tree_filename,$VARS{submission_time},$completion_time,$run_date\n";	
    @ans = cp_rasmol_gradesPE_and_pipe::create_consurf_pipe_new($VARS{working_dir},$IN_pdb_id_capital,$FORM{chain},\@header_pipe,$VARS{pipeFile},$identical_chains,"partOfPipe",$VARS{working_dir},$VARS{run_number},$msa_filename,$msa_query_seq_name,$tree_filename,$VARS{submission_time},$completion_time,$run_date);
    unless ($ans[0] eq "OK") {&exit_on_error('sys_error',"cp_rasmol_gradesPE_and_pipe::create_consurf_pipe_new FAILED: @ans");}	
    # Add the PDB data to the pipe
    print LOG "Calling: cp_rasmol_gradesPE_and_pipe::add_pdb_data_to_pipe($VARS{working_dir}/$VARS{pdb_file_name},$VARS{pipeFile})\n";
    @ans=cp_rasmol_gradesPE_and_pipe::add_pdb_data_to_pipe("$VARS{working_dir}/$VARS{pdb_file_name}",$VARS{pipeFile});
    unless ($ans[0] eq "OK") {&exit_on_error('sys_error',"cp_rasmol_gradesPE_and_pipe::add_pdb_data_to_pipe FAILED: @ans");}
    	
    if (!-e "$VARS{pipeFile}" or -z "$VARS{pipeFile}") {
    exit_on_error('sys_error', "create_pipe_file : Did not create the FGiJ output");
    }
}

#----------------------------------------------
sub stop_reload {
##---------------------------------------------
    my $OutHtmlFile = shift;
    sleep 5;
    open OUTPUT, "<$OutHtmlFile";
    flock OUTPUT, 2;
    my @output = <OUTPUT>;
    flock OUTPUT, 8;
    close OUTPUT;
    open OUTPUT, ">$OutHtmlFile";
    
    foreach my $line (@output){    
        if ($line eq "include (\"/var/www/html/ConSurf/php/templates/output_header.tpl\");\n"){
            print OUTPUT "include (\"/var/www/html/ConSurf/php/templates/output_header_no_refresh.tpl\");\n";
        }
	else
	{
	    print OUTPUT $line;
	}	
    }
    close (OUTPUT);
 }
#---------------------------------------------
sub replace_TmpFactor_Consurf {
# This Will create a File containing th ATOMS records with the ConSurf grades instead of the TempFactor column 
	my $chain;
    	unless ($FORM{chain} =~ /none/i){
        	$chain = $FORM{chain};}
    	else{
        	$chain = " ";}
	print LOG  "Calling: cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurf($chain,$VARS{working_dir}/$VARS{pdb_file_name},$VARS{gradesPE},$VARS{ATOMS_with_ConSurf_Scores},$VARS{ATOMS_with_ConSurf_Scores_isd});\n";
	my @ans=cp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurf($chain,"$VARS{working_dir}/$VARS{pdb_file_name}",$VARS{gradesPE},$VARS{ATOMS_with_ConSurf_Scores},$VARS{ATOMS_with_ConSurf_Scores_isd});
	unless ($ans[0] eq "OK") {&exit_on_error('sys_error',"ccp_rasmol_gradesPE_and_pipe::ReplaceTempFactConSurf FAILED: @ans");}
}
#---------------------------------------------
#---------------------------------------------
