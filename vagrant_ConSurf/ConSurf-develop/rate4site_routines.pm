#!/usr/bin/perl -w

package rate4site_routines;

use lib "/bioseq/bioSequence_scripts_and_constants";
use GENERAL_CONSTANTS;
use lib "/bioseq/ConSurf";
use CONSURF_CONSTANTS;

my %ColorScale = (0 => 9,
                1 => 8,
                2 => 7,
                3 => 6,
                4 => 5,
                5 => 4,
                6 => 3,
                7 => 2,
                8 => 1);

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
        $print_to_log = "rate4site_routines::check_if_rate4site_failed : the file $res_flag does not exsits. \n";
        $error_found = "yes";
    }
    elsif (-e $res_flag && -z $res_flag) #1
    {
        $print_to_log = "rate4site_routines::check_if_rate4site_failed : the file $res_flag was found to be of size 0. \n";        
        $error_found = "yes";
    }
    if(-e $res_flag){
        unless (open R4S_RES, $res_flag){
            $print_to_log = "rate4site_routines::check_if_rate4site_failed : can not open file: $res_flag. aborting.\n";
            $error_found = "yes";
        }
        while(<R4S_RES>){
            if(/In the tree file there is the name: (.+) that is not found in the sequence file/){
                $print_to_html.= "The sequence name $1 was found in the tree file, but was not found in your MSA.<br>\nPlease correct your tree file, so it will include the same names as they appear in the MSA and re-run your query.<br>\n";
                $error_found = "yes";
                $print_to_log = "rate4site_routines::check_if_rate4site_failed : sequence name $1 was found in the tree file, was not found in the MSA\n";
                last;
            }
        }
        close R4S_RES;
    }
    if (-e $r4s_log && !(-z $r4s_log)) #2,3
    {
        unless (open R4SLOG, $r4s_log) {
            $print_to_log = "rate4site_routines::check_if_rate4site_failed : can not open file: $r4s_log. aborting.\n";
            $error_found = "yes";
        }
        while (<R4SLOG>)
        {
            if (/^.Process_id= (\d+)/){
                $r4s_process_id = $1;
            }
            if ($_ =~ m/likelihood of pos was zero/){
                $print_to_log = "rate4site_routines::check_if_rate4site_failed : the line: \"likelihood of pos was zero\" was found in $r4s_log.\n";
                $error_found = "yes";
                last;
            }
            if ($_ =~ m/rate of pos\:\s\d\s=/){ #if we see this line, we change the flag{                
                $return_ans = "no";
                last;
            }
            if($_ =~ m/The amino-acid sequences contained the character: (.+)/){
                $print_to_html .= "The illegal character $1 was found in your MSA. Please make sure your MSA contains only the following characters:<br />\nA, B, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, X, Z, -";
                $print_to_log = "rate4site_routines::check_if_rate4site_failed : illegal character $1 was found in the MSA\n";
                $error_found = "yes";
                last;
            }
            if($_ =~ m/Could not find a sequence that matches the sequence name/){
                my $seq_name = <R4SLOG>;
                $print_to_log = "rate4site_routines::check_if_rate4site_failed : the submitted query sequence name $seq_name was not found in MSA";
                $print_to_html .= "The query sequence name $seq_name you have submitted was not found in the uploaded MSA. Please note that the sequence name should be exactly as it is in the MSA file<br>";
                $error_found = "yes";
                last;
            }
            if($_ =~ m/The sequence name: (.+)was found in the tree file but not found in the sequence file/){
                my $seq_name = $1;
                $print_to_log = "rate4site_routines::check_if_rate4site_failed : the sequence name $1 was found in the tree file, but not in the MSA";
                $print_to_html .= " The tree file is inconsistant with the uploaded MSA. The sequence: \'$1\' was found in the tree file, but was not found in the MSA.<br>";
                $error_found = "yes";
                last;
            }
            if ($_ =~ m/Bad format in tree file/){
                $print_to_log = "rate4site_routines::check_if_rate4site_failed : ";
                $print_to_html .= " There is an error in the tree file format. Please check that your tree is in the <a href = \"".GENERAL_CONSTANTS::CONSURF_TREE_FAQ."\">requested format</a> and reupload it to the server.<br>";
                $error_found = "yes";
                last;
            }
            if ($_ =~ m/not all sequences are of the same lengths/){
                $print_to_log = "rate4site_routines::check_if_rate4site_failed : problem with the MSA : not all sequences are of the same lengths";
                $error_found = "yes";
                last;
            }
        }
        close R4SLOG;
    }
    $return_ans = "yes" if ($error_found eq "yes");
    return ($return_ans, $print_to_log, $r4s_process_id, $print_to_html);
}

1;