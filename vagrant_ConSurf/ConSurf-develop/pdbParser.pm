#!/usr/bin/perl -w
##############################################################################
# pdbParser.pm
#
#   This module is used for simple parsing of pdb files.
#   After parsing a pdb file via read(filename) function, several items can
#   be retrieved via get("item name")
#
#   Available item names:
#       "HEADER.classification" - as in pdb specification
#       "HEADER.depDate" - as in pdb specification
#       "HEADER.idCode" - as in pdb specification
#       "NMR" - if defined with value 1 means that more than one NMR model was found.
#               in this case the pdb file is shorten to include only the first model
#       "SEQRESx.raw" - amino acid list of chain x in aa one letter format (SEQRES)
#       "SEQRES_chains" - array of chains, according to SEQRES
#       "MODIFIED_COUNTx.raw" - number of modified residues to chain x (SEQRES)
#       "MODIFIED_LISTx.raw" - list of modified residues for chain x (SEQRES)
#       "ATOM_chains" - array of chains, according to ATOM
#       "ATOMx.raw" - amino acid list of chain x in aa one letter format (ATOM)
#       "ATOM.xMODIFIED_X" - if defined with value 1 means that in chain x
#               there was a residue which was translated to x (ATOM) - not really
#               possible, becaue it should be under HETATM and not under ATOM
#
#
#
#   Note:
#       The read() function does not clear previous data. It is best to
#       create a new pdbParser for each file parsed.
#
#   Changelog:
#       2007-Nov-10, Ofir Goldenberg - Created
#       2008-Nov-1,  Elana Erez - added ATOM fields
##############################################################################

package pdbParser;

use strict;

# global vars
my %aa;
my %modified_residues;
my $nucleic_acid_flag = "!";

##############################################################################
# new() 
#   initialize object and amino acid hash conversion table
##############################################################################
sub new {
    my ($class_name) = @_;

    my ($self) = {};
    bless ($self, $class_name);
    $self->{'_created'} = 1;
 
    # the conversion table
    # !!! keys in lower case !!!
    $aa{"ala"} = "A"; # Alanine
    $aa{"arg"} = "R"; # Arginine
    $aa{"asn"} = "N"; # Asparagine
    $aa{"asp"} = "D"; # Aspartic acid
    $aa{"cys"} = "C"; # Cysteine
    $aa{"gln"} = "Q"; # Glutamine
    $aa{"glu"} = "E"; # Glutamic acid
    $aa{"gly"} = "G"; # Glycine
    $aa{"his"} = "H"; # Histidine
    $aa{"ile"} = "I"; # Isoleucine
    $aa{"leu"} = "L"; # Leucine
    $aa{"lys"} = "K"; # Lysine
    $aa{"met"} = "M"; # Methionine
    $aa{"phe"} = "F"; # Phenylalanine
    $aa{"pro"} = "P"; # Proline
    $aa{"ser"} = "S"; # Serine
    $aa{"thr"} = "T"; # Threonine
    $aa{"trp"} = "W"; # Tryptophan
    $aa{"tyr"} = "Y"; # Tyrosine
    $aa{"val"} = "V"; # Valine

# known items in SEQRES that will cause the file to be skipped
    $aa{"a"} = $nucleic_acid_flag; 
    $aa{"t"} = $nucleic_acid_flag;
    $aa{"c"} = $nucleic_acid_flag;
    $aa{"g"} = $nucleic_acid_flag;
    $aa{"u"} = $nucleic_acid_flag; # old style rna
    $aa{"da"} = $nucleic_acid_flag; 
    $aa{"dt"} = $nucleic_acid_flag;
    $aa{"dc"} = $nucleic_acid_flag;
    $aa{"dg"} = $nucleic_acid_flag;
    $aa{"du"} = $nucleic_acid_flag;
    $aa{"di"} = $nucleic_acid_flag;
    $aa{"5cm"} = $nucleic_acid_flag;

    $modified_residues{"mse"} =	"met";
    $modified_residues{"asn"} =	"asn";
    $modified_residues{"mly"} =	"lys";
    $modified_residues{"hyp"} =	"pro";
    $modified_residues{"ser"} =	"ser";
    $modified_residues{"thr"} =	"thr";
    $modified_residues{"cme"} =	"cys";
    $modified_residues{"cgu"} =	"glu";
    $modified_residues{"sep"} =	"ser";
    $modified_residues{"kcx"} =	"lys";
    $modified_residues{"mle"} =	"leu";
    $modified_residues{"tpo"} =	"thr";
    $modified_residues{"cso"} =	"cys";
    $modified_residues{"ptr"} =	"tyr";
    $modified_residues{"dle"} =	"leu";
    $modified_residues{"llp"} =	"lys";
    $modified_residues{"dva"} =	"val";
    $modified_residues{"tys"} =	"tyr";
    $modified_residues{"aib"} =	"ala";
    $modified_residues{"ocs"} =	"cys";
    $modified_residues{"nle"} =	"leu";
    $modified_residues{"mva"} =	"val";

    return $self;
}

##############################################################################
# read(filename) 
#   Reads file and extracts needed data
#
# Params:
#   filename - path to plain text pdb file (unzipped)
##############################################################################
sub read {
    my ($self, $file) = @_;
    my $line;
    my $result = "1";
    my $last_record_found = "";
    my @chains = (); # stores the SEQRES fields chains
    my @ATOM_chains = (); # stores the ATOM fields chains    
    my $pdb_line = 0;
    
    # open file
    open (PDBFILE, $file) or return 0;    

    # save filename we just read
    $self->{'_filename'} = $file;
    LINE: while ($line = <PDBFILE>)
    {
        $pdb_line++;
		chomp $line;
        # HEADER record
        if ($line =~ /^HEADER.*/)
        {
            # save data
            $self->{"HEADER.classification"} = substr($line,10,40);
            $self->{"HEADER.depDate"} = substr($line,50,8);
            $self->{"HEADER.idCode"} = substr($line,62,4);
        }
        # SEQRES record
        elsif ($line =~ /^SEQRES.*/)
        {
            
            # get chain id and prev data
            my $chainID = substr($line,11,1);
            my $seq = $self->{"SEQRES$chainID.raw"};
            #print "seq is: $seq\n";
			my $modified_residue_counter = $self->{"MODIFIED_COUNT$chainID.raw"};
			my $modified_residue_list = $self->{"MODIFIED_LIST$chainID.raw"};
			
			# skip to next line if this is a nucleic acid chain
            if (defined($seq) and substr($seq,0,1) eq $nucleic_acid_flag)
            {
            	next;
            }
            
            # check if this chain was already processed
            my $chainfound = 0;
            foreach my $id (@chains)
            {
                if ($id eq $chainID) {$chainfound = 1;}
            }
            if (!$chainfound)
            {
                push(@chains,$chainID);
				$modified_residue_counter = 0;
            }
            
            # convert to one letter format
            foreach my $aminoName (split(" ",substr($line,19,51)))
            {
                # try to convert using regular amino acid hash
                my $shortName = $aa{lc($aminoName)};
                
                # check if residue is identified
                if (!defined($shortName))
                {
					# count this modified residue
					$modified_residue_counter++;
					
					# try to check for known modifed residues
					my $std_residue_name = $modified_residues{lc($aminoName)};
					
					# check if residue is identified
					if (!defined($std_residue_name))
					{
						# set residue name to X
						$shortName = "X";
						
						# add message to front of modified residue list
						my $modified_changed_to_X_msg = "Modified residue(s) in this chain were converted to the one letter representation 'X'\n";
						if (!defined($modified_residue_list))
						{
							$modified_residue_list = $modified_changed_to_X_msg;
						}
						elsif ($modified_residue_list !~ /Modified residue/)
						{
							$modified_residue_list = $modified_changed_to_X_msg.$modified_residue_list;
						}						
					}
					else
					{
						$shortName = $aa{$std_residue_name};
						my $modified_residue_name = uc($aminoName);
						
						# add to modified residue list
						if (!defined($modified_residue_list))
						{
							$modified_residue_list = $modified_residue_name."->".uc($std_residue_name)."\n";
						}
						elsif ($modified_residue_list !~ /$modified_residue_name/)
						{
							$modified_residue_list .= $modified_residue_name."->".uc($std_residue_name)."\n";
						}
					}
				}
		        
				# check if we found nucleic acid
				elsif ($shortName eq $nucleic_acid_flag)
		        {
		           	$seq = $nucleic_acid_flag."101 this chain contains nucleic acid";
		            last;
		        }
                
                # add one letter code to sequence
                $seq .= $shortName;
            }
            # save data
            $self->{"SEQRES$chainID.raw"} = $seq;
			$self->{"MODIFIED_COUNT$chainID.raw"} = $modified_residue_counter;
			$self->{"MODIFIED_LIST$chainID.raw"} = $modified_residue_list;
        }
        elsif($line =~ /^MODEL\s+2/){
            # if it NMR model: we copy only the first model to the pdb file.
            $self->{"NMR"} = 1;
            close PDBFILE;
            
            open PDBFILE, "<$file";
            my @orig_pdb = <PDBFILE>;
            close PDBFILE;
            
            open (PDBFILE, ">".$file) or return 0;
            COPIED_LINE: foreach my $orig_line (@orig_pdb){
                print PDBFILE $orig_line;
                last COPIED_LINE if ($orig_line =~ /^ENDMDL/);
            }
            close PDBFILE;
            
            last LINE;
        }
        elsif ($line =~ /^ATOM.*/)
        {            
            # extract atom data
            my $Res3L = substr ($line, 17, 3);
            $Res3L = lc $Res3L;
            my $ATOM_chain = substr ($line, 21, 1);  # $_chain            
            my $residue_number = substr ($line, 22, 4); # $number
            $residue_number =~  s/\s+//g;
            # seems like these fields are important only for protein explorer
            #my $atomNo .= substr ($line, 7, 5);
            #my $iCode = substr ($line, 26, 2);
            #if ($iCode !~ /[A-Z]/i){
            #    $chainInsCode = '';
            #}
            # get chain previous data
            my $last_residue_num = $self->{"ATOM_RES_NUM$ATOM_chain.raw"};
            if (!defined $last_residue_num){$last_residue_num="NA";}
            my $ATOM_seq = $self->{"ATOM$ATOM_chain.raw"};
            # check if this chain was already processed
            my $chainfound = 0;
            foreach my $id (@ATOM_chains)
            {
                if ($id eq $ATOM_chain) {$chainfound = 1;}
            }
            if (!$chainfound)
            {
                push(@ATOM_chains,$ATOM_chain);
            }
            # update ATOM sequence only if we read the next residue (not the next atom)
            if (($last_residue_num eq "NA") or $last_residue_num < $residue_number){
				if ($last_residue_num eq "NA") {$last_residue_num=0;}
                while ($residue_number!=$last_residue_num+1){ # For Disorder regions
                     $ATOM_seq=$ATOM_seq."X";
					 $last_residue_num++;
                }
                # convert to one letter format
                my $resShortName = $aa{$Res3L};
                if (!defined($resShortName))
                {
                    $resShortName = "X";
                    $self->{"ATOM.$ATOM_chain"."MODIFIED_X"} = 1;
                }
                $ATOM_seq .= $resShortName;                
                $self->{"ATOM_RES_NUM$ATOM_chain.raw"} = $residue_number;
                $self->{"ATOM$ATOM_chain.raw"} = $ATOM_seq;
            }
           
        }              
    }
    
    # save chain list
    #$self->{"chains"} = [@chains];  # OLD
    $self->{"SEQRES_chains"} = [@chains];
    $self->{"ATOM_chains"} = [@ATOM_chains];
    
    # finish up
    close PDBFILE;
    return $result;
    
}
# just return the $nucleic_acid_flag
sub nucleic_acid_flag
{
    return $nucleic_acid_flag;
}

# just return the $nucleic_acid_flag
sub skip_file_code
{
    return $nucleic_acid_flag;
}

# get a parsed value
sub get {
    my ($self, $key) = @_;

    return $self->{$key};
}

1;
