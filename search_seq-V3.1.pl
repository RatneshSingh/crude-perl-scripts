#!/usr/bin/perl -w
use strict;
use Getopt::Std;

## find directory for RS_Subs.pm
#use FindBin qw($Bin); # finds the directory of the script
#use lib "$Bin";		  # sets the found directory as lib to search for .pm
#require "RS_Subs.pm";

#This script Version 3.0 can be used to pull sequences from a fasta based on key words written in pattern file seperated by new line.
#any white sapces in the pattern will be removed so dont use pattern having white spaces in between. white spaces at end wont be a
#problem. It reads whole sequences in memory so not advised for very long sequence files e.g >100 MB on computers having small memory.
# vesion  is modified to remove other non recognizable signs from pattern and sequence header .
# version 6 modified for adding X in the last of pattern and header to avoid similar wrong sequence because of similarity.
# 09/25/2013: modified script to save multiple sequences in case multiple hits found for a pattern.
our($inputfile,$header,@sequence,$sequence,$outputfile,$patternfile,$pattern,$count,@pattern,$ta,$seq_hash);
our($opt_s,$opt_o,$opt_l,$opt_m,$opt_t,$opt_d,$opt_c,$opt_h);

$opt_d="FALSE";
$opt_c=1;

getopts('s:o:l:m:t:d:c:h');


my$usage= "
This script will read sequence names and pick out sequence related sequence from given fasta file.
pattern can be provided as list of names stored in a file or typed manually when asked.

\n\nusage : perl script_name [options]

options:
-s	Sequence file containing sequences in fasta format
-l	List file containing pattern seperated by new line.
-t	Table of output file name and sequences to put in those files.
	The name of out file and sequences should be seperated by semicolon, File name first in row.
	e.g. filename1;seq1;seq2;seq3......
		 filename2;seq11;seq12;seq13......

-o	Output file to store resulting sequences.
-m	Search mode to be used. Use following options:[exact]
	'match' to use pattern matching mode. Whole word match required.
	'pmatch' to use partial pattern matching. get sequences with matching pattern.
	'exact' to use exact match mode.
-d	delimiter to split the sequence header and use -c column for search
-c	use this column number for search purpose.
-h	Print help and exit.
";
if($opt_h){print $usage; exit}
if(!$opt_s){print "\n\nSequence file is missing\n\n",$usage}

#opening inputfile containing sequences in fasta format.
if ($opt_s){$inputfile=$opt_s;}else{print "Enter input file containing sequences\n"; $inputfile=<STDIN>;}
chomp ($inputfile);
open INPUT,"<$inputfile" or die "Cannot open $inputfile.....\n\n";


#opening outputfile
if ($opt_o){$outputfile=$opt_o;} elsif(!$opt_o && $opt_t){print "Results will be saved in individual files listed in table."  } else{print "Enter output file name\n"; $outputfile=<STDIN>;}

if($outputfile){
chomp($outputfile) ;
open OUT,">$outputfile" or die "cannot create $outputfile.....\n\n" ;
}




#opening file containing TA or pattern to search seperated by new line.

if ($opt_l && !$opt_t){$patternfile=$opt_l;}
elsif($opt_l && $opt_t){die "Cannot use -t and -l options togather\n";}
elsif(!$opt_l && $opt_t){$patternfile='NA';}
else{print "Enter pattern file containing pattern seperated by newline \nor Press ENTER for manual entry of pattern\n"; $patternfile=<STDIN>;}
chomp ($patternfile);# if $opt_l;

# set search mode to default 'exact' if not requested
if(!$opt_m){$opt_m='exact';}




#make an array from the input patterns search. CHANGE THE PATTERN RECOGNITION HERE TO FIND SPECIAL KIND OF PATTERN
if ($patternfile ne '' && $patternfile ne 'NA') {

	open PATTERN,"<$patternfile" or die "Cannot find pattern file $patternfile \n\n";
	print "\n\nReading pattern from pattern file....Plz wait...\n";
	while(<PATTERN>){
	chomp;
	if($_ =~ /^\s*$/){next;}

	$pattern=$_ ;
	chomp($pattern);
	$pattern=~s/\s+$//g;
	push (@pattern,$pattern);

   }
}

# for seperating multiple sequences in multiple files
elsif($opt_t){
	print "Reading sequence file.........\n";
	$seq_hash=ReadFasta($inputfile,$opt_d,$opt_c);
	my@seqs_not_found;
	open INFILE,"$opt_t" or die "Cannot open table $opt_t\n";
	my $tcount_found=my $tcount_notfound=0;
	while(<INFILE>){

		next if /^\s*$/;
		(my$outputfile,my@seqlist)=split(/;/,$_);
		my $count_found=my $count_notfound=0;
		$outputfile=~s/\s+//g;
		#my@seqlist=split(/;/,$seqnameslist);
		#open OUT2,">$outputfile.fa";


		foreach my$name(@seqlist){
			$name=~s/^\s*//;
			$name=~s/\s*$//;
			chomp($name);
			next if $name=~/^\s*$/;

			if($opt_m =~/match/i){
				my $ref_key_array= searchPattern($name,$seq_hash);

				if(scalar @{$ref_key_array} > 0){
					open OUT2,">$outputfile.fa";
					foreach my$header(@{$ref_key_array}){
						if($$seq_hash{$header}){
							print OUT2 ">$header\n$$seq_hash{$header}\n";
							#~ print " Found Match for $ta\n";
						}
					}
			            #print " Found Match for $name\n";
		            $count_found++;$tcount_found++;
					close OUT2;
		        }
		        else{
					$count_notfound++;$tcount_notfound++;
					print"Couldnot found Match for $name\n";
					push(@seqs_not_found,$name);
		        }


			}

			else{
				my($header,$sequence)= searchExact($name,$seq_hash);

		        if ($header){
					open OUT2,">$outputfile.fa";
					print OUT2">$header\n$sequence\n";
		            #~ print " Found Match for $ta\n";
		            $count_found++;$tcount_found++;
					close OUT2;
		        }
		        else{ $count_notfound++;$tcount_notfound++;
					 push(@seqs_not_found,$name);
		        	 print"Couldnot found Match for $name\n";
		        }

			}

		}

		print "Saved $count_found sequences in $outputfile.fa\n\n";
		print "Number of sequences not found: $count_notfound\n" if $count_notfound>0;


		#close OUT2;
	}
		close INFILE;
		print "\n**** Total Number of sequences found: $tcount_found\n";
		print "Total Number of sequences not found: $tcount_notfound\n";# if $tcount_notfound>0;
		print "\nThese sequences were not found :\n",join "\t",@seqs_not_found if $tcount_notfound>0;

sleep(1);
exit;
}

else {
	print "Enter the names of sequence seperated by newline. press enter again when done\n";
	my $name;
	do{
		$name = <STDIN>;
		$name=~s/^\s*//;
		$name=~s/\s*$//;
		#last if $name=~/^\s*$/;
		push(@pattern,$name)if $name!~/^\s*$/;
	}while($name!~/^\s*$/)
}
#*************************************************************************************************************************
#to check if pattern is properly formed and working or not.

$count=@pattern; #just to check if array is formed or not
#~ print "Elements in pattern array \n@pattern\n";
print "Done........\nNumber of Pattern in pattern file =$count\n\n\n";

$seq_hash=ReadFasta($inputfile,$opt_d,$opt_c);

#*************************************************************************************************************************
#read all the pattern or key words into array @pattern. and search for the matching one in abovemade hash.
my $count_found=my $count_notfound=my$total_found=0;
print"Looking for patterns...Plz wait....\n";
foreach $ta(@pattern){##
  chomp($ta);

	if($opt_m =~/match/i){
			my $ref_key_array= searchPattern($ta,$seq_hash);

			if(scalar @{$ref_key_array} > 0){
				foreach my$header(@{$ref_key_array}){
					if($$seq_hash{$header}){print OUT ">$header\n$$seq_hash{$header}\n";
						#~ print " Found Match for $ta\n";
						$total_found++;
					}
				}
				$count_found++;
			}
			else{ $count_notfound++;
			print"Couldnot found Match for $ta\n";
			}
	}

	else{
			my($header,$sequence)= searchExact($ta,$seq_hash);

	        if ($header&&$sequence){print OUT">$header\n$sequence\n";
	            #~ print " Found Match for $ta\n";
	            $count_found++;
				$total_found++;

	        }
	        else{ $count_notfound++;
	        print"Couldnot found Match for $ta\n";
	        }

	}


}

print "\nNumber of Patterns match found for = $count_found\n";
print "\n\t**-->$total_found sequences saved for $count_found patterns\n";
print "Number of Patterns match not found for = $count_notfound\n";


#*************************************************************************************************************************
#exit the program and close all files.
close INPUT;
close OUT;
close PATTERN;
exit;






#######################################################
#-------------------------------------Start ReadFasta---------------------------------------+

sub ReadFasta{ # to read fasta format files into hash. returns hash.

	my $seqfile=shift;
	my $header_delim=shift;
	my $header_col=shift;


	if ($header_delim ne "FALSE"){print "***********\n\nUse delimiter to split sequence headers for search:$header_delim";
								  print "***********\n\nUsing only Column\#:$header_col for searching \n\n***********\n";
							}
	else{print "***********\n\nUsing full sequence name for searching \n\n***********\n";}

	my ($header,@sequence);
	my $changedheaders='FALSE';
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "reading Sequences from fasta file.....Plz wait...\n";
	my%seq_hash;
	#$seq_hash{'RS_Concatenated'}="";

	$/="\n>";    # Change record seperator to read Fasta
	my$last_N=1;
	while(<FASTA>){
    	chomp;
    	($header,@sequence)=split("\n",$_);

    	$header=~s/>//;						# Remove Leading > from Header
    	$header=~s/\s*$//;					# Remove trailing spaces from header
    	$header=~s/^\s*//;					# Remove Leading spaces from Header


		if($header_col =~ /\d+/i && $header_delim !~/FALSE/i){

			my@header_pieces=split(/\Q$header_delim/,$header);
			$header=$header_pieces[$header_col-1];

		}

    	my$sequence= join("",@sequence);
    	$sequence=~ s/\s//g;
    	$sequence=~s/\n//g;


    	if($header=~/^\s*$/){next;}

    	# Make sure no duplicate header exists else they will be overwritten and lost during hashing due to same key.
    	if(!exists $seq_hash{$header}){
    		$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
			#$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
		}
		else {

			# find a uniq header name by adding a number at the end. If header still exists, increase the number by one
			while(exists $seq_hash{$header}){$header=$header.$last_N;$last_N++;}

			$seq_hash{$header}=$sequence;
			#$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
		    $changedheaders='TRUE'

		}
	}

	my @seq_count=keys (%seq_hash);
	my $seq_count=@seq_count;

	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.
    print "\nFew sequences had identical names. Such names has been changed by adding numbers in the name\n" if $changedheaders eq 'TRUE';
	return(\%seq_hash);

}

#-------------------------------------End ReadFasta---------------------------------------+




##SUBROUTINE TO PRSETHROUGH HASH LOKING FOR PATTERN AND RETURN THE SEQUENCE IF PATTERN MATCHES
sub searchPattern {
      my $pattern=shift;
	  my$seqhash=shift;
	  my@seq_keys;
	  my$found=0;
      #my ($key);
      #$my $pattern=$pattern[0];
      chomp ($pattern);
	  my $mod_pattern=clean_string($pattern);
      foreach my $key (keys %{$seqhash}){
            my $mod_key=clean_string($key);
            if ($mod_key =~ /$mod_pattern/ && $opt_m =~/pmatch/i){ #use this line for matching pattern
             	print "For pattern:\t$pattern\t\tPicked seq-->\t$key\n" if $found==0;
				print "\t\t\t\t\t--->Mulitple match found\t-->\t$key\n" if $found==1;
				$found=1;
               	push(@seq_keys,$key);
            }
			elsif (($mod_key =~ /^\s*$mod_pattern\s*$/ || $mod_key =~ /^$mod_pattern[\D]+/ || $mod_key =~ /[\D]+$mod_pattern[\D]+/||$mod_key =~ /[\D]+$mod_pattern[\D]*$/) && $opt_m =~/match/i){ #use this line for matching pattern
			#elsif (($mod_key =~ /^$mod_pattern[\b]+/ || $mod_key =~ /[\b]+$mod_pattern[\b]+/||$mod_key =~ /[\b]+$mod_pattern[\b]*$/) && $opt_m =~/match/i){ #use this line for matching pattern
             	print "For pattern:\t$pattern\t\tPicked seq-->\t$key\n" if $found==0;
				print "\t\t\t\t\t--->Mulitple match found\t-->\t$key\n" if $found==1;
				$found=1;
               	push(@seq_keys,$key);
            }

      }
	  return (\@seq_keys);
}


###########################################################################################
#SUBROUTINE TO PRSETHROUGH HASH LOKING FOR PATTERN AND RETURN THE SEQUENCE IF PATTERN MATCHES
sub searchExact {
      my $pattern=shift;
	  my$seqhash=shift;
		$pattern=~s/>//;						# Remove Leading > from Header
		$pattern=~s/\s*$//;					# Remove trailing spaces from header
		$pattern=~s/^\s*//;					# Remove Leading spaces from Header

      chomp ($pattern);
	  if(exists $$seqhash{$pattern}){return ($pattern,$$seqhash{$pattern});}
	  else{print "Cannot find $pattern\n ";return 0}
}

#########################################################################################
#################################

sub clean_string{
	my$mod_pattern=shift;
	$mod_pattern=~s/\s*$//g;					# Remove trailing spaces from header
	$mod_pattern=~s/^\s*//g;					# Remove Leading spaces from Header
#	$mod_pattern=~s/([^A-Za-z0-9\s\-\+\=\_])/\\$1/g; # not working. Escapes special characters
	$mod_pattern=~s/([^A-Za-z0-9\s\-\+\=\_])/_/g;
	return($mod_pattern);
}