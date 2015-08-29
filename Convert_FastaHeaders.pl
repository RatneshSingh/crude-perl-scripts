#!/usr/bin/perl -w
use strict;

#This script Version 1.4 can be used to pull sequences from a fasta based on key words written in pattern file seperated by new line. 
#any white sapces in the pattern will be removed so dont use pattern having white spaces in between. white spaces at end wont be a 
#problem. It reads whole sequences in memory so not advised for very long sequence files e.g >100 MB on computers having small memory.
#vesion  is modified to remove other non recognizable signs from pattern and sequence header .
# version 6 modified for adding X in the last of pattern and header to avoid similar wrong sequence because of similarity.

my($inputfile,$header,@sequence,$sequence,$outputfile,$patternfile,$pattern,$count,@pattern,$ta,%seq_hash);

print " \n****Welcome!!! This script will read pattern to search from the given file and pick out sequence containing pattern from another input file. can be used AS MANUAL MODE***  \n\n\nusage : perl script_name sequence_file outputfile\n\n";

#opening inputfile containing sequences in fasta format.
if ($ARGV[0]){$inputfile=$ARGV[0];}else{print "Enter input file containing sequences\n"; $inputfile=<STDIN>;}
chomp ($inputfile);
open INPUT,"<$inputfile" or die "Cannot open $inputfile.....\n\n";


#opening outputfile
if ($ARGV[1]){$outputfile=$ARGV[1];}else{print "Enter output file name\n"; $outputfile=<STDIN>;}
chomp($outputfile);
open OUT,">$outputfile" or die "cannot create $outputfile.....\n\n";


#opening file containing TA or pattern to search seperated by new line.
if ($ARGV[2]){$patternfile=$ARGV[2];}else{print "Enter pattern file containing pattern seperated by newline \nor Press ENTER for manual entry of pattern\n"; $patternfile=<STDIN>;}
chomp ($patternfile);


#make an array from the input patterns search. CHANGE THE PATTERN RECOGNITION HERE TO FIND SPECIAL KIND OF PATTERN
if ($patternfile ne '') {
	
	open PATTERN,"<$patternfile" or die "Cannot find pattern file $patternfile \n\n";
	print "\n\nReading pattern from pattern file....Plz wait...\n";
	while(<PATTERN>){
	chomp;
	$pattern=$_ if($_ ne /\s*$/);
	chomp($pattern); 
	$pattern=~s/\s+$//g;
	push (@pattern,$pattern);
	next;                
                }
                }
else {
	print "Enter the name of sequence you want to sarch for\n";
	my $name = <STDIN>;
	push(@pattern,$name);
	
}
#*************************************************************************************************************************
#to check if pattern is properly formed and working or not.

$count=@pattern; #just to check if array is formed or not
#~ print "Elements in pattern array \n@pattern\n";
print "Done........\nNumber of Pattern in pattern file =$count\n\n\n";
#~ print "$pattern[0],$pattern[1],$pattern[2],$pattern[3],$pattern[4]";

#*************************************************************************************************************************

#Read database sequence in to a hash %seq_hash{}, remove ">" from header. remove any white spaces and newline from sequence.
print "reading Sequences from input file.....Plz wait...\n";
$/="\n>";

while(<INPUT>){#
    chomp;
    ($header,@sequence)=split("\n",$_);
    $header=~s/>//;
    $sequence= join("",@sequence);
    $sequence=~ s/\s//g;
    $sequence=~s/\n//g;
    #print "$header\n";
    #feed headers and sequences in hash.
    $seq_hash{$header}=$sequence;

    #print "$seq_hash{$header}\n\n";
}#
my @seq_count=keys (%seq_hash);
my $seq_count=@seq_count;

print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
@seq_count=();

#*************************************************************************************************************************
#read all the pattern or key words into array @pattern. and search for the matching one in abovemade hash.
my $count_found=my $count_notfound=0;
print"Looking for patterns...Plz wait....\n";
foreach $ta(@pattern){##
  chomp($ta);
#  my ($header,$sequence)= searchPattern($ta);
  
  my($header,$sequence)= searchExact($ta);
  
        if ($header){print OUT">$header\n$sequence\n"; 
            #~ print " Found Match for $ta\n";
            $count_found=$count_found+1;
        }
        else{ $count_notfound=$count_notfound+1; 
        print"Couldnot found Match for $ta\n";
        }
}
print "\nNumber of Patterns match found for = $count_found\n";
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
	
	my $seqfile=shift(@_);
	
	
	my ($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "reading Sequences from input file.....Plz wait...\n";
	my%seq_hash=();
	$/="\n>";    # Change record seperator to read Fasta

	while(<FASTA>){
    	chomp;
    	($header,@sequence)=split("\n",$_);
    	
    	$header=~s/>//;						# Remove Leading > from Header
    	$header=~s/\s*$//;					# Remove trailing spaces from header
    	$header=~s/^\s*//;					# Remove Leading spaces from Header
    	
    	my$sequence= join("",@sequence);
    	$sequence=~ s/\s//g;
    	$sequence=~s/\n//g;
    	if($header=~/^\s*$/){next;}
    	$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
		
	}
	
	my @seq_count=keys (%seq_hash);
	my $seq_count=@seq_count;

	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.

	return(%seq_hash);

}
#-------------------------------------End ReadFasta---------------------------------------+




#SUBROUTINE TO PRSETHROUGH HASH LOKING FOR PATTERN AND RETURN THE SEQUENCE IF PATTERN MATCHES
sub searchPattern {
      my @pattern=@_;
      #my ($key);
      my $pattern=$pattern[0];
      chomp ($pattern);
	 			   my $mod_pattern=$pattern;
	  				  $mod_pattern=~ s/\[/-/g;
                      $mod_pattern=~ s/\]/-/g;
                      $mod_pattern=~ s/\+/-/g;
                      $mod_pattern=~ s/\|/-/g;
                      $mod_pattern=~ s/\s//g;
                      $mod_pattern=~ s/\./-/g;
                      $mod_pattern=~ s/aa/--/;
                      $mod_pattern=~ s/\*/-/g;
                      $mod_pattern=~ s/\?/-/g;
                      $mod_pattern=~ s/#/-/g; 
 
#					  $pattern=$pattern.'X';
      foreach my $key (keys %seq_hash){
                      my $mod_key=$key;
                      $mod_key=~ s/\|/-/g;
                      $mod_key=~ s/\[/-/g;
                      $mod_key=~ s/\]/-/g;
                      $mod_key=~ s/\+/-/g;
		      		  $mod_key=~ s/aa/--/;
		      		  $mod_key=~ s/\./-/g;
		      		  $mod_key=~ s/\s//g;
		      		  $mod_key=~ s/\*/-/g;
		      		  $mod_key=~ s/\?/-/g;
					  $mod_key=~ s/#/-/g;
#					  $mod_key=$mod_key.'X';
		      		  
		      		  
#print "header-- $mod_key\n pattern-- $pattern\n";
#                    if ($mod_key =~ /$pattern/){ #use this line for matching pattern
 					if ($mod_key eq $mod_pattern){		#use this line for exact pattern
                      	print "For pattern:\t$pattern\t\tPicked seq-->\t$key\n";
                      	return ($key,$seq_hash{$key});
                      }
      next;
                                      }
            }
            

###########################################################################################
#SUBROUTINE TO PRSETHROUGH HASH LOKING FOR PATTERN AND RETURN THE SEQUENCE IF PATTERN MATCHES
sub searchExact {
      my @pattern=@_;
      #my ($key);
      my $pattern=$pattern[0];
      chomp ($pattern);
	  if(exists $seqhash{$pattern}){return ($pattern,$seqhash{$pattern})};
	  else{}
}
      
#########################################################################################
