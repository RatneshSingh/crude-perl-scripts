#!/usr/bin/perl -w
use strict;
use Getopt::Std;



# Define usage here.
my $usage = "\n\nThis script will parse output file from 'LTRFinder' program and present the results in table format
\nusage: perl scriptname -i inputfile -o outputfile\n\n
";

print "\n\n#########################################################################\n\n$usage";
# collect arguments and process to check if all the required items are provided
our($opt_i,$opt_o);
getopt('io');

die "Input file not found\n $usage\n" if !defined $opt_i;
$opt_o = 'parsed_'.$opt_i if !defined $opt_o;

# Open input file and output file
print "\n-------------------------------------------------------------------------\nOpening Input file: $opt_i\n";
open INFILE,"$opt_i" or die "\n Could not find file: $opt_i;\n";

print "\nPreparing output file to save parsed data: $opt_o\n";
open OUT,">$opt_o" or die "\n\n!!!Cannot create output file. Please check your permissions or make sure file is not in use!!!\n";
print OUT"Seq_Name\tLoc_Start\tLoc_end\tLoc_length\tLoc_strand\tScore\tStatus\t5'-LTR\t3'-LTR\t5'-TG\t3'-CA\tTSR\tSharpness\tStrand\tPBS\tPPT";


#Reading input file line by line

my$name;
my$lines;
my$count;
my$notfound;

while(<INFILE>){
	
	chomp($_);	
	$_=~s/\s+$//g;

	# Find line with Sequence name and store name in $name
	if($_ =~ /\&gt\;Sequence:\s*([\w\W]+)$/){ 
		my$seqName=$1; 
		#print "***************\n\nSeqName:\t$seqName\n"; 
		chomp($seqName);
		print OUT"\n$seqName";
		$name=$seqName;
		$count++;
	}
	
	# Looking for the multiple entries from one sequence
	elsif($_ =~/^\[\d+\]/ && $_  !~ /^\[1\]/){
		print OUT"\n$name"; 
		#print "$count2:in side of extension: $_\n $name\n";
	}
#	elsif($_ =~/^\[\d+\]/) {if ($_  !~ /^\[1\]/){print OUT"\n$name"; $count2++; print "$count2:in side of extension: $_\n";}}

	#Parsing for other items
	elsif($_ =~ /Location :\s*([\w\W]+)$/){ 
		my$Location=$1;
		(my$LocStart,my$div,my$LocEnd,my$div2,my$LocLen,my$LocStrand)=split(/\s+/,$Location);
		$LocStart=~s/\D//g;
		$LocEnd=~s/\D//g;
		$LocLen=~s/\D//g;
		$LocStrand =~s/[^\+\-]//g;
		if($LocStrand=~ /\+/){$LocStrand='+';} 
		elsif($LocStrand=~ /\-/){$LocStrand='-';} 
		else{print "Cannot determine strand for $name\n;"}
		
		#print "LocStart:\t$LocStart\nLocend:\t$LocEnd\nLoclength:\t$LocLen\nLocStrand:\t$LocStrand\n";
		chomp($LocStart,$LocEnd,$LocLen,$LocStrand);		
		print OUT"\t$LocStart\t$LocEnd\t$LocLen\t$LocStrand";

	}
	
	

	elsif($_ =~ /Score\s*:\s*([\w\W]+)$/){ my$Score=$1; 
		#print "Score:\t$Score\n"; 
		chomp($Score);
		print OUT"\t$Score";
		$lines++;
	}
	
	elsif($_ =~ /Status\s*:\s*([\w\W]+)$/){ my$Status=$1; 
		#print "Status:\t$Status\n"; 
		chomp($Status);
		$Status=~s/[^\d]//g;	
		print OUT"\t$Status";
	}
	

	elsif($_ =~ /5'-LTR\s*:\s*([\w\W]+)$/){ my$LTR5=$1; 
		#print "5'-LTR:\t$LTR5\n"; 
		chomp($LTR5);		
		print OUT"\t$LTR5";
	}
	

	elsif($_ =~ /3'-LTR\s*:\s*([\w\W]+)$/){ my$LTR3=$1; 
		#print "3'-LTR:\t$LTR3\n"; 
		chomp($LTR3);	
		print OUT"\t$LTR3";
	}
	
	elsif($_ =~ /5'-TG\s*:\s*([\w\W]+)$/){ my$TG5=$1; 
		#print "5'-TG:\t$TG5\n"; 
		chomp($TG5);	
		print OUT"\t$TG5";
	}
	
	elsif($_ =~ /3'-CA\s*:\s*([\w\W]+)$/){ my$CA3=$1; 
		#print "3'-CA:\t$CA3\n"; 
		chomp($CA3);	
		print OUT"\t$CA3";
	}
	
	elsif($_ =~ /TSR\s*:\s*([\w\W]+)$/){ my$TSR=$1; 
		#print "TSR:\t$TSR\n"; 
		chomp($TSR);	
		print OUT"\t$TSR";
	}
	

	elsif($_ =~ /Sharpness\s*:\s*([\w\W]+)$/){ my$Sharpness=$1; 
		#print "Sharpness:\t$Sharpness\n"; 
		chomp($Sharpness);	
		print OUT"\t$Sharpness";
	}
	

	elsif($_ =~ /Strand\s*([+-])/){ my$Strand=$1; 
		#print "Strand:\t$Strand\n"; 
		chomp($Strand);	
		print OUT"\t$Strand";
	}
	
	
	elsif($_ =~ /PBS\s*:\s*([\w\W]+)$/){ my$PBS=$1; 
		#print "PBS:\t$PBS\n";
		chomp($PBS);	 
		print OUT"\t$PBS";
	}
	

	elsif($_ =~ /PPT\s*:\s*([\w\W]+)$/){ my$PPT=$1; 
		#print "PPT:\t$PPT\n";
		chomp($PPT);	 
		print OUT"\t$PPT";
	}
	
	elsif($_=~/No LTR Retrotransposons Found/){$notfound++;}
	
	# If nothing found go next line.
	else{next;}
	
}


print "\nDone parsing......\n\nParsed $lines entries for $count Sequences.\n$notfound sequences did not show any entry for LTRs\n";
print "\n\n#########################################################################\n\n";
close INFILE;
close OUT;
exit;
