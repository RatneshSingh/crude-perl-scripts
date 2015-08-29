#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_s,$opt_l,$opt_n,$opt_p,$opt_a,$opt_b,$opt_w,$opt_f,$opt_c,$opt_o,$opt_v);
my(%OEratio_bin,$fh);
$opt_f='indiv';
$opt_b=100;
$opt_w=100;
$opt_l=10000;
$opt_n=100;
$opt_a='no';
$opt_c='all';
$opt_v='no';
getopt('snlpabwfcov');

# Write usage rules and options
my$usage="
****************************************************************
This script calculates the Observed/Expected 
ratios of all the dinucleotide pair in a given sequence

Optionally this can be used to calculate the OE ratio
of randomally selected fragments of specified length
from the given sequence for statistical test purposes

Usage:

=>To calculate OE ratio for individual files 
(multiple sequences in each file will be merged)

	perl\tscript\tfile1\tfile2\tfile.........

  
=> to calculate the OE ratio for multiple sequences in a file.
	*for random option, Multiple sequences in sequence file will be joined into
	*one large molecule an will be used for selecting random fragments.

	perl script [options...........]
	
-s	file containing multiple sequences
-f	indiv.	 Calculate OE ratio for individual sequences in file [default].
	random.  Calculate OE ratio in randomaly selected sequence fragments for statistical test purposes.
	bin.     Calculate OE ratio with bin of sliding window.
	randbin  Calculate OE ratio in randomaly selected sequence fragments with bin and sliding window.
-l	Length of randomly selected sequence [10000]
-n	Number of random sequences [100]
-a	yes|no Remove ambigous bases from sequence [no]
-b	Bin size to calculate O/E ratio [100]
-w	Sliding distnce for bin [100].
-c	all|Any Dinucleotide eg CG. dinucleotide to calculate frequency for [all]
-o	output filename to save results [STDOUT]
	if using bin, leave the -o option as program will create informative outputfile name by itself.
-v	Yes|no. Verbose [No]
*****************************************************************
\n\n"
;

print $usage;


#if(!defined @ARGV){print "\nParameters missing!!!\n"; die;}

#my%count;
my@basepairs;
my%basepair;
my@bases=("A","T","G","C");

#create all possible dinucleotide pairs of A,T,G and C and save in hash
foreach my$base1(@bases){foreach my$base2(@bases){my$basep=$base1.$base2;$basepair{$basep}=0;}}


#transfer above created pairs in an array to make sure they appear in same sequence when called
foreach(keys%basepair){push(@basepairs,$_);}


# change the pairs to calculate OE ratio if provided by the user
if(lc$opt_c ne 'all'){$opt_c=~s/\W//g;@basepairs=() if $opt_c ne "";push(@basepairs,$opt_c) if $opt_c ne "";}



# Check if ambigous bases need to be removed. default is No
my $ambigous=0;
if(lc$opt_a eq 'yes'||lc$opt_a eq 'true'){$ambigous=1}
print "Not removing ambigous bases from sequences before OE ratio calculations.\n\n" if $ambigous==0;
print "Removing ambigous bases from sequences before OE ratio calculations.\n\n" if $ambigous==1;



#open output file to save results
if(defined $opt_o){open $fh,'>',"$opt_o" or die "\n******Error opening output file\n";}
elsif(lc$opt_f eq 'bin'){
	my@inputfile=split(/\//,$opt_s);
	
	my$outputfile="OEratioFor.".$opt_c."basepairs.in.".$inputfile[-1].".bin".$opt_b.".Window".$opt_w.".table";
	open $fh,'>',$outputfile or die "\n******Error opening output file:$outputfile\n";}
else{$fh=*STDOUT;}


#create titles for printing.
my$headers="";
for(my$i=0;$i<=scalar@basepairs-1;$i++){$headers.="\t$basepairs[$i]";}

#print {$fh}"$headers\n";



my $count;
# for individual sequences in a file or multiple files.
if(lc$opt_f eq 'indiv'){
	
	# for individuals files with one sequence in each file provided on command line without any flag
	if(!defined $opt_s && $ARGV[0]){
		print {$fh}"$headers\n";
		foreach my$file(@ARGV){
			if($file=~/^-/){next}
			print {$fh}"$file";
			my$sequence=ReadFasta($file,$ambigous);
			my$OEratio=count_OE($sequence,$file,\@basepairs);
			print {$fh}"$OEratio"
		}	
	}
	
	# When one file with multiple sequence is provided
	elsif($opt_s){
		
		my$sequences= ReadFastaTohash($opt_s,$ambigous);
		print {$fh}"$headers\n";
		foreach my$file(keys %{$sequences}){
			#print {$fh}"$file";
			my$sequence=$$sequences{$file};
			my$OEratio=count_OE($sequence,$file,\@basepairs);
			print {$fh}"$OEratio"
		}	
	}

	

}

# for randomly selected sequence 
elsif(lc$opt_f eq 'random'){
	
	print {$fh}"\t";
	for(my$i=0;$i<=scalar@basepairs-1;$i++){print "\t$basepairs[$i]";print {$fh}"\t$basepairs[$i]";}

	my$sequence1= ReadFasta($opt_s,$ambigous); # join multiple sequence to one
		
	for(my$i=1;$i<=$opt_n;$i++){
		my$length=length($sequence1);
		my$random=int(rand($length-$opt_l-1));
		my$sequence=substr($sequence1,$random,$opt_l);
		count_OE($sequence,$i);
	}


	
}
elsif(lc$opt_f eq 'bin'){
		
	die "Cannot find sequence file" if !defined $opt_s;
	
	#open $fh,">OEratioFor.$opt_c.basepairs.in.$opt_s.bin$opt_b.Window.$opt_w.OEratio.table";
	print {$fh}"\t";
	for(my$i=0;$i<=scalar@basepairs-1;$i++){print {$fh}"\tOE ratio for $basepairs[$i] base pair\n";}
		
	my$sequences= ReadFastaTohash($opt_s,$ambigous);
	foreach my$file(keys %{$sequences}){
		my$sequence=$$sequences{$file};
		count_OE_bin($sequence,$file,$opt_b,$opt_w); # adds values to global hash %OEratio_bin
	}
	
	# print OE ratios of all the base pairs for each sequence in a file
	#-------------------------------------------------------------------------------------
	# collect position information from
	my@filenames;
	my@positions;
	my%position;
	foreach my$filename(keys %OEratio_bin){push(@filenames,$filename)}
	foreach my$filename(@filenames){
		foreach my$base(@basepairs){
			foreach(keys %{$OEratio_bin{$filename}{$base}}){
				$position{$_}=0;	
			}
		}
	}
	
	foreach(keys %position){push(@positions,$_)}
	@positions= sort { $a <=> $b }@positions;
	#-------------------------------------------------------------------------------------
	
	# print information 
	foreach my$basepair(@basepairs){
		print {$fh}"Observed/Expected ratio of basepair: $basepair\nBin size:$opt_b\nwindow size:$opt_w\n"; 
		foreach my$posit(@positions){my$median_posit=$posit+($opt_b/2);print {$fh}"\t$median_posit";}
		print {$fh}"\n";
		foreach my$filename(@filenames){
			print {$fh}"$filename($basepair)";
			foreach my$position(keys %{$OEratio_bin{$filename}{$basepair}}){
				
				
				printf {$fh}"\t%0.4f",$OEratio_bin{$filename}{$basepair}{$position};
				
			}
			print {$fh}"\n";
		}
		print {$fh}"\n\n\n";
	}
	
	
	
}
else {print "\nPlease use 'True' or 'False' options only for flag '-r'\n or ";}

print "\n\n";
#printf OUT"\n\nFrequency of Bases:\nA\t%.3f\nT\t%.3f\nG\t%.3f\nC\t%.3f",$freq{'A'},$freq{'T'},$freq{'G'},$freq{'C'};	




#-------------------------------------Start ReadFasta---------------------------------------+

sub ReadFasta{ # to read fasta format files into hash. joins all the sequences in one file into one sequence and returns hash.
	
	my $seqfile=shift(@_);
	my $amb=shift(@_);

	
	
	my ($header,@sequence);
	my$sequence2=();
	my$sequence=();
	chomp $seqfile;
	open FASTA,"$seqfile";
#	print "reading Sequences from input file.....Plz wait...\n";
	my%seq_hash=();
	$/="\n>";    # Change record seperator to read Fasta

	while(<FASTA>){
		chomp;
		($header,@sequence)=split("\n",$_);
		
		$header=~s/>//;						# Remove Leading > from Header
		$header=~s/\s*$//;					# Remove trailing spaces from header
		$header=~s/^\s*//;					# Remove Leading spaces from Header
		
		$sequence= join("",@sequence);
		$sequence=~ s/\s//g;
		$sequence=~s/\n//g;
		if($amb==1){$sequence=~s/[^ATGC]//gi;}	# remove ambigous bases if asked for
		if($header=~/^\s*$/){next;}
			$sequence2.=$sequence;
	}
	
	$/="\n";    							# Record seperator set back to default newline.
	$sequence2=~s/\s+//g;
	return($sequence2);

}
#-------------------------------------End ReadFasta---------------------------------------+

sub ReadFastaTohash{
	
	 my$seqfile=shift(@_);
	 my $amb=shift(@_);

	 my($header,@sequence);
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
    	if($amb==1){$sequence=~s/[^ATGC]//gi;}
    	if($header=~/^\s*$/){next;}
    	$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
		
	}
	
	my @seq_count=keys (%seq_hash);
	my $seq_count=@seq_count;

	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.

	return(\%seq_hash);
}

sub count_OE{

my $sequence=shift;
my $file=shift;
my$basepair=shift;

my%freq;
my%pair_count;
my$print_line;
my$seq_length=length$sequence;
		my$GCpcent=sprintf "%.2f", GC_content($sequence);
		$freq{'A'}=base_frequency($sequence,'A');
		$freq{'T'}=base_frequency($sequence,'T');
		$freq{'C'}=base_frequency($sequence,'C');
		$freq{'G'}=base_frequency($sequence,'G');
	
		#printf "\n%s\t%.3f",$file,$GCpcent;
		#if($opt_p){printf OUT"\n%s\t%.3f",$file,$GCpcent;}
		#$print_line.="\n$file\t$GCpcent";
		#Calculate pair frequency of pairs sent as array reference $basepair	
		for(my$i=0;$i<=scalar@{$basepair}-1;$i++){
			$pair_count{$$basepair[$i]}{'base_freq'}=sprintf "%.2f", base_frequency($sequence,$$basepair[$i]);
			#printf "\t%.3f",$pair_count{$$basepair[$i]}{'base_freq'};
			#if($opt_p){printf OUT"\t%.3f",$pair_count{$$basepair[$i]}{'base_freq'};}
			#$print_line.="\t$pair_count{$$basepair[$i]}{'base_freq'}";
			
			#$count=0;
		}

		#Calculate Obs/Exp ratios

		#printf "\n%s\tObs/Exp ratio",$file;
		#printf OUT"\n%s\tObs/Exp ratio",$file;
		$print_line.="\n$file";
		for(my$i=0;$i<=scalar@{$basepair}-1;$i++){
#			#print "\tpair:$basepairs[$i]";
			#$count++ while $sequence=~/$basepairs[$i]/gi;
			(my$b1,my$b2)=split(//,$$basepair[$i]);
							
			my$baseOEratio=sprintf "%.2f", $pair_count{$$basepair[$i]}{'base_freq'}/($freq{$b1}*$freq{$b2});
			my$baseOEratio_method2=$pair_count{$$basepair[$i]}{'base_freq'}/((($freq{$b1}+$freq{$b2})/2)^2); # expected frequency is [(G%+C%)/2]to the power2. Based on http://www.jstage.jst.go.jp/article/gi/25/1/53/_pdf
			#printf"\t%.3f",$baseOEratio;
			#printf"\t%.3f",$baseOEratio_method2;
			#printf OUT"\t%.3f",$baseOEratio;
			$print_line.="\t$baseOEratio";
		
			$count=0;
		}
	
$print_line.="\n";
return ($print_line);
}	




sub count_OE_bin{ # returns the OE ratios of sequence in sliding window instead of printing it.

my $sequence=shift;
my $filename=shift;
my $bin=shift;
my $window=shift;
my%pair_count;
my$print_line;
my%freq;
my$seq_length=length$sequence;
		
	my$GCpcent=GC_content($sequence);
	$freq{'A'}=base_frequency($sequence,'A');
	$freq{'T'}=base_frequency($sequence,'T');
	$freq{'C'}=base_frequency($sequence,'C');
	$freq{'G'}=base_frequency($sequence,'G');
	
	#Calculate frequency	
	for(my$i=0;$i<=scalar@basepairs-1;$i++){
		for(my$j=0;$j<$seq_length;$j+=$opt_w){
			my$sub_length=$opt_b;
			if(($j+$opt_b) > $seq_length){$sub_length=$seq_length-$j}
			$OEratio_bin{$filename}{$basepairs[$i]}{$j} = OEratio(substr($sequence,$j,$sub_length),$basepairs[$i]);
			print "\nOEratio is 0 for $basepairs[$i] at position $j plus $opt_w of sequence:$filename\n" if ($OEratio_bin{$filename}{$basepairs[$i]}{$j}<=0 && lc$opt_v eq 'yes');
		}
	}
}






sub OEratio{ # calculate OE ratio of a dinucleotide for given sequence 
	
	my$sequence=shift;
	my$basepair=shift;
	#my$filename=shift;
	$basepair=~s/\W//g;
	#print "\nprint from line 370:$sequence\n";
	return(0) if $sequence=~/^\s*$/; # return zero if sequence is empty
	return(0) if $basepair=~/^\s*$/; # return zero if basepair is empty
	
	my$pair_count=0;
	my%bcount;
	my%bfreq;
	my$bfreq_all=1;
	#print "\nprint from line 372:$basepair\n";
	$pair_count++ while $sequence=~/$basepair/gi;
	#print "\nprint from line 374:$pair_count\n";
	my$pair_freq = $pair_count/length$sequence;
	
	
	print "\nObserved freq of pair $basepair is 0 for sequence:\n$sequence\n" if ($pair_freq<=0 && lc$opt_v eq 'yes');
	return(0) if $pair_freq<=0;
	#(my$b1,my$b2)=split("",$basepair);
		
	my@b=split("",$basepair);
	
	for(my$i=0;$i<scalar@b;$i++){
		next if $b[$i]=~/^\s*$/;
		$bcount{$b[$i]}++ while $sequence=~/$b[$i]/gi;
		$bfreq{$b[$i]}=$bcount{$b[$i]}/length$sequence;
		$bfreq_all*=$bfreq{$b[$i]};
	}
	print "\nExpected freq of pair $basepair is 0 for sequence:\n$sequence\n" if ($bfreq_all<=0 && lc$opt_v eq 'yes');
	return(0) if $bfreq_all<=0;
	
	
	my$OEratio=$pair_freq/$bfreq_all;
	#print "OERatio:$OEratio\n";
	return($OEratio);
}



sub base_frequency{ # o calculate frequency of 
	
	my$sequence=shift;
	my$base=shift;
	my$Anumber=0;
	$Anumber++ while $sequence=~/$base/gi;
	return 0 if $Anumber==0;
	return 0 if length$sequence==0;
	my$frequency=$Anumber/length$sequence;
	return $frequency;
	
}

sub GC_content{
	my$sequence=shift;
	my$GCnumber=$sequence=~tr/GCgc//;
	my$GCcontent=$GCnumber/length$sequence;
	return($GCcontent);
}

sub roundoff{
	
	
}