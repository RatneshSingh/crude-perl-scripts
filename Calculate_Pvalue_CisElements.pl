#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
use List::Util qw(sum);
use Statistics::Descriptive;

getopt('mnsdo');

my$usage="\n\n---------------------------------------------------------------------------------------
perl script {options.....}
-s 'real sequences SignalScan results'
-d 'All the promoters signal scan results(formatted by Process Webscan script)' 
-m 'number of datasets to test' [1000]
-n 'number of sequences in each dataset'
-o 'output file to store results'
---------------------------------------------------------------------------------------\n\n";
print "$usage\n";
our($opt_m,$opt_n,$opt_s,$opt_d,$opt_o,%count,%header,%hit);

open OUT,">$opt_o" or die "Cannot open file $opt_o\n$usage";
printf OUT"%-15s %-10s %-10s %-10s %-10s %-10s\n","Cis_element","True_Hits","Expec_Hit","Var","Std","Z-score";

#Assign default value to $opt_m if not provided by the user
$opt_m=1000 if !defined $opt_m;
	# Read sequence signal scan to get values of motifs and their frequencies.
	if(defined $opt_s){open RSCAN,"$opt_s";}
	else{die "\nCannot find file containing Signal Scan results to test\n$usage";}
	$/="\n________________________";
	my @section=();
	while(<RSCAN>){push(@section,$_);}
	close RSCAN; # close file to reduce stress on computer.
	#push(@filenames,$opt_s);
	print "Reading file and creating list of motifs:$opt_s....";
	#------------------------------------------------------------------------

	#extract sequence name from file.
	my@sequenceheader=split(/\n/,$section[0]);
	my$seqname=();
	foreach(@sequenceheader){if(/>/){$_=~s/>//;($seqname)=split(/\s+/,$_);$seqname=~s/\s*//;last;}}
	#push(@seqnames,$seqname);
	$section[1]=~s/_+||-----+[\w\s\W\S]+//g;
	my@motif_line=split(/\n/,$section[1]);
	my%motifs;
	foreach(@motif_line){
		next if $_=~/^\s*$/;
		$_=~s/^\s+//;
		my@motif_info=();
		@motif_info=split(/\s+/,$_);
		$motif_info[0]=~s/\s+//;
		my$motif_name=$motif_info[0];
		next if ($motif_name=~/^\s*$/);
		$count{$motif_name}{$opt_s}++ if ($motif_name=~/$motif_info[0]/);
		#$count{$motif_name}{$seqname}++ if ($motif_name=~/$motif_info[0]/); # use this for real seq name rather than filenames.
		$motifs{$motif_name}=();
	}

print "...Done\nCreating headers of motifs read.....";
	push(my @motif_in_seq,keys%motifs);
	print"\nNumber of Motifs found in sequence file:";
	print scalar@motif_in_seq;
	print"\n";

	
	$/="\n";
	open DATASET,"$opt_d";
	my$first_line=<DATASET>;
	my @header=split(/\t/,$first_line);
	#shift @headers;
	#print "@header\n";
	for(my$i=0;$i<=scalar@header-1;$i++){$header{$header[$i]}=$i;};
	#print "\nheader: @header\n";
	close DATASET;
print "Done\nCreating table from formatted dataset file......";	

#--------------------------------------------------------------------------
	my%file_value;
	my$i=1;
	$/="\n";
	open DATASET,"$opt_d";
	while(<DATASET>){
		my@hit=split(/\t/,$_);
		my$hit=join("\t",@hit);
		$file_value{$hit[0]}=$hit;
		$i++;
		#print "reading $i line in dataset\n"; 
	}

close DATASET;
print "...Done\n";

push(my@dataset_range,keys%file_value);
my$range=scalar@dataset_range;
$range--;
print "Total number of datasets:$range\n";
#--------------------------------------------------------------------------

#	my $range=33518; # manual entry of random number generator limit.
	# pick values for each motifs in true file from n random sequences; 
	foreach (keys %motifs){
		my$motif=$_;
		#print "Collecting frequencies for $motif\n";
		my@hit_values=();
		for(my$i=1;$i<=$opt_m;$i++){
		loop:	my$random= int(rand($range)); # generate a random number
			goto loop if $random ==0; # discard if random number is zero. refind random number.
			my$filename=$random.'.txt';
			my $hit_value=read_hitValues($filename,$header{$motif});
			push(@hit_values,$hit_value);
			#print"Reading $i th time\n";
		}
	
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@hit_values); 

		my $mean = $stat->mean();
		my $var  = $stat->variance();
		my $std  = $stat->standard_deviation();
		$Statistics::Descriptive::Tolerance = 1e-180;
		my $Zscore=($count{$motif}{$opt_s} - $mean)/$std if ($std!=0);
		
		if(defined $Zscore){
		printf OUT"%-15s%10.2f%10.2f%10.2f%10.2f%10.2f\n",$motif,$count{$motif}{$opt_s},$mean,$var,$std,$Zscore;
		printf "%-15s%10.2f%10.2f%10.2f%10.2f%10.2f\n",$motif,$count{$motif}{$opt_s},$mean,$var,$std,$Zscore;	
		}
		else{
		$Zscore=' NA';
		printf OUT"%-15s%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",$motif,$count{$motif}{$opt_s},$mean,$var,$std,$Zscore;
		printf "%-15s%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",$motif,$count{$motif}{$opt_s},$mean,$var,$std,$Zscore;
		}

	}

#--------------------------------------------------------------------------------

sub read_hitValues{
	
	(my$filename,my$position)=@_;
	my %hit=();

	$/="\n";
	if (exists$file_value{$filename}){
		my@hit=split(/\t/,$file_value{$filename});
		for(my$i=0;$i<=scalar@hit;$i++){$hit{$hit[0]}{$i}=$hit[$i];}
		return($hit{$filename}{$position});
		@hit=();
		
	}
	else{print"\nCould not find $filename\n";}
}
#---------------------------------------------------------------------------------
