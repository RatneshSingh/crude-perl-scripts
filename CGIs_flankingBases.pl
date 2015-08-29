#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_t,$opt_s,$opt_o,$opt_b,$opt_l);

getopt('tsobl');

# define usage 
my$usage = '

This script is to find flanking bases 
of a CpG island predicted by CpG island finder 
programs
Usage: perl script options
-s	sequence file
-t	table output of newcpgreport
-o	output
-b	relative base to start[1]
-l	relative base to end [6]
';

#check for options and initiate files
if (defined $opt_s) {if(-e $opt_s){;}else{die "Cannot find sequence file\n $usage";}}
if (!defined $opt_s) {die "Cannot find sequence file\n $usage";}

if(defined $opt_t){open TABLE,"$opt_t" ;} else{die "Cannot find output from newcpgreport\n $usage";}
if (defined $opt_o){open OUT,">$opt_o";}

# set defaulst
my$start_base=1 if !defined $opt_b;
my$end_base=6 if !defined $opt_l;


#read seq file to a hash.
#my%seq;
my$seq_length;
my%seq= ReadFasta($opt_s);


#Read table and count percent of ATGC at +/- 6 positions
my(%ATGC);
my$count=0;
my$CGI_Total=0;
my$total_count=0;
$/="\n";
TABLE: while (<TABLE>){
			# next if header line or an empty line found.
			if($_=~/Seq/||$_=~/start/||$_=~/end/ ||$_=~/^\s+$/){next;}
	
			# split line and catch values.
			my($seq,$start,$end,$size,$Sum_GC,$percent_GC,$Obs_Exp)=split(/\t+/,$_);
			chomp($seq,$start,$end,$size,$Sum_GC,$percent_GC,$Obs_Exp);
			$CGI_Total+=$size;
			$seq=~s/\s+//g;
			$total_count++;

			#correct string identifier. Firste element is not 1 its 0;
			my$R_start_base=$end-1;
			my$L_start_base=$start-1;
	
	
			# Find actual name to locate sequence in fasta file.
			#	foreach my$seqName(keys %seq){
			#	foreach my$seqName(keys %seq){

		
			#duplicate name to avoid changes in original name
			#		my$seqName2=$seqName;
			#		my$seq2=$seq;
		
			# truncate name as longer names causes problems.
			#$seqName2=substr($seqName,0,5);
			#$seq2=substr($seq,0,5);
		
			# remove non word characters and spaces from names to ease comparison
			#		$seqName2=~s/\W/_/g;
			#		$seq2=~s/\W/_/g;

			#		$seqName2=~s/\s+//g;
			#		$seq2=~s/\s+//g;
		
			#compare and assign the full name to seq variable
			#		if($seqName2=~/$seq2/){$seq=$seqName;last;}
			#		if($seqName=~/$seq/){$seq=$seqName;last;}

			#	}
	
			# check for available bases at all the position to compare.
			#if (($end+$end_base)>length($seq{$seq})){}
	
			# if sequence exists count A,T,G,andC
			if(exists $seq{$seq}){
			#	if($seq{$_}){	
		
		
			#if(($R_start_base+$end_base)>length($seq{$seq})){$end_base=length($seq{$seq})-$R_start_base;}
			if(($R_start_base+$end_base)>length($seq{$seq})){next;}				
			#ATGC at + 6 positions in the flanking region
			for(my$i=$start_base;$i<=$end_base;$i++){
			
				# pick base at position $i
				my$nt=substr($seq{$seq},$R_start_base+$i,1);
			
				# Identify found base and count
				#ATGC at +1 to +6 positions in the flanking region
				if($nt eq q(A)||$nt eq q(a)){$ATGC{'A'}{$i}++;}
				elsif($nt eq q(G)||$nt eq q(g)){$ATGC{'G'}{$i}++;}
				elsif($nt eq q(T)||$nt eq q(t)){$ATGC{'T'}{$i}++;}
				elsif($nt eq q(C)||$nt eq q(c)){$ATGC{'C'}{$i}++;}
				#anything other than ATGC will skew the results. sostart next.
				# comment following line if you dont care.
				else{next TABLE;}
		}	
		
		
		if(!$opt_l){$end_base=6;}
		else{$end_base=$opt_l};

		# use this when need to avoid ends to cause trouble. 
		# this will only count available places
#		if(($L_start_base-$end_base)<=0){$end_base=$L_start_base;}
		
		# Use this if dont want to skew reeuslt due to unavailable bases to compare
		if(($L_start_base-$end_base)<=0){next;}
		#ATGC at - 6 to -1 positions in the flanking region
		for(my$i=$start_base;$i<=$end_base;$i++){

			my$nt1=substr($seq{$seq},$L_start_base-$i,1);
			
			# Identify found base and count
			#ATGC at +1 to +6 positions in the flanking region
			if($nt1 eq q(A)||$nt1 eq q(a)){$ATGC{'A'}{-$i}++;}
			elsif($nt1 eq q(G)||$nt1 eq q(g)){$ATGC{'G'}{-$i}++;}
			elsif($nt1 eq q(T)||$nt1 eq q(t)){$ATGC{'T'}{-$i}++;}
			elsif($nt1 eq q(C)||$nt1 eq q(c)){$ATGC{'C'}{-$i}++;}
			#anything other than ATGC will skew the results. sostart next.
			# comment following line if you dont care.
			else{next TABLE;}

		}
		$count++;

	}

	else{print "Cannot find sequence $seq in file\n";}
	
}

# Number of CpG islands processed
if($opt_o){print OUT"\nTotal number of CpG islands processed: $count\n";
print OUT"\nTotal number of CpG islands used in calculations: $total_count\n";
print OUT"\nTotal length of all the CGIs: $CGI_Total\n";
#print OUT"Total length of all the sequences: $seq_length\n\n\n";
}

print "\nTotal number of CpG islands processed: $count\n";
print "\nTotal number of CpG islands used in calculations: $total_count\n";
print "\nTotal length of all the CGIs: $CGI_Total\n";
#print "Total length of all the sequences: $seq_length\n";



# Calculate percents of ATGC at all the positions in both sides of flanking regions

my(%R_Npcent);

for (my$i=$start_base;$i<=$end_base;$i++){

	$R_Npcent{$i}{'A'}=$ATGC{'A'}{$i}*100/$count;	
	$R_Npcent{$i}{'G'}=$ATGC{'G'}{$i}*100/$count;
	$R_Npcent{$i}{'T'}=$ATGC{'T'}{$i}*100/$count;
	$R_Npcent{$i}{'C'}=$ATGC{'C'}{$i}*100/$count;
	
	$R_Npcent{-$i}{'A'}=$ATGC{'A'}{-$i}*100/$count;	
	$R_Npcent{-$i}{'G'}=$ATGC{'G'}{-$i}*100/$count;
	$R_Npcent{-$i}{'T'}=$ATGC{'T'}{-$i}*100/$count;
	$R_Npcent{-$i}{'C'}=$ATGC{'C'}{-$i}*100/$count;

}

# print the results out. its tricky thats why looks complex
# print headers for results.
if($opt_o){	print OUT"\t";}
else{	print "\t";}


for(my$i=$end_base;$i>=$start_base;$i--){
		if($opt_o){	print OUT"\t-$i";}
		else{	print "\t-$i";}

}



for(my$i=$start_base;$i<=$end_base;$i++){
		if($opt_o){	print OUT"\t$i";}
		else{	print "\t$i";}

}



# Print % values of A,T,G,and C

foreach my$char(keys%ATGC){

	if($opt_o){print OUT"\n$opt_s\t$char";}
	else{print "\n$opt_s\t$char";}
	
	for(my$i=$start_base;$i<=$end_base;$i++){
		if($opt_o){	printf OUT"\t%.2f",$R_Npcent{$i}{$char};}
		else{	printf "\t%.2f",$R_Npcent{$i}{$char};}

	}
	
	for(my$i=$end_base;$i>=$start_base;$i--){
		if($opt_o){	printf OUT"\t%.2f",$R_Npcent{-$i}{$char};}
		else{	printf "\t%.2f",$R_Npcent{-$i}{$char};}

	}

}

close TABLE;
close OUT;









#-------------------------------------Start ReadFasta---------------------------------------+

sub ReadFasta{ # to read fasta format files into hash. returns hash.
	
	my $seqfile=shift(@_);
	
	
	my ($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "reading Sequences from input file.....Plz wait...\n";
	my%seq_hash=();
	my$Total_length=0;
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
		$Total_length+=length$sequence;
	}
	
	my @seq_count=keys (%seq_hash);
	my $seq_count=@seq_count;

	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.
	return(%seq_hash);
	%seq_hash=();

}
#-------------------------------------End ReadFasta---------------------------------------+