#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_s,$opt_o,$opt_b,$opt_l);

getopt('sobl');

# define usage 
my$usage = '

This script is to find flanking bases 
of a CpG island predicted by CpG island finder 
programs
Usage: perl script options
-s	sequence file
-o	output
-b	relative base to start[1]
-l	relative base to end [6]
';

#check for options and initiate files
if (defined $opt_s) {if(-e $opt_s){;}else{die "Cannot find sequence file\n $usage";}}
if (!defined $opt_s) {die "Cannot find sequence file\n $usage";}
if (defined $opt_o){open OUT,">$opt_o";}

# set defaulst
$opt_b=1 if !defined $opt_b;
$opt_l=6 if !defined $opt_l;


#read seq file to a hash.

my%seq= ReadFasta($opt_s);


#Read table and count percent of ATGC at +/- 6 positions
my(%ATGC);
my$count=0;

$/="\n";

foreach my$seqname(keys%seq){
MATCH:	while($seq{$seqname}=~m/CG/gi){
		
		# location of CG match on the string
		my$Cposition=$-[0];
		my$Gposition=$+[0];
		
		#check if sequence has any ambigous character. ignore if it does
		
		
		
		if($Cposition-$opt_l<0){next MATCH;}
		if(($Gposition-1+$opt_l)>length($seq{$seqname})){next MATCH;}
		for(my$i=$opt_b;$i<=$opt_l;$i++){
			next MATCH if (substr($seq{$seqname},$Cposition-$i,1)=~/[^ATGCatgc]/); 
		}
		for(my$i=$opt_b;$i<=$opt_l;$i++){
			next MATCH if (substr($seq{$seqname},$Gposition-1+$i,1)=~/[^ATGCatgc]/); 
		}

		
		for(my$i=$opt_b;$i<=$opt_l;$i++){
			# pick base at position $i
			my$nt=substr($seq{$seqname},$Cposition-$i,1);
			
			# Identify found base and count
			#ATGC at +1 to +6 positions in the flanking region
			if($nt eq q(A)||$nt eq q(a)){$ATGC{'A'}{-$i}++;}
			elsif($nt eq q(G)||$nt eq q(g)){$ATGC{'G'}{-$i}++;}
			elsif($nt eq q(T)||$nt eq q(t)){$ATGC{'T'}{-$i}++;}
			elsif($nt eq q(C)||$nt eq q(c)){$ATGC{'C'}{-$i}++;}
			#anything other than ATGC will skew the results. sostart next.
			# comment following line if you dont care.
			else{next TABLE;}
		}	
		

		for(my$i=$opt_b;$i<=$opt_l;$i++){

			my$nt1=substr($seq{$seqname},$Gposition-1+$i,1);
			
			# Identify found base and count
			#ATGC at +1 to +6 positions in the flanking region
			if($nt1 eq q(A)||$nt1 eq q(a)){$ATGC{'A'}{$i}++;}
			elsif($nt1 eq q(G)||$nt1 eq q(g)){$ATGC{'G'}{$i}++;}
			elsif($nt1 eq q(T)||$nt1 eq q(t)){$ATGC{'T'}{$i}++;}
			elsif($nt1 eq q(C)||$nt1 eq q(c)){$ATGC{'C'}{$i}++;}
		}
		$count++;

	}
}

# Number of CpG islands processed
if($opt_o){print OUT"\nTotal number of CpG processed: $count\n";
}

print "\nTotal number of CpG pairs processed: $count\n";


# Calculate percents of ATGC at all the positions in both sides of flanking regions

my(%R_Npcent);

for (my$i=$opt_b;$i<=$opt_l;$i++){

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
else{print "\t";}


for(my$i=-$opt_l;$i<$opt_b-1;$i++){
		if($opt_o){	print OUT"\t$i";}
		else{	print "\t$i";}

}

if($opt_o){	print OUT"\tC\tG";}
		else{	print "\tC\tG";}


for(my$i=$opt_b;$i<=$opt_l;$i++){
		if($opt_o){	print OUT"\t$i";}
		else{	print "\t$i";}

}

# end printing headers


# start printing % values of A,T,G,and C

foreach my$char(keys%ATGC){

	if($opt_o){print OUT"\n$opt_s\t$char";}
	else{print "\n$opt_s\t$char";}
	
	for(my$i=-$opt_l;$i<$opt_b-1;$i++){
		if($opt_o){	printf OUT"\t%.2f",$R_Npcent{$i}{$char};}
		else{	printf "\t%.2f",$R_Npcent{$i}{$char};}

	}
	
	if($opt_o){	print OUT"\t\t";}
		else{	print "\t\t";}

	
	for(my$i=$opt_b;$i<=$opt_l;$i++){
		if($opt_o){	printf OUT"\t%.2f",$R_Npcent{$i}{$char};}
		else{	printf "\t%.2f",$R_Npcent{$i}{$char};}

	}

}

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