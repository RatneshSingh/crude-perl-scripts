#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_s,$opt_l,$opt_n);
getopt('sln');

my$usage='
usage: perl script options
-s sequence 
-l length to chop[500]
-n %tolerance to N [5]
';

$opt_n=5 if !defined $opt_n;



die "Cannot find Seqenuce file\n\n$usage" if !defined $opt_s;


my%seq=ReadFasta($opt_s);

#my$out_0='0GC_'.$opt_s;
my$out_20='20GC_'.$opt_s;
my$out_35='35GC_'.$opt_s;
my$out_50='50GC_'.$opt_s;
#my$out_70='70GC_'.$opt_s;

my$GC50=my$GC35=my$GC20=0;
#open OUT0,">$out_0";
open OUT20,">$out_20";
open OUT35,">$out_35";
open OUT50,">$out_50";
#open OUT70,">$out_70";


#my(@bin_50_more,@bin_35_50,@bin_20_35,@bin_0_20);

LOOP: foreach my$name(keys%seq){
	
	
	my$length=length$seq{$name};
	my$GC=($seq{$name}=~tr/GCgc//);
	my$N=($seq{$name}=~tr/Nn//);
	my$GC_percent;
	if($N){my$N_percent=$N*100/$length; next LOOP if $N_percent>$opt_n;}
	if($GC){$GC_percent=$GC*100/$length;}
	else{next;};
	
#	if($GC_percent>70){print OUT70"\n>$name\n$seq{$name}";$GC70++}
	if($GC_percent>50){print OUT50"\n>$name\n$seq{$name}";$GC50++}
	elsif($GC_percent>35){print OUT35"\n>$name\n$seq{$name}";$GC35++}
	elsif($GC_percent>20){print OUT20"\n>$name\n$seq{$name}";$GC20++}
	#else{print OUT0"\n>$name\n$seq{$name}";$GC0++}

}


#print "\n$GC70 Fragments with Gc content >70%";
print "\n$GC50 Fragments with Gc content 50-70%";
print "\n$GC35 Fragments with Gc content 35-50%";
print "\n$GC20 Fragments with Gc content 20-35%\n";
#print "\n$GC0 Fragments with Gc content 0-20%\n";



###########################################################
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

