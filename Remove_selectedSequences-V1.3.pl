#!/usr/bin/perl -w
 use warnings;
 use strict;
 use Getopt::Std;
#read file for the sequences and load it in hash



our($opt_s,$opt_l,$opt_m,$opt_d,$opt_c);

$opt_m='exact';
$opt_c=1;
getopt('slmdc');

my$usage="perl script -s sequence_file    -l exclusion_list    -m match_type

-m	match|exact
optional:
-d	delimiter to split head and use column -c as sequence name.[space]
-c	Use this column as sequence name after splitting on delimiter -d.[1]

";

open(FASTAFILE,$opt_s) or die"Cannot find the file \n$usage";
open(LISTEXCLUDE,$opt_l) or die"Cannot find the file \n$usage";

print "\n\nArguments provided:\n$opt_s\tAs Fasta file\n$opt_l\tAs Exclusion List\n";

open(OUT2,">$opt_s.removed.fasta");
open(OUT3,">$opt_s.include.fasta");

my (@exclusion_list,@header_list,%exclusion_list);
#define record seperator ($/) as "\>"

while(<LISTEXCLUDE>){
	chomp $_;
	next if $_=~/^\s*$/;
	$_=~s/([^A-Za-z0-9\s\-\+\=\_])/_/g;

	push(@exclusion_list,$_) if ($_ ne '');
	$exclusion_list{$_}=$_ if ($_ ne '');

}

my$list_length=@exclusion_list;
print "\n\nNumber of sequences to Remove is : $list_length\n";

$/="\n>";

#seperate header and sequence in a hash called $header and @seqarray
my$included=0;
my$excluded=0;
while (<FASTAFILE>){
	my%seq;

	chomp;
	#split record into header and sequence and feed to $header and @seq
	(my$header,my@seq)= split (/\n/,$_);
	$header=~s/>//;
	#$header=~s/\\/_/g;
	my@headsplit=split(/\s+/,$header) if $opt_d && $opt_d =~/space/i;
	@headsplit=split(/$opt_d/,$header) if $opt_d && $opt_d !~/space/i;
	$header=$headsplit[$opt_c-1] if $opt_d;
	$header=~s/\s*$//g;

	#join elements of @seq array to remove any spaces.
	my$seq=join ("",@seq);
	chomp ($header,$seq);
	my$mod_header=$header;
	$mod_header=~s/([^A-Za-z0-9\s\-\+\=\_])/_/g;
	$seq{$header}{'modheader'}=$mod_header;
	$seq{$header}{'seq'}=$seq;

	if (lc$opt_m eq 'exact' && exists $exclusion_list{$mod_header}){
		print OUT2">$header\n$seq\n";$excluded++;
		delete $exclusion_list{$mod_header};
	}
	elsif(lc$opt_m eq 'match' && $seq{$header}{'modheader'} =~ /$exclusion_list{$mod_header}/){
		print OUT2">$header\n$seq\n";$excluded++;
		delete $exclusion_list{$mod_header};

	}else{
		print OUT3 ">$header\n$seq\n";
		$included++;
	}

}

print "\n\nTotal number of sequences removed are: $excluded\nAre saved in file named:$opt_s.removed.fasta\n\n";
print "Total number of sequences remaining in the file: $included\nAre saved in file named:$opt_s.include.fasta\n\n";

print "following sequences could not be removed:\n" if scalar keys %exclusion_list > 0;



foreach (keys %exclusion_list){

	print $_."\n";
}



close FASTAFILE;
close OUT2;
close OUT3;
exit;
