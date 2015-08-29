#!/usr/bin/perl -w
#This script will read through hmmsearch output file and collect
#domain infrormation from it into different file.keep this script 
#in the folder containing hmmsearch output and run from there.
#Version:1.0
#Author: Ratnesh Singh


use strict;
my(@lines,@parts);
#open the directory containing files(current directory".").
opendir(DIR, ".");
my$pattern=$ARGV[0] if defined $ARGV[0];
if(!$ARGV[0]){
#	print "what pattern to look for file name\n"
	$pattern='hmmsearch';
}



#read filesnames (readdir(DIR)  having .txt pattern in an array
my @files = grep(/$pattern/,readdir(DIR));
closedir(DIR);
#print "files to be read @files\n";
open(OUT,">hmmsearch_to_table-Result.out");
open(OUT2,">hmmsearch_to_table-Domain.out");

print OUT"Filename\tsequence\tScore\tEvalue\n";



#processing files and printing information one after one.
$/='Parsed for domains';
foreach my $file (@files) {
print "file being read: $file\n";
	open(FILE,"$file") or die "Can't open file $file";
#$/='Parsed for domains';
	while(<FILE>){
#print"default variable contains: $_\n";
			chomp;
#print"\n This is the fragment read in file:$_";
			my$parts=$_;
			push(@parts,$parts);
			}
			
			
			my$first=shift(@parts);
			@parts=();
#			print"array parts assigned:$first";
			@lines=split(/\n/,$first);
#			print"\nlines fed in to an array:@lines";
			for(my$i=1;$i<=21;$i++){
				shift(@lines);
						}
			my$len=@lines;
#			print"\narray after stripping header:@lines";
			for(my$j=0;$j<=$len-1;$j++){
				$lines[$j]=~s/\s+/\t/g;
				my($name,$score,$evalue,$N)=split(/\t/,$lines[$j]);
#				$name=~s/\s*/\t/g;
#				$score=~s/\s*/\t/g;
#				$evalue=~s/\s*/\t/g;
#				$N=~s/\s*/\t/g;
#				my($name,$score,$evalue,$N)=split(/\s*/,$lines[$j]);
				print OUT"$file\t$name\t$score\t$evalue\n";
								
			}

			
			}
			
######################################################################################
# to parse domain information to table

######################################################################################			

$/="\n";

foreach my $file (@files) {

print "file being read: $file\n";
	{
    local( $/, *FILE ) ;
    open( FILE, $file ) or die "sudden flaming death\n";
    my $text = <FILE>;
	#print $text;
	my $Pfam=();
	my @fragment= split(/Alignments(\s)*of(\s)*top-scoring(\s)*domains:/,$text);
	
	if($fragment[0]=~/HMM file:[\s]+(PF[\d]+)/){$Pfam=$1}
	
	$fragment[0]=~ s/hmmsearch [\s\S\w\W\d\D]*?Parsed[\s\S\w\W\d\D]*?E-value[-\s]+//;
	
	
	
	print OUT2"$fragment[0]";
			
	}
	
}
			