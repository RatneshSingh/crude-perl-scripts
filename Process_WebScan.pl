#/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

getopt('no');

my($opt_n,$opt_o,$n);

$n=$opt_n if defined $opt_n;
if(!defined $opt_n){print "number of files to be processed:\n";$n=<STDIN>;}

if($opt_o){open OUT, ">$opt_o";}
else{open OUT, ">formated_motifs.txt";}



#------------------------------------------------------------------------
#read file to seperate name sequence and patterns.
my(@filenames,@seqnames,%count,%motifs);
for(my$i=1;$i<=$n;$i++){
	my$scanfile=$i.'.txt';
	open SCANFILE,"$scanfile" or die "Cannot fine $scanfile\n";
	$/="\n________________________";
	my @section=();
	while(<SCANFILE>){push(@section,$_);}
	close SCANFILE; # close file to reduce stress on computer.
	push(@filenames,$scanfile);
	print "File being read:$scanfile\n";
	#------------------------------------------------------------------------

	#extract sequence name from file.
	my@sequenceheader=split(/\n/,$section[0]);
	my$seqname=();
	foreach(@sequenceheader){if(/>/){$_=~s/>//;($seqname)=split(/\s+/,$_);$seqname=~s/\s*//;last;}}
	push(@seqnames,$seqname);
	#print "Sequence Name:$seqname\n";

	#count the number of motifs and create motif-file specific hash

	# remove leading and trailing ends away from pattern.
	$section[1]=~s/_+||-----+[\w\s\W\S]+//g;
	#print "Motif Hits:\n$section[1]\n";
	my@motif_line=split(/\n/,$section[1]);

	foreach(@motif_line){
		next if $_=~/^\s*$/;
		$_=~s/^\s+//;
		my@motif_info=();
		@motif_info=split(/\s+/,$_);
		#print "Motif info:@motif_info\n";
		$motif_info[0]=~s/\s+//;
		my$motif_name=$motif_info[0];
		#print "Motif name:$motif_name\n";
		next if ($motif_name=~/^\s*$/);
		$count{$motif_name}{$scanfile}++;# if ($motif_name=~/$motif_info[0]/);
		#$count{$motif_name}{$seqname}++ if ($motif_name=~/$motif_info[0]/); # use this for real seq name rather than filenames.
		$motifs{$motif_name}=();
		#print "$count{$motif_name}{$scanfile}\n";	
	}
}

# print to check if it counted the motifs.

print "Writing Cis-elements summary in file formated_motifs.txt\n";
my @header=();
push(@header,keys %motifs);
my $header=join("\t",@header);
print OUT"\t$header";
foreach my $files(@filenames){
	print OUT"\n$files";
	foreach my$motifs(@header){
		$count{$motifs}{$files}=0 if !defined $count{$motifs}{$files};
		print OUT"\t$count{$motifs}{$files}";}
}

exit;







