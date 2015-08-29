
#! /usr/bin/perl -w;
use strict;

$/="\n\/\/";
while(<>){
	my$fragment=();
	my$host=();
	my$country=();
	my$sequence=();
	my$id=();
	if(/(FEATURES[\w\W]+ORIGIN)/){$fragment=$1;}
	
	if($fragment=~/host=([\"\w\d\s-]+)\n/){$host=$1;}
	elsif($fragment=~/isolate=([\"\w\d\s-]+)\n/){$host=$1;}
	elsif($fragment=~/cultivar=([\"\w\d\s-]+)\n/){$host=$1;}
	if($fragment=~/country=([\"\w\d\:\s\-\,]+)\n/){$country=$1;}
	if($fragment=~/translation=([\"\w\s]+)\n/){
		$sequence=$1;
		$sequence=~s/\s//g;
		
	}
	if($fragment=~/protein_id=([\"\w\d\.]+)\n/){$id=$1;}
	
	print "$id\t$host\t$country\t$sequence\n";	
}