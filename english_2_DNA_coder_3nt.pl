#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
our($opt_c,$opt_d,$opt_f,$opt_e);
getopt('cfde');
my(%ABC_2_DNA,%DNA_2_ABC);

if($opt_c){ABC_2_DNA();}
if($opt_d){DNA_2_ABC();}
if(!$opt_c &&!$opt_d){print "please provide -c or -d option on the command line\n"; die;}

if($opt_f){open FILE,"$opt_f" or die "Cannot open file\n";}
else{print "Provide file to use with -f option\n";die;}
my$DNA_line=();


if($opt_c){
	while(<FILE>){
		my$line=$_;
		my@line=split(//,$_);
		foreach my$char(@line){
		my$code= code($char);
		$DNA_line.=$code;
		$code=();	
		}
	}

	print "$DNA_line\n";
}

else{

	while(<FILE>){
		my$line=$_;
		for(my$i=0;$i<=length($line)-1;$i=$i+3){
			my$char=substr($line,$i,3);
			my$code= decode($char);
			$DNA_line.=$code;
			$code=();	
			}
	}
	print "$DNA_line\n";

}
	

####################################################################

sub code{
	my$char=shift@_;
	chomp $char;
	if($ABC_2_DNA{$char}){return $ABC_2_DNA{$char}; }
	elsif($char=~/\s*$/){return 'ACG';}
	else{print "No code is assigned for $char\n";return 'NNN';}
}

sub decode{
	my$char1=shift@_;
	chomp $char1;
	if($DNA_2_ABC{$char1}){return $DNA_2_ABC{$char1}; }
	elsif($char1=~/'ACG'/){return "\n";}
	else{return "\n"; print "No code is assigned for $char1\n";}
}

sub ABC_2_DNA{

%ABC_2_DNA=(
'A'=>'AGA',
'B'=>'CAA',
'C'=>'AAC',
'D'=>'CCC',
'E'=>'CTG',
'F'=>'GGC',
'G'=>'TGA',
'H'=>'TCA',
'I'=>'CGT',
'J'=>'ACA',
'K'=>'AAA',
'L'=>'ATA',
'M'=>'CCG',
'N'=>'TGT',
'O'=>'TTT',
'P'=>'TCT',
'Q'=>'TAA',
'R'=>'AAT',
'S'=>'CAT',
'T'=>'CGC',
'U'=>'GCA',
'V'=>'ACC',
'W'=>'TTC',
'X'=>'GTG',
'Y'=>'ATG',
'Z'=>'GAA',
' '=>'AAG',
','=>'ACT',
'_'=>'AGC',
';'=>'AGG',
'='=>'AGT',
'+'=>'ATC',
'-'=>'ATT',
'?'=>'CAC',
'.'=>'CAG',
'('=>'CCA',
')'=>'CCT',
'z'=>'CGA',
'k'=>'CGG',
'q'=>'CTA',
'w'=>'CTC',
'd'=>'CTT',
'o'=>'GAC',
'r'=>'GAG',
'a'=>'GAT',
'p'=>'GCG',
'e'=>'GCT',
'm'=>'GCC',
'c'=>'GGA',
's'=>'GGG',
'j'=>'GGT',
't'=>'GTA',
'f'=>'GTC',
'n'=>'GTT',
'b'=>'TAC',
'x'=>'TAG',
'i'=>'TAT',
'v'=>'TCG',
'g'=>'TCC',
'y'=>'TGC',
'l'=>'TGG',
'u'=>'TTA',
'h'=>'TTG');


}

sub DNA_2_ABC{


%DNA_2_ABC=(
'AGA'=>'A',
'CAA'=>'B',
'AAC'=>'C',
'CCC'=>'D',
'CTG'=>'E',
'GGC'=>'F',
'TGA'=>'G',
'TCA'=>'H',
'CGT'=>'I',
'ACA'=>'J',
'AAA'=>'K',
'ATA'=>'L',
'CCG'=>'M',
'TGT'=>'N',
'TTT'=>'O',
'TCT'=>'P',
'TAA'=>'Q',
'AAT'=>'R',
'CAT'=>'S',
'CGC'=>'T',
'GCA'=>'U',
'ACC'=>'V',
'TTC'=>'W',
'GTG'=>'X',
'ATG'=>'Y',
'GAA'=>'Z',
'AAG'=>' ',
'ACT'=>',',
'AGC'=>'_',
'AGG'=>';',
'AGT'=>'=',
'ATC'=>'+',
'ATT'=>'-',
'CAC'=>'?',
'CAG'=>'.',
'CCA'=>'(',
'CCT'=>')',
'CGA'=>'z',
'CGG'=>'k',
'CTA'=>'q',
'CTC'=>'w',
'CTT'=>'d',
'GAC'=>'o',
'GAG'=>'r',
'GAT'=>'a',
'GCG'=>'p',
'GCT'=>'e',
'GCC'=>'m',
'GGA'=>'c',
'GGG'=>'s',
'GGT'=>'j',
'GTA'=>'t',
'GTC'=>'f',
'GTT'=>'n',
'TAC'=>'b',
'TAG'=>'x',
'TAT'=>'i',
'TCG'=>'v',
'TCC'=>'g',
'TGC'=>'y',
'TGG'=>'l',
'TTA'=>'u',
'TTG'=>'h');








}