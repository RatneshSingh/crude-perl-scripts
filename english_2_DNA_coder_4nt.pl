#!/usr/bin/perl

# Created 03/10/2010
# Author Ratnesh Singh
#
use warnings;
use strict;
use Getopt::Std;
our($opt_c,$opt_d,$opt_f,$opt_o);
getopt('cfdo');
my(%ABC_2_DNA,%DNA_2_ABC);


my $usage="\n\n
==========================================================
This program can be used for coding text 
in to a DNA and decoding a DNA sequence 
back to original text.
----------------------------------------------------------
usage: this_program.exe -c/-d True -f input_file -o output
where:-
-c	code text into DNA.
-d	decode DNA into text
-f	Input file containing text/DNA seq
-o	output file to store result
==========================================================
\n\n";

if($opt_c){ABC_2_DNA();}
elsif($opt_d){DNA_2_ABC();}
else{if(!$opt_c &&!$opt_d){print "please provide -c or -d option on the command line\n$usage"; die;}
}
if($opt_f){open FILE,"$opt_f" or die "Cannot open file\n$usage";}
else{print "Provide file to use with -f option\n$usage";die;}
my$DNA_line=();

open OUT,">$opt_o" if defined $opt_o;

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

	if($opt_o){print OUT"$DNA_line\n";}
	else{print "$DNA_line\n";}
}

else{
	my @whole_line=();
	while(<FILE>){push(@whole_line,$_);}
	
		
		my$line=join("",@whole_line);
		$line=~s/\s+//g;
		for(my$i=0;$i<=length($line)-1;$i=$i+4){
			my$char=substr($line,$i,4);
			my$code= decode($char);
			$DNA_line.=$code;
			$code=();	
			}
	
	if($opt_o){print OUT"$DNA_line\n";}
	else{print "$DNA_line\n";}

}
	

####################################################################

sub code{
	my$char=shift@_;
	#chomp $char;
	#print "$char is assigned: $ABC_2_DNA{$char}\n";

	if($ABC_2_DNA{$char}){#print "$char is assigned: $ABC_2_DNA{$char}\n";
	return $ABC_2_DNA{$char}; 
	}
	elsif($char=~/\t/){#print "$char is assigned:TTAG\n";
	return 'TTAG';
	}
	elsif($char=~/^$/){#print "$char is assigned:TACG\n";
	return 'TACG';
	}
	#elsif($char=~/\t/){print "$char is assigned:GAAC\n";;return('GAAC');}
	elsif($char=~/\\/){#print "$char is assigned:CCGA\n";
	return 'CCGA';
	}
	elsif($char=~/0/){#print "$char is assigned:CCGA\n";
	return 'GTGA';
	}
	else{print "No code is assigned for $char\n";
	return 'NNN';
	}
	
	}

sub decode{
	my$char1=shift@_;
	#chomp $char1;
	#print "$char1 is assigned: $DNA_2_ABC{$char1}\n";
	if($DNA_2_ABC{$char1}){return $DNA_2_ABC{$char1}; }
	elsif($char1=~/TTAG/){return "\t";}
	#elsif($char1=~/GAAC/){return "\t";}
	elsif($char1=~/TACG/){return "\n";}
	elsif($char1=~/CCGA/){return qq(\\);}
	elsif($char1=~/GTGA/){return "0";}
	else{print "No code is assigned for $char1\n";return " ";}
	
}

sub ABC_2_DNA{

%ABC_2_DNA=(
'A'=>'AGAA',
'B'=>'CAAT',
'C'=>'AACG',
'D'=>'CCCG',
'E'=>'CTGT',
'F'=>'GGCA',
'G'=>'TGAG',
'H'=>'TCAT',
'I'=>'CGTG',
'J'=>'ACAC',
'K'=>'AAAT',
'L'=>'ATAA',
'M'=>'CCGG',
'N'=>'TGTC',
'O'=>'TTTT',
'P'=>'TCTG',
'Q'=>'TAAA',
'R'=>'AATA',
'S'=>'CATG',
'T'=>'CGCC',
'U'=>'GCAT',
'V'=>'ACCG',
'W'=>'TTCT',
'X'=>'GTGG',
'Y'=>'ATGC',
'Z'=>'GAAT',
" "=>'AAGA',
','=>'ACTC',
'_'=>'AGCT',
';'=>'AGGG',
'='=>'AGTT',
'+'=>'ATCG',
'-'=>'ATTC',
'?'=>'CACC',
'.'=>'CAGA',
'('=>'CCAC',
')'=>'CCTG',
'z'=>'CGAG',
'k'=>'CGGT',
'q'=>'CTAT',
'w'=>'CTCG',
'd'=>'CTTC',
'o'=>'GACA',
'r'=>'GAGC',
'a'=>'GATG',
'p'=>'GCGT',
'e'=>'GCTA',
'm'=>'GCCT',
'c'=>'GGAG',
's'=>'GGGC',
'j'=>'GGTA',
't'=>'GTAT',
'f'=>'GTCT',
'n'=>'GTTG',
'b'=>'TACA',
'x'=>'TAGC',
'i'=>'TATT',
'v'=>'TCGA',
'g'=>'TCCG',
'y'=>'TGCC',
'l'=>'TGGT',
'u'=>'TTAA',
'h'=>'TTGG',
'{'=>'TTGC',
'}'=>'CTGA',
'`'=>'GATC',
'~'=>'CTTG',
'!'=>'GGAT',
'|'=>'AATC',
q(\\)=>'CCGA',
'>'=>'ACGT',
'<'=>'TGCG',
'*'=>'GGGG',
'&'=>'CCCC',
'%'=>'AAAA',
'#'=>'ACTG',
'['=>'GTAG',
']'=>'GTGT',
':'=>'TCGT',
'1'=>'GACG',
'2'=>'GAGT',
'3'=>'GACT',
'4'=>'CACA',
'5'=>'ATGA',
'6'=>'TGCA',
'7'=>'AGTA',
'8'=>'TATA',
'9'=>'CGTA',
#0=>'GTGA',
q(/)=>'GGCT',
q($)=>'ATGT',
q(')=>'GCTG',
q(")=>'GCGC',
q(@)=>'AGTC',
q(^)=>'GACC',
''=>'TATG');


}

sub DNA_2_ABC{


%DNA_2_ABC=(
'AGAA'=>'A',
'CAAT'=>'B',
'AACG'=>'C',
'CCCG'=>'D',
'CTGT'=>'E',
'GGCA'=>'F',
'TGAG'=>'G',
'TCAT'=>'H',
'CGTG'=>'I',
'ACAC'=>'J',
'AAAT'=>'K',
'ATAA'=>'L',
'CCGG'=>'M',
'TGTC'=>'N',
'TTTT'=>'O',
'TCTG'=>'P',
'TAAA'=>'Q',
'AATA'=>'R',
'CATG'=>'S',
'CGCC'=>'T',
'GCAT'=>'U',
'ACCG'=>'V',
'TTCT'=>'W',
'GTGG'=>'X',
'ATGC'=>'Y',
'GAAT'=>'Z',
'AAGA'=>" ",
'ACTC'=>',',
'AGCT'=>'_',
'AGGG'=>';',
'AGTT'=>'=',
'ATCG'=>'+',
'ATTC'=>'-',
'CACC'=>'?',
'CAGA'=>'.',
'CCAC'=>'(',
'CCTG'=>')',
'CGAG'=>'z',
'CGGT'=>'k',
'CTAT'=>'q',
'CTCG'=>'w',
'CTTC'=>'d',
'GACA'=>'o',
'GAGC'=>'r',
'GATG'=>'a',
'GCGT'=>'p',
'GCTA'=>'e',
'GCCT'=>'m',
'GGAG'=>'c',
'GGGC'=>'s',
'GGTA'=>'j',
'GTAT'=>'t',
'GTCT'=>'f',
'GTTG'=>'n',
'TACA'=>'b',
'TAGC'=>'x',
'TATT'=>'i',
'TCGA'=>'v',
'TCCG'=>'g',
'TGCC'=>'y',
'TGGT'=>'l',
'TTAA'=>'u',
'TTGG'=>'h',
'TTGC'=>'{',
'CTGA'=>'}',
'GATC'=>'`',
'CTTG'=>'~',
'GGAT'=>'!',
'AATC'=>'|',
'CCGA'=>q(\\),
'ACGT'=>'>',
'TGCG'=>'<',
'GGGG'=>'*',
'CCCC'=>'&',
'AAAA'=>'%',
'ACTG'=>'#',
'GTAG'=>'[',
'GTGT'=>']',
'TCGT'=>':',
'GACG'=>'1',
'GAGT'=>'2',
'GACT'=>'3',
'CACA'=>'4',
'ATGA'=>'5',
'TGCA'=>'6',
'AGTA'=>'7',
'TATA'=>'8',
'CGTA'=>'9',
#'GTGA'=>0,
'GGCT'=>q(/),
'ATGT'=>q($),
'GCTG'=>q('),
'GCGC'=>q("),
'AGTC'=>q(@),
'GACC'=>q(^),
'TATG'=>'');
	}