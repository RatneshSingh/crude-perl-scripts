#!/usr/bin/perl -w
use strict;
use Getopt::std;
use RS_Subs;



our($opt_s,$opt_w,$opt_m,$opt_n,$opt_c,$opt_o,$opt_g,$opt_f);
getopt('swmncogf');
my$usage='usage:perl script options.....
options:
   -s	sequence  seqall    Nucleotide sequence(s) filename and optional
                            format, or reference (input USA)
   -w	window    integer   [100] The percentage CG content and the
                            Observed frequency of CG is calculated
                            within a window whose size is set by this
                            parameter. The window is moved down the
                            sequence and these statistics are calculated
                            at each position that the window is moved
                            to. (Integer 1 or more)
   -m	minlen    integer   [200] This sets the minimum length that a
                            CpG island has to be before it is reported.
                            (Integer 1 or more)
   -n	minoe     float     [0.6] This sets the minimum average observed
                            to expected ratio of C plus G to CpG in a
                            set of 10 windows that are required before a
                            CpG island is reported. (Number from 0.000
                            to 10.000)
   -c	minpc     float     [50.] This sets the minimum average
                            percentage of G plus C a set of 10 windows
                            that are required before a CpG island is
                            reported. (Number from 0.000 to 100.000)
  -o	outfile   outfile   [*.cpgplot] This sets the name of the file
                            holding the report of the input sequence
                            name, CpG island parameters and the output
                            details of any CpG islands that are found.
  -g	graph     xygraph   [$EMBOSS_GRAPHICS value, or x11] Graph type
                            (ps, hpgl, hp7470, hp7580, meta, cps, x11,
                            tek, tekt, none, data, xterm)
  -f	outfeat   featout   [unknown.gff] File for output features

';
print $usage;
if(!defined$opt_s){print"Enter the file containing sequenc(s)\n";$opt_s=<STDIN>;}
if(!defined$opt_w){$opt_w=100;}
if(!defined$opt_m){$opt_m=200;}
if(!defined$opt_n){$opt_n=0.6;}
if(!defined$opt_c){$opt_c=50;}
#if(!defined$opt_o){$opt_o=100;}
if(!defined$opt_g){$opt_g=q(ps);}
#if(!defined$opt_f){$opt_f=100;}

my%seq= ReadFasta($opt_s);

open TEMP,">temp_seq.txt";

foreach $seq_header(keys %seq){

	print TEMP">$seq_header\n$seq{$seq_header}";
	$opt_o="$seq_header.fasta";
	$opt_f="$seq_header.gff";
	my$run_cpgplot="cpgplot -sequence temp_seq.txt -window $opt_w -minlen -minoe $opt_n -minpc $opt_c -outfile $opt_o -graph $opt_g -outfeat $opt_f";
	`$run_cpgplot`;
	
}

exit;