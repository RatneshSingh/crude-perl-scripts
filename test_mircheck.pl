#!/usr/bin/perl -w
use lib "/home/rhoades/_miRcheck"; # edit to directory that contains miRcheck.pm 
use miRcheck;
$hitfile =$ARGV[0]; # list of sequences and genomic coorinates

#parameters for miRNA checking
$WIN = 240; #nt added to each side of miRNA
$MAX_UNPAIR = 6; #max # unpaired bases in putative mir
$MAX_STAR_UNPAIR = 6;  #max # unpaired bases in putative mir*
$MAX_SIZEDIFFERENCE = 3; # max size difference between mir and mir*
$MAX_MIR_GAP = 2; # longest acceptable run of unpaired bases in mir
$MAX_STAR_GAP = 3; # longest acceptable run of unpaired bases in mir*
$MIN_FBACK_SIZE = 60; # shortest acceptable length of folback including mir and mir*
$MAX_MIR_AS_BULGE = 1; # maximum total # assymetrically unpaired bases in mir
$MIN_UNPAIR = 1; # minimum number of unpair bases in acceptable extended mirs/mir*s
$BP_EXTENSION = 3; # number of nt to extend mirs and mir*s

while (@ARGV){
    $thisarg = shift @ARGV;
    if ($thisarg eq "-win") {$WIN=shift @ARGV;}
    if ($thisarg eq "-unpair" ) {$MAX_UNPAIR=shift @ARGV;} 
    if ($thisarg eq "-star_unpair" ) {$MAX_STAR_UNPAIR=shift @ARGV;}
    if ($thisarg eq "-size_diff" ) {$MAX_SIZEDIFFERENCE=shift @ARGV;}
    if ($thisarg eq "-mir_bulge") {$MAX_MIR_GAP=shift @ARGV;}
    if ($thisarg eq "-star_bulge") {$MAX_STAR_GAP=shift @ARGV;}
    if ($thisarg eq "-fback_min") {$MIN_FBACK_SIZE=shift @ARGV;}
    if ($thisarg eq "-ass") {$MAX_MIR_AS_BULGE=shift @ARGV;}
    if ($thisarg eq "-min_unpair") {$MIN_UNPAIR=shift @ARGV;}
    if ($thisarg eq "-bp_ext") {$BP_EXTENSION=shift @ARGV;}
}

open(F,$hitfile);
while(<F>){
    if (/(\S+)\s+(\S+)\s+(\S+)\s+(\d+)/){ 
     # assumes input file contains an ID, a "genomic" sequence, a potential miRNA sequence (a substring of the "genomic" sequence,
     # and the 1st coordinate of potential miRNA within	
	$ID = $1;
	$genomic_seq = $2;
	$miR_seq = $3;
	$miR_start = $4;
	$miR_stop = $miR_start + length($miR_seq)-1;
	$filename = substr($ID,0,12);
	($fold,$E) = ();
	($fold,$E) = RNAfold("$ID",$genomic_seq,1);
	if (not($filename eq $ID)) {system("mv $filename\_ss.ps $ID\_ss.ps");} #rename postscript file, if name truncated by RNAfold

	 #judge quality of foldback
        $fback = '';
        ($fback,$fback_start,$fback_stop)  = miR_check($fold,$miR_start,$miR_stop,
						       "-unpair",$MAX_UNPAIR,
						       "-star_unpair",$MAX_STAR_UNPAIR,
						       "-size_diff",$MAX_SIZEDIFFERENCE,
						       "-mir_bulge",$MAX_MIR_GAP,
						       "-star_bulge",$MAX_STAR_GAP,
						       "-fback_min",$MIN_FBACK_SIZE,
						       "-ass",$MAX_MIR_AS_BULGE,
						       "-min_unpair",$MIN_UNPAIR,
						       "-bp_ext",$BP_EXTENSION
						       );
	print "$ID $fback $fback_start $fback_stop\n";
    }
}

	
	
