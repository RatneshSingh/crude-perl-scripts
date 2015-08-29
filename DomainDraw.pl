#! /usr/bin/perl -w
use strict;
use GD::Simple;
#use Math::Round;
#use Bio::Graphics;
#use Bio::SearchIO;
my(%namesofseqs);
	
	my$maxbarlen=8000;
	my$maxbarwidth=50;
	my$marginx=50;
	my$marginy=500;
	my$distancebetweenbars=-15;
	my$domainraiseup=0;
	my$domainraisedown=0;
	my$maxseqlength=8070272;
	my$pixelperbase=$maxbarlen/$maxseqlength;
	my$distancefromtop=0;
	#$baseperpixel=$maxseqlength/$maxbarlen;
	
	#determine the length of the first molecule.
	
	#$moleculelength=$pixelperbase*$maxseqlength;
	
	
open FILE,"$ARGV[0]";
while(<FILE>){
	next if $_=~/^\s*$/;
	(my$nameseq,my$start,my$end)=split(/\s+/,$_);
	chomp($nameseq);
	$nameseq=~s/\s+//g;
	if($nameseq eq ''){next;}else{$namesofseqs{$nameseq}{'name'}=$nameseq;}
	#$namesofseqs{$nameseq}{}
	
	
}	

close FILE;
my$num_of_bars=scalar(keys%namesofseqs);
my$graphics_width=($num_of_bars*$maxbarwidth)+(($num_of_bars-1)*$distancebetweenbars)+($distancefromtop*2)+($marginy*2);
my$graphic_length=$maxbarlen+($marginx*2);
# create a new image
my$img2= GD::Simple->new($graphic_length,$graphics_width);



#print length %namesofseqs;
#my$distancebetweenbars=50;

foreach my$value (keys %namesofseqs){
	print "keys: $value\n";
	
	$img2->bgcolor('green');
	$img2->fgcolor('black');
	$img2->rectangle($marginx,$distancefromtop+$marginy,$maxbarlen+$marginx,$marginy+$maxbarwidth+$distancefromtop);
	
	
	$namesofseqs{$value}{'y1'}=$marginy+$distancefromtop-$domainraiseup;
	$namesofseqs{$value}{'y2'}=$marginy+$maxbarwidth+$distancefromtop+$domainraisedown;
	
	print $namesofseqs{$value}{'y1'}."\n";
	print $namesofseqs{$value}{'y2'}."\n";

	
	
	$distancefromtop=$distancefromtop+$distancebetweenbars+$maxbarwidth;
	
}








open FILE,"$ARGV[0]";
print "\n Drawing coordinates from file $ARGV[0]\n ";
while(<FILE>){
	
	(my$name,my$start,my$end)=split(/\s+/,$_);
	chomp($name);
	$start=~s/\D+//g;
	$end=~s/\D+//g;
	my$positionx1=int(($start*$pixelperbase)+0.5)+$marginx;
	#print $positionx1;
	my$positiony1= $namesofseqs{$name}{'y1'};#7;
	my$positionx2= int(($end*$pixelperbase)+0.5)+$marginx;
	my$positiony2= $namesofseqs{$name}{'y2'};#53;
	#$seqname=$namesofseqs{$name};
		#draw red rectangles on the gray bars
		$img2->bgcolor('red');
		$img2->fgcolor(undef);
		$img2->rectangle($positionx1,$positiony1,$positionx2,$positiony2);
	
	
}

print "\n Finished drawing coordinates from file $ARGV[0]\n ";

	
# Draw gene on the the top of CpG island bars

my$above=5; # distance of gene above CpGI bar. 
my$color='blue';
my$direction='arrow';
$img2->penSize(1);




open FILE2,"$ARGV[1]";

print "\nDrawing coordinates from file $ARGV[1]\n ";
while(<FILE2>){
	# ChrX-pseudomolecule	PYhCpXYh18_69D24	Exon 1	3472866-3473204	minus
	#print "$_\n";
	#$_=~s/"Exon "/"Exon"/;
	(my$Genomic,my$gene,my$exon,my$coordinates,my$strand)=split(/\t+/,$_);
	#print "Genomic:$Genomic\tGene:$gene\tExon:$exon\tCoord:$coordinates\tStrand:$strand\n";
	my($start,$end)=split(/-/,$coordinates);

	chomp($Genomic);
	#print "Coordinates\t$coordinates\tstart:$start\tEnd:$end\n";
	$start=~s/\D+//g;
	$end=~s/\D+//g;
	my($positionx1,$positionx2,$positiony1,$positiony2,$old_positionx1,$old_positionx2,$old_positiony1);
	
	
	if($exon=~/Exon/){
		$positionx1=int(($start*$pixelperbase)+0.5)+$marginx;
		#print $positionx1;
		$positiony1= $namesofseqs{$Genomic}{'y1'}-$above;#7;
		$positionx2= int(($end*$pixelperbase)+0.5)+$marginx;
		#$positiony2= $namesofseqs{$Genomic}{'y2'}-$above-5;#53;
		$color='red';

	}
	
	if($exon=~/gene/){
		
		#print $positionx1;
		$positionx1= int(($start*$pixelperbase)+0.5)+$marginx;
		$positiony1= $namesofseqs{$Genomic}{'y1'}-$above-2;#7;
		$positionx2= int(($end*$pixelperbase)+0.5)+$marginx;
		#$positiony2= $namesofseqs{$Genomic}{'y2'}-$above-5;#53;
		$color='blue';
			
		#$old_positionx1=$positionx1 if !defined $old_positionx1;
#		$old_positionx2=$positionx2 if !defined $old_positionx2;
#		$old_positiony1=$positiony1 if !defined $old_positiony1;
#
#
#		if($positionx1==$old_positionx2){
#			$positionx1=int(($start*$pixelperbase)+0.5)+$marginx;
#			$positiony1= $old_positiony1-2;#7;
#			$positionx2= int(($end*$pixelperbase)+0.5)+$marginx;
#			
#			$old_positionx1=$positionx1;
#			$old_positionx2=$positionx2;
#			$old_positiony1=$positiony1;
#	
#		}
#
#		$old_positionx1=$positionx1;
#		$old_positionx2=$positionx2;
#		$old_positiony1=$positiony1;
	}
	

	#$seqname=$namesofseqs{$name};
		#draw red rectangles on the gray bars
		$img2->bgcolor(undef);
		$img2->fgcolor($color);
		#$img2->penSize(1,1);
		#$img2->rectangle($positionx1,$positiony1,$positionx2,$positiony2);
		
		
		
		
		if($exon eq 'gene'){$img2->bgcolor(undef);
		$img2->fgcolor('black');
		$img2->moveTo($positionx1-2,$positiony1-2);
    	#$img2->font('Times:italic');
    	$img2->fontsize(20);
    	$img2->angle(-90);
		$img2->string($gene);
		}
		$img2->bgcolor(undef);
		$img2->fgcolor($color);

		$img2->angle(0);
		$img2->moveTo($positionx1,$positiony1);
    	$img2->lineTo($positionx2,$positiony1);
    		
    		print "\nDrawing $exon from $positionx1 to $positionx2";
	
}
 
 print "\n Finished drawing coordinates from file $ARGV[1]";

 
 
 
 
    
   
   open OUT2,">test2.png";
   print OUT2 $img2->png(0);
   close OUT2;
   open OUT2,">test2.jpg";
   print OUT2 $img2->jpeg([100]);
   close OUT2;

   
