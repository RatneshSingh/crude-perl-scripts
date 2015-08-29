#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

our($opt_s,$opt_o);
getopt('so');
my $usage = "Usage $0 -s <infile.embl> -o <outfile.fa>\n";

my $infile = $opt_s or die "No infile name given\n$usage\n";
#my$infile="test_embl2fasta.embl";
my $outfile =$opt_o or die "No outfile $opt_s name given\n$usage\n";
$/='//';
open (IN,'<',$infile) || die "unable to open $infile: $!\n";
open(OUT,'>',$outfile) || die "Unable to open $outfile: $!\n";

print "\n
####################################################
Converting embl file in RepeatMasker Fasta format\n
####################################################
\n
";
while(<IN>){
   chomp;
   $_=~s/\/\///g;
   next if $_=~/^\s*$/;
   #print "\n\n********Sending New entry to subroutine******** \n$_";
   my($header,$seq)=embl2fasta_repeatmasker(\$_);
   #print "\n\n********Recieved from subroutine******** \nHeader:$header\n\nseq:$seq\n End of New entry**********\n\n";
   print OUT">$header\n$seq\n";
     
}
close(IN);
close(OUT);




sub embl2fasta_repeatmasker{
   my$aRef_embl=shift;
   my@line=split(/\n/,$$aRef_embl);
   
   my($current_id,$current_species,$current_type,$current_subtype,$current_searchstages,$current_bufferstages,$seq,$current_RepbaseID)=" ";
   my$count=0;
   my$newcount=0;
   foreach my$lines(@line){
      $count++;
      next if $lines=~/^\s+$/;
      $lines=~s/^\s+//;
      if($lines=~/^ID\s+([\.\(\)\w\d\-\_]+)\s+/){$current_id = $1;chomp $current_id}
      if($lines=~/^CC\s+Species:\s*([\w\d\.\-\,]+)/){$current_species = $1;chomp $current_species}
      if($lines=~/^CC\s+Type:\s*([\w\d\?\_]+)/){$current_type = $1;chomp $current_type}
      if($lines=~/^CC\s+SubType:\s*([\w\d\?\_\-]+)/){$current_subtype = $1;chomp $current_subtype}
      if($lines=~/^CC\s+SearchStages:\s*([\d\,]+)/){$current_searchstages = $1;chomp $current_searchstages}
      if($lines=~/^CC\s+BufferStages:\s*([\d\,]+)/){$current_bufferstages = $1;chomp $current_bufferstages}
      if($lines=~/^DE\s+RepbaseID:\s*([\w\W]+)/){$current_RepbaseID = $1;chomp $current_RepbaseID}
      if($lines=~/^SQ\s+/){$newcount=$count }
   }
   
   if($newcount!=0){for(my$i=$newcount;$i<scalar@line;$i++){$seq.=$line[$i];}}
   else {$seq='NNNNN'}
   $seq=~s/\W+|\d+|\s+//g;
   my$header=$current_id.
                   '#'.
                   type_subtype($current_type,$current_subtype).
                   '@'.$current_species.
                   research_stages($current_searchstages).
                   Repbase($current_RepbaseID);
   return($header,$seq);
   
   
}


sub type_subtype{
   my($type,$subtype)=@_;
   if(!$subtype||$subtype=~/^\s*$/){return "$type "}
   else{return "$type/$subtype "}
    
}

sub research_stages{
 my$resstag=shift;
 if(!$resstag || $resstag=~/^\s*$/){return "  [S:]"}
 else {return "  [S:$resstag]"}
}

sub Repbase{
   my$Repbase=shift;
   if(!$Repbase||$Repbase=~/^\s*$/){return ""}
 else {return " RepbaseID: $Repbase"}
   
}
