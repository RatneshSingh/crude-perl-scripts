#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

our($qseq,$qtable,@target_fasta,@blast_tables,$best_by_col,$best_by_min,$help,$blast,$dbtype,%nam_hash,%nam_hash_val,$delim,$col,%blastinfo,$diamond);

my$border="\n#########################################\n";
my$sborder="\n*****************************************\n";
my $usage= <<EOF;

This script takes either fasta sequences to blast against defined multiple databases and
create a mapping table for each query seq to database seq. Alternatively it can take
multiple blast result in table format and map it to the first column names in qtable to create a
mapped qtable. qtable could be Expression ratio table with names in first column and
expression data in other columns. Mixing of sequence and blast table is allowed as long
as they are provided by appropriate flags.

Usage: $0  options...
options:
  -qseq|qs     Query sequence in fasta format if you want to run blast with other databases.
  -qtable|qt   Table (eg. expression ratio)  in which to add new names.
              First column is used as query name to map onto other tables/databases.
  -tseq|ts     Target sequence(s) to be used as blast database to blast -qseq sequences.
  -ttable|tt   Target blast result(s) from blast against databases.
  -bcol|bc     To Pick best hit from blast result, use this column.
  -bmin|bmn    best hit is with minimum value in column -bcol  [max is best hit].
              If not provided, ttable is assumed to contain top hits only.
  -blast|bt    Blasttype to use if blast run is requested.
              e.g:
	      blastn   : Q: Nuc -> Db: Nuc
	      blastx   : Q: Nuc -> Db: Prot
	      blastp   : Q: Prot -> Db: Prot
	      tblastn   : Q: Prot-> Db: Nuc
	      tblastx   : Q: Nuc -> Db: Nuc
  -diamond      Use diamond instead of blast to run database search.
                currently only blastx or blastp is supported. Make sure db is protein
  -delim 		delimiter to split Query name and col. overrides -blast
  -col 		column to use as name for matching after splitting the name.
  -help|h      Print this help and exit.



EOF



Getopt::Long::GetOptions(
  "qseq|qs=s"=>\$qseq,
  "qtable|qt=s"=>\$qtable,
  "tseq|ts=s"=>\@target_fasta,
  "ttable|tt=s"=>\@blast_tables,
  "bcol|bc=i"=>\$best_by_col,  ### if to use filtering, use this column.
  "bmin|bmn"=>\$best_by_min,   ### best hit is with minimum value in column $best_by_col
  "blast|bt=s"=>\$blast,
  "diamond|dia"=>\$diamond,
  "delim|d=s"=>\$delim,
  "col|c=i"=>\$col,
  "help|h"=>\$help
);


die "\n$usage\n" if $help;
die "\n No options provided. See Usage instructions\n$usage\n" if (!$qseq && !$qtable);
## process if options are provided as commaseperated list.
@target_fasta = split(/,/,join(',',@target_fasta)) if @target_fasta;
@blast_tables = split(/,/,join(',',@blast_tables)) if @blast_tables;
my $qtype = seqtype($qseq) if $qseq;


### run blast if fasta sequences are provided.
if (scalar@target_fasta > 0) {
  foreach my$sfile(@target_fasta){
    ## find the type of blast needed to run if not provided by the user.
    my$stype;
    if ($blast) {
      print "\nWill run blast type: $blast"
    }else{

      ## if provided $sfile is database then use extensions to decide seqtype.
      if ( -e $sfile ) {
        $stype=seqtype($sfile) if  -T $sfile;
      }elsif(-e "$sfile.nal" || -e "$sfile.nsq" ){
        $stype="nucl";
      }elsif(-e "$sfile.pal" || -e "$sfile.psq"){
        $stype="prot";
      }elsif(-e "$sfile.dmnd"){
        $stype="prot";
      }

      $blast="blastn"  if ($qtype eq "nucl" && $stype eq "nucl");
      $blast="blastx"  if ($qtype eq "nucl" && $stype eq "prot");
      $blast="blastp"  if ($qtype eq "prot" && $stype eq "prot");
      $blast="tblastn" if ($qtype eq "prot" && $stype eq "nucl");
    }
    ## Start creating blast command to run
    my@dbname=split /\//,$sfile;
    my@qfilename=split /\./,$qseq;
    my$outfile=join ".",@qfilename[ 0 .. $#qfilename-1 ],"vs",$dbname[-1],$diamond?"diamond":$blast,"table";

    ## run database search if output file does not exists.
    if( -e $outfile ){
	      print "$sborder Result from previous $blast exists. \nRename or remove the file $outfile to rerun the blast.... Skipping blast run $sborder";
    }else{

        if ( $diamond && ($blast eq 'blastp' || $blast eq 'blastx') ) {
          besthit_diamond($blast,$qseq,$sfile,$outfile) || print "\nCould not finish diamond. Quiting.";
        }elsif( $diamond && $stype eq "nucl" ) {
          print "\n Error diamond program can only be run on protein databases. " ;
        }elsif($blast) {
          besthit_blast($blast,$qseq,$sfile,$outfile) || print "\nCould not finish blastrun. Quiting.";
        }
    }

    push(@blast_tables,$outfile) if -e $outfile;
    $blastinfo{$outfile}{db}=$dbname[-1];
  }
}







### process blast results if blast tables are provided.

if (scalar@blast_tables > 0) {
	print "$sborder Processing blast files $sborder";
  foreach my$bfile(@blast_tables){
	  print "Precessing : $bfile\n";
    open BLAST,"$bfile"  or die "Unable to open file $bfile. Quiting";
    while (<BLAST>){
      #my($query,$subject,@rest);
      #($query,$subject,@rest)=split(/\s+/,$_);
      my@rest=split(/\s+/,$_);
      $rest[0] = getname($rest[0],$delim,$col) if $delim;
      ### set current value if not assigned yet.
      $nam_hash{$rest[0]}{$bfile}     ||= $rest[1];

      ## test if values need to be changed
      if ($best_by_col){
        $nam_hash_val{$rest[0]}{$bfile} ||= $rest[$best_by_col-1];
        if ($best_by_min) {
          $nam_hash{$rest[0]}{$bfile}       = $rest[$best_by_col-1]  <  $nam_hash_val{$rest[0]}{$bfile} ? $rest[1] : $nam_hash{$rest[0]}{$bfile};
          $nam_hash_val{$rest[0]}{$bfile}   = $rest[$best_by_col-1]  <  $nam_hash_val{$rest[0]}{$bfile} ? $rest[$best_by_col-1] : $nam_hash_val{$rest[0]}{$bfile};
        }else{
          $nam_hash{$rest[0]}{$bfile}       = $rest[$best_by_col-1]  >  $nam_hash_val{$rest[0]}{$bfile} ? $rest[1] : $nam_hash{$rest[0]}{$bfile};
          $nam_hash_val{$rest[0]}{$bfile}   = $rest[$best_by_col-1]  >  $nam_hash_val{$rest[0]}{$bfile} ? $rest[$best_by_col-1] : $nam_hash_val{$rest[0]}{$bfile};
        }
      }
    }
  }
}

#print join("\t","AthHomolog","GenbankID","GenBankAccession","RefSeqID","RefSeqAccession");
print "$sborder Saving Mapping results:$sborder";
if ($qtable) {
open OUT, ">$qtable.blastmapped.table";
open RATIO,"$qtable";
  while(<RATIO>){
    my($bl,@rest)=split(/\s+/,$_);
    $bl=~s/\|[^\s]+//g;
   $bl=getname($bl,$delim,$col) if $delim;
   if ($. == 1) {
     #print OUT "$bl\t";
     print OUT join "\t", $blastinfo{$_}{db}?$blastinfo{$_}{db}:$_, ""  foreach ($bl,@blast_tables);
     print OUT join "\t",@rest,"\n";
     next;
   }

    print OUT "$bl\t";
    foreach my$fblast(@blast_tables){print OUT join "\t",$nam_hash{$bl}{$fblast}?$nam_hash{$bl}{$fblast}:"NoMatch",""}
    print OUT join "\t",@rest,"\n";
    #print "\n$bl";
  }
  print "\nMapping results are saved as: $qtable.blastmapped.table $border";
  print "$border RUN FINISHED $border";
}
else{
	open OUT, ">$qseq.blastmapped.table";
    print OUT join "\t", $blastinfo{$_}{db}?$blastinfo{$_}{db}:$_, ""  foreach ($qseq,@blast_tables);
    foreach my$qnames( keys %nam_hash ){
      print OUT "\n$qnames\t";
      foreach my$fblast1(@blast_tables){print OUT join "\t",$nam_hash{$qnames}{$fblast1}?$nam_hash{$qnames}{$fblast1}:"NoMatch","";}
    }
    print "Mapping results are saved as: $qseq.blastmapped.table";
    print "$border RUN FINISHED $border";
}




### subroutines

sub getname{
	my$name=shift;
	my$lim=shift;
	my$ncol=shift;
	$name=~s/^\s*//g;
	$name=~s/\s*$//g;
	my@frags=split /\Q$lim/,$name;
	#print "\nreturned $frags[$ncol-1]  for  $name";
	return ($frags[$ncol-1])

}

## runs blast program and returns 1 if successful.
sub besthit_blast{
  my$btype=shift;
  my$query=shift;
  my$db=shift;
  my$out=shift;

  if (-e "$db.nsq" ||  -e "$db.nal" || -e "$db.pal" || -e "$db.psq") {
    print "Blast database for $db exists.  Not creating blastdb"
  }
  else{
    make_blast_db($db)

  }

  my $cmd = join " ",$btype, ' -outfmt 6 ',' -num_threads 10 ', ' -max_target_seqs 1 ', ' -max_hsps 1 ', ' -xdrop_gap 1000 '," -query $query ", " -db $db ", " -out $out ", " -evalue 1e-5 ";

  if (system($cmd) == 0){return 1}else{return 0}

}


sub besthit_diamond{
  my$btype=shift;
  my$query=shift;
  my$db=shift;
  my$out=shift;

  ### inform user if exists result
  if ( -e $out) {
    print "\n Existing result from previous run exists. Please rename or remove to run diamond again.\n NOT RUNNING DIAMOND ON THIS SET OF QUERY AND DATABASE\n";
    return 1;
  }



  if ( $btype ne 'blastx' && $btype ne 'blastp' ){print "\nblasttype could only be blastx or blastp. \nCurrent blast type is $btype\n"; return 0;}

  ## check if database exists.
  if (-e "$db.dmnd") {
    print "\n Diamond database $db.dmnd found. Skipping database creation step.";
  }elsif (print "\nCreating diamond database\n" && make_diamond_db($db)){
    print "\nCreated diamond database.\n"
  }else{ print "\nUnable to create diamond database."; return 0; }

  ### cmd creation for running diamond and converting output to table format.
  my $cmd = join " ","diamond "," $btype ", ' --max-target-seqs 1 ', " --query $query ", " --db $db ", " --daa $out.daa ", " --evalue 0.00001 ", " --sensitive ", " --gapped-xdrop 1000 ";
  my $cmd2= join " ","diamond ", " view ", " --daa $out.daa ", " --out $out ";


  ### if esists $out.daa do not run diamond and convert it to table.
  if (-e "$out.daa") {
    print "\nResult from previous diamond run exists. \nPlease remove or rename $out.daa to run diamond again.\nConverting $out.daa to table format.\n$cmd2";
    return 1 if system($cmd2) ==0;
  }


  print "\nRunning $cmd\n";
  if (system($cmd)==0 ){
    print "\n Finished Running cmd:\n";
    if (print "\nRunning $cmd2\n" && system($cmd2)==0 ){
      print "\n Finished converting result to table format.\n";
      return 1;
    }else{ print "\n failed to convert diamond result to Table format.\n $cmd2 failed."; return 0}
  }else{
    print "\n Failed to run diamond \n $cmd failed."; return 0;
  }

}


sub make_diamond_db{
my$dbfile=shift;
return 1 if (-e "$dbfile.dmnd");
my$cmd="diamond makedb --in $dbfile -d $dbfile  ";
if (system($cmd)== 0 ){return 1} else { return 0}
}

sub make_blast_db{
my$dbfile=shift;
return 1 if (-e "$dbfile.nsq" ||  -e "$dbfile.nal" || -e "$dbfile.pal" || -e "$dbfile.psq");
my$cmd="makeblastdb -in $dbfile -dbtype ".seqtype($dbfile);
if(system($cmd)==0){return 1}else{return 0}
}



sub seqtype{
  my$seqfile=shift;
  my$samp_seq="";
  return "bin" if -B $seqfile;
  open(DBTYPE,"$seqfile") or die "\n Unable to open $seqfile for guessing sequence type\n";
  while (<DBTYPE>) {
    next if m/^>/;
    last if length($samp_seq) > 1000;
    $samp_seq.= substr($_, 0, 100 < length($_)? 100 : length($_));
 }
  close DBTYPE;

  my $nucs = $samp_seq =~ tr/ATGCNXatgcnx//;

  return "nucl" if ($nucs/length($samp_seq)) > 0.8;
  return "prot";

}
