use warnings;
use strict;
use Getopt::Long;
use File::Basename;
#use Sys::CpuAffinity;

our($fasta,$out,$args,$help,$database,$dir,@args,$rem,$outfmt);

$outfmt=0;
$args="";
my$numcpu = 20;
my$prg='est2genome';
my$option=GetOptions(
      'program|prg|p:s'=>\$prg,
      'numcpu|num_threads|n:i'      => \$numcpu,
      'sequence|query|s:s'      => \$fasta,
      'database|d|db:s'  =>\$database,
      'output|out|o:s'      => \$out,
	  'outfmt|of=s'=>\$outfmt,
      'dir|folder|f:s'      =>\$dir,
      'args|a:s'      => \@args,
      'remove|r|rt' =>\$rem,
      'help|h'       =>\$help
);
### create directory to collect blast results.
$dir=join "","parralal.","$prg.",localtime,".tmp";

my $usage="$0 -s sequencefile -db database  [options]

-p model to run in exonerate [est2genome]
    ungapped	The simplest type of model, used by default. An appropriate model with be selected automatically for the type of input sequences provided.
    ungapped:trans	This ungapped model includes translation of all frames of both the query and target sequences. This is similar to an ungapped tblastx type search.
    affine:global	This performs gapped global alignment, similar to the Needleman-Wunsch algorithm, except with affine gaps. Global alignment requires that both the sequences in their entirety are included in the alignment.
    affine:bestfit	This performs a best fit or best location alignment of the query onto the target sequence. The entire query sequence will be included in the alignment, but only the best location for its alignment on the target sequence.
    affine:local	This is local alignment with affine gaps, similar to the Smith-Waterman-Gotoh algorithm. A general-purpose alignment algorithm. As this is local alignment, any subsequence of the query and target sequence may appear in the alignment.
    affine:overlap	This type of alignment finds the best overlap between the query and target. The overlap alignment must include the start of the query or target and the end of the query or the target sequence, to align sequences which overlap at the ends, or in the mid-section of a longer sequence.. This is the type of alignment frequently used in assembly algorithms.
    est2genome	This model is similar to the affine:local model, but it also includes intron modelling on the target sequence to allow alignment of spliced to unspliced coding sequences for both forward and reversed genes. This is similar to the alignment models used in programs such as EST_GENOME and sim4.
    ner	NERs are non-equivalenced regions - large regions in both the query and target which are not aligned. This model can be used for protein alignments where strongly conserved helix regions will be aligned, but weakly conserved loop regions are not. Similarly, this model could be used to look for co-linearly conserved regions in comparison of genomic sequences.
    protein2dna	This model compares a protein sequence to a DNA sequence, incorporating all the appropriate gaps and frameshifts.
    protein2dna:bestfit	NEW: This is a bestfit version of the protein2dna model, with which the entire protein is included in the alignment. It is currently only available when using exhaustive alignment.
    protein2genome	This model allows alignment of a protein sequence to genomic DNA. This is similar to the protein2dna model, with the addition of modelling of introns and intron phases. This model is simliar to those used by genewise.
    protein2genome:bestfit	NEW: This is a bestfit version of the protein2genome model, with which the entire protein is included in the alignment. It is currently only available when using exhaustive alignment.
    coding2coding	This model is similar to the ungapped:trans model, except that gaps and frameshifts are allowed. It is similar to a gapped tblastx search.
    coding2genome	This is similar to the est2genome model, except that the query sequence is translated during comparison, allowing a more sensitive comparison.
    cdna2genome	 This combines properties of the est2genome and coding2genome models, to allow modeling of an whole cDNA where a central coding region can be flanked by non-coding UTRs. When the CDS start and end is known it may be specified using the --annotation option (see below) to permit only the correct coding region to appear in the alignemnt.
    genome2genome	This model is similar to the coding2coding model, except introns are modelled on both sequences. (not working well yet) 


-n NumCpuToUse [20]
-o Outputfile[$fasta.exonerate.out]
-of outputformat []
    showalignment 
    showsugar
    showcigar
    showvulgar    
    showquerygff
    showtargetgff[1]

-f directoryToSaveRunIntermediateResults[exonerate.tmp]
-a ExtraArgumentforBlastn. Dont put hyphens.eg. -a 'evalue 1e-10' -a 'culling_limit 1' etc..
-r do not remove temp folder
";
die "$usage" if $help;

## join extra arguments for blast
unshift(@args," ") if @args>0;  ### add empty element as first in array so that all the arguments get '-' infront of them
$args=join(" -",@args);

## check if necessary parameters are given.
open(FASTA,"$fasta") or die "Unable to find sequence file $fasta\n$!\n$usage";
die "Provide database file to blast against\n$usage\n" if !$database;

## check if database is formatted. format if it is not. 
## create temp folder for saving divided sequences
mkdir $dir;
my$seqfilename=basename($fasta);
my$dbfilename=basename($database);
$out=$out?$out:"$seqfilename.vs.$dbfilename.exonerate.out";


my($seqhash)=ReadFasta($fasta);

my @seq_count=(keys %$seqhash);
my $seq_count=@seq_count;



my@outfile;
my$seq2write=my$start=0;
for my$i(0..$numcpu) {
  open(TMP,">$dir/$seqfilename.$i") or die "Could not create file. $!";
  push(@outfile,"$dir/$seqfilename.$i");
  $seq2write+= int($seq_count/$numcpu) == $seq_count/$numcpu?int($seq_count/$numcpu):int($seq_count/$numcpu)+1;
  my$lim=$seq_count>$seq2write?$seq2write:$seq_count;
  if ($i == $numcpu) {
    $lim=$seq_count;
  }

  for(my$j=$start;$j < $lim;$j++){
    print TMP ">$seq_count[$j]\n$$seqhash{$seq_count[$j]}\n"
  }
  $start=$seq2write;
}

## empty $seqhash to
undef %$seqhash;

my@cmds;
my@PIDs;
for (0..$numcpu) {
  $cmds[$_]="nice exonerate --model $prg  --maxintron 10000 --score 1000 --softmasktarget 1    --showalignment 1  --showsugar 0  --showcigar 0 --verbose 0  --showquerygff 0   --showtargetgff 1  $dir/$seqfilename.$_   $database 1> $dir/$seqfilename.$_.out  2> $dir/exonerate.$_.log";
}

my $Tempfile = 'job_controller.'.$$;

for (0..$numcpu) {
if ( $PIDs[$_] = fork() ) {

} elsif ( defined $PIDs[$_] ) {
    print "Started blastn run:$cmds[$_]\n";
    exec( "$cmds[$_]" );

} else {
    warn "Something broke: can't fork()!\n";
}
}
my$count=1;
print "Waiting for blastn runs to finish\n";
sleep(1);
print $count++." of ",$numcpu+1," blast runs are Done \n" while wait != -1 ;  # avoid zombies
print "\nAll blast runs are done\nCombining outputs of all the blast runs.\n";
sleep(5);



system("rm $out");
open(OUT,">>$out") or die "Cannot open outfile $out";
for my$filenum(0..$numcpu) {
 open my $TmpFH,"$dir/$seqfilename.$filenum.out" or die "Can't open \"$dir/$seqfilename.$filenum.out\": $!";
  #print "\nProcessing file num $filenum\n";
 if ($outfmt == 5) {
	while (<$TmpFH>){
	  if ($filenum == 0) {
	   next if $_ =~ /\<\/BlastOutput_iterations\>|\<\/BlastOutput\>/;
	   print OUT $_;
	   #print "\nNor Processing file num $filenum in eq 0\n";
	  }

	  elsif ($filenum > 0 && $filenum < $numcpu) {
	   next if $_ =~ /\<[\/]*BlastOutput|Parameters|\<\?xml|\<\!DOCTYPE/;
	   print OUT $_ ;
	   #print "\nNor Processing file num $filenum  in gt 0\n";
	  }
	  elsif($filenum == $numcpu){
	   next if $_ =~ /\<BlastOutput|<\/BlastOutput_param|<[\/]*Parameters|\<\?xml|\<\!DOCTYPE/;
	   print OUT $_ ;
	   #print "\nNor Processing file num $filenum  in eq 4\n";

	  }
	  else{
		print "\nNor Processing file num $filenum else loop\n";
	  }
	}
 }
 else{
  print OUT "\n".join("",<$TmpFH>) ;

 }

 close $TmpFH;
 unlink "$dir/$seqfilename.$filenum.out" if !$rem;
 unlink "$dir/$seqfilename.$filenum" if !$rem;
 print "Finished processing $dir/$seqfilename.$filenum\n";
}


print "\nremoving temp folder\n" if !$rem;
sleep(5);
system("rm -r $dir") if!$rem;
print "Exonerate run is finished\nResults are saved in file:\n$out\n\n";


















### subroutines
sub ReadFasta{ # to read fasta format files into hash. returns hash.

	my $seqfile=shift;


	my ($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "reading Sequences from input file.....Plz wait...\n";
	my%seq_hash=();
	#$seq_hash{'RS_Concatenated'}="";

	$/="\n>";    # Change record seperator to read Fasta
	my$last_N=1;
	while(<FASTA>){
    	chomp;
    	($header,@sequence)=split("\n",$_);

    	$header=~s/>//;						# Remove Leading > from Header
    	$header=~s/\s*$//;					# Remove trailing spaces from header
    	$header=~s/^\s*//;					# Remove Leading spaces from Header

    	my$sequence= join("",@sequence);
    	$sequence=~ s/\s//g;
    	$sequence=~s/\n//g;


    	if($header=~/^\s*$/){next;}

    	# Make sure no duplicate header exists else they will be overwritten and lost during hashing due to same key.
    	if(!exists $seq_hash{$header}){
    		$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
			#$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
		}
		else {

			# find a uniq header name by adding a number at the end. If header still exists, increase the number by one
			while(exists $seq_hash{$header}){$header=$header.$last_N;$last_N++;}

			$seq_hash{$header}=$sequence;
			#$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;


		}
	}

	my @seq_count=keys (%seq_hash);
	my $seq_count=@seq_count;

	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.

	return(\%seq_hash);

}

###############################################
