use warnings;
use strict;
use Getopt::Long;
use File::Basename;
#use Sys::CpuAffinity;

our($fasta,$out,$args,$help,$database,$dir,@args,$rem,$outfmt);
$dir="blastn.tmp";
$outfmt=0;
$args="";
my$numcpu = 20;
my$prg='blastn';
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

my $usage="$0 -s sequencefile -db database  [options]
-p Program to run[blastn]
-n NumCpuToUse [20]
-o Outputfile[$fasta.blastn]
-of outputformat [0]
     0 = pairwise,
     1 = query-anchored showing identities,
     2 = query-anchored no identities,
     3 = flat query-anchored, show identities,
     4 = flat query-anchored, no identities,
     5 = XML Blast output,
     6 = tabular,
     7 = tabular with comment lines,
     8 = Text ASN.1,
     9 = Binary ASN.1,
    10 = Comma-separated values,
    11 = BLAST archive format (ASN.1)
    12 = JSON Seqalign output
-f directoryToSaveResults[blastn.tmp]
-a ExtraArgumentforBlastn. Dont put hyphens before flag.eg. -a 'evalue 1e-10' -a 'culling_limit 1' etc..
-r do not remove temp folder
";
die "$usage" if $help;

## join extra arguments for blast
unshift(@args," ") if @args>0;  ### add empty element as first in array so that all the arguments get '-' infront of them
$args=join(" -",@args);

## check if necessary parameters are given.
open(FASTA,"$fasta") or die "Unable to find sequence file $fasta\n$!\n$usage";
die "Provide database file to blast against\n$usage\n" if !$database;


## create temp folder for saving divided sequences
mkdir $dir;
my$seqfilename=basename($fasta);
my$dbfilename=basename($database);
$out=$out?$out:"$seqfilename.vs.$dbfilename.blastn";


my($seqhash)=ReadFasta($fasta);

my @seq_count=(keys %$seqhash);
my $seq_count=@seq_count;
my@outfile;
my$seq2write=my$start=0;
for my$i(0..$numcpu) {
  open(TMP,">$dir/$seqfilename.$i") or die "Could not create file. $!";
  push(@outfile,"$dir/$seqfilename.$i");
  $seq2write+= int($seq_count/$numcpu);
  my$lim=$seq_count>$seq2write?$seq2write:$seq_count;
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
  $cmds[$_]="$prg  -query $dir/$seqfilename.$_  -db $database   -out $dir/$seqfilename.$_.out -outfmt \"$outfmt\"   $args > $dir/blastn.$_.log 2>&1";
  #$cmds[$_]="ls" ### for testing script witjout running blast everytime.
  ## extra options sometimes useful.
  ##-culling_limit 1 -num_threads 30 -evalue 1e-10 -outfmt '6 std qcovhsp qcovs qlen slen sscinames scomnames stitle'
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
sleep(1);
print "Waiting for blastn runs to finish\n";
print $count++." of ",$numcpu+1," blast runs are Done \n" while wait != -1 ;  # avoid zombies
print "\nAll blast runs are done\nCombining outputs of all the blast runs.\n";
sleep(5);



system("rm $out");
open(OUT,">>$out") or die "Cannot open outfile $out";
for my$filenum(0..$numcpu) {
 open my $TmpFH,"$dir/$seqfilename.$filenum.out" or die "Can't open \"$dir/$seqfilename.$filenum.out\": $!";
  #print "\nProcessing file num $filenum\n";
  ## need special way of combining files if output is XML format. 
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
system("rm -r $dir") if !$rem;
print "Blast run finished\nResults are saved in file:\n$out\n\n";


















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
