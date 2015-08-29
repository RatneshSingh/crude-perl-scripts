#!usr/bin/perl -w
use strict;
use Getopt::Std;


our($opt_d,$opt_b,$opt_l,$opt_o,$opt_q);
getopt('dbloq');


my($fastafile,$blast_table,$qlength,$sstart1,$send1,$extend,$output, %seq, $header, $sequence,@sequence,$extra,$key,$revcom,%lengths);
print "\n**********************************************************************\n\n
This program will read blast (blastn) tabular output file and sequence(subject) file. 
it will chop subject on the basis of given alignment cordinates.
chopped subject sequence will be of the length of query assigned by the user.
it will reverse compliment the sequences which are in opposite direction
use it interactively or on the command line.
\nfor command line use => 'program_name' -d sequence_file -b blast_output_file  -l extra_nt -o output_file -q list_of_query_length

\n***********************************************************************\n\n  ";

#*********************************************************************************************************
# open fasta file containing subject sequences used in the blast as database..
if($opt_d){	$fastafile=$opt_d;}

else{
	print "Type the name of file containing sequences in fasta format\n";
	$fastafile=<STDIN>;
}
if ($fastafile){
	chomp($fastafile);
	open (FASTAFILE,$fastafile) or die"Cannot open the file $fastafile";
}
#************************************************************************************************************

# open blast output file in tabular format. here it is tblastn where one query is blasted against a set of subject sequences. query is  aminoacid file and subject nucleotide.
if($opt_b){
	$blast_table=$opt_b;
}
else{
	print "Type the name of file containing tabular blast output\n";
	$blast_table=<STDIN>;
}
chomp($blast_table);

open (BLASTFILE,$blast_table) or die"Cannot open the file $blast_table";

#***************************************************************************************************************
# open file containing sequence length information or ask for query length manually.

if($opt_q){open QUERYLENGTH,$opt_q;}
else{
	print "Please provide the length of the query sequence or provide the list with -q option at the command line for multiple query lengths\n Type 0 if you do not want to use query length extension.\n";
	$qlength=<STDIN>;
	chomp($qlength);
}



#*****************************************************************************************************************


#ask for the extra length needed to be cut on both sides of sequences
if($opt_l){	$extra=$opt_l;}
else{
	print "Type the extra length you want to added in sequence on both sides\n";
	$extra=<STDIN>;
}
chomp($extra);



#*****************************************************************************************************************

# Ask for output file name where choped sequences will be saved.
if($opt_o){$output=$opt_o;}
else{
	print "\n Type output file name\n";
	$output=<STDIN>;
}

chomp($output);

open (NEW_FASTA,">$output");
open (TSD,">$output.extra.txt");
#******************************************************************************************************************
#read query length file to get the query's length

if($opt_q){
	print "Reading query lengths from file $opt_q\n";
	$/="\n"; # define record seperator to read the file.
	while(<QUERYLENGTH>){#1

		my($query,$qlength)=split(":n:",$_);
		$query=~s/>//;
		$query=~s/^\s+//;
		my@query=split(/\s+/,$query);
		$query=$query[0];
		$qlength=~ s/\D//;
		$lengths{$query}=$qlength;
		#print "query:$query\t length :$qlength\t length_hash:$lengths{$query}\n";
		#print "\n\n";
		
	}#1
}



#******************************************************************************************************************
#read blast database sequences from file in a hash table named  %seq
print "Reading database sequences in memory...Plz wait......";
$/="\n>"; # define record seperator for fasta sequence file.
while(<FASTAFILE>){#1

	($header,@sequence)=split("\n",$_);
	$header=~s/>//;
	$sequence=join("",@sequence);
	$sequence=~ s/\s//;
	$sequence=~ s/>//;
	$seq{$header}=$sequence;
	$sublength{$header}=length($sequence);
	#print keys %seq;
	#print "\n\n";
	
}#1


#*****************************************************************************************************************
# parse blast tabular output to find out the cordinates to chop the sequence

$/="\n";   # define record seperator for blast file. current record seperator is "\n>" as defined for sequence file.
my $count=0;
print".....Done.\nProcessing blast file.....Plz wait.... ";
while(<BLASTFILE>){#1
	chomp;
	next if($_ =~ /^\s*$/);

	my ($query,$subject,$identity,$length,$mismatch,$gap,$qstart,$qend,$sstart,$send,$evalue,$bitscore)=split(/\t/,$_);
	chomp ($query,$subject,$identity,$length,$mismatch,$gap,$qstart,$qend,$sstart,$send,$evalue,$bitscore);
	#print "$subject,$identity,$length,$mismatch,$gap,$qstart,$qend,$sstart,$send,$evalue,$bitscore\n\n";
	my $extend=0;
	#calculate starting and end point for chopping subject sequence.

	if($opt_q){


		if ($sstart<$send){#2
			$sstart1=$sstart-$extra-$qstart-1;
			# check if query value becomes less than 0
			if($sstart1<0){$sstart1=0;}
			
			my$send_extend=$lengths{$query}-$qend;
			#print" send:$send_extend\tqlength:$lengths{$query}\tqend:$qend\n";
			
			$send1=$send+$extra+$send_extend;
			# test if subject end extends beyond its actual length
			if($send1>$sublength{$subject}){$send1=$sublength{$subject};}
			
			$extend=$send1-$sstart1;
			$revcom=0;
		}#2

		elsif($sstart>$send) {#3
			my$send_extend=$lengths{$query}-$qend;
			
			my	$sstart_N=$send-$extra-$send_extend-1 ;
			if($sstart_N<0){$sstart_N=0;}
			
			my	$send_N=$sstart+$extra+$qstart;
			if($send_N>$sublength{$subject}){$send_N=$sublength{$subject};}
			
			$extend=$send_N-$sstart_N;
			$sstart1=$sstart_N;
			$send1=$send_N;
			$revcom=1;
		}#3
	}
	else{
	
		if ($sstart<$send){#2
			
			$sstart1=$sstart-$extra-$qstart-1;
			# check if query value becomes less than 0
			if($sstart1<0){$sstart1=0;}

			my$send_extend=0;
			
			if($qlength==0){$send_extend=0}
			else{$send_extend=$qlength-$qend}
			
			$send1=$send+$extra+$send_extend;
			# test if subject end extends beyond its actual length
			if($send1>$sublength{$subject}){$send1=$sublength{$subject};}
			
			$extend=$send1-$sstart1;
			$revcom=0;
		}#2

		elsif($sstart>$send) {#3
			my$send_extend=$lengths{$query}-$qend;
			
			my	$sstart_N=$send-$extra-$send_extend-1 ;
			if($sstart_N<0){$sstart_N=0;}
			
			my	$send_N=$sstart+$extra+$qstart;
			if($send_N>$sublength{$subject}){$send_N=$sublength{$subject};}
			
			$extend=$send_N-$sstart_N;
			$sstart1=$sstart_N;
			$send1=$send_N;
			$revcom=1;
		}#3
	
	}

	#*****************************************************************************************************************
	# chop sequence and write in to the out put file.
	#print "\n\n\n$subject---subject befor entering subroutine ";
	my ($header,$choped_sequence)= search_seq($subject,$sstart1,$extend,$revcom);
	#print "\n$subject---subject after calling subroutine is \n\n\n";

	if($header){#6
		print NEW_FASTA ">$query:$header:$identity:$length:$mismatch:$gap:$qstart:$qend:$sstart:$send:$evalue:$bitscore:sstart-$sstart1:send-$send1:extend-$extend \n$choped_sequence\n";

		my $head= substr($choped_sequence,0,20);

		my $seq_len= length($choped_sequence);
		my $tail_start= $seq_len-20;
		my $tail= substr($choped_sequence,$tail_start,20);
		print TSD "\n$head\t$tail\t$header";
		$count++;
		#print "processing line : $count\n"
				
	}#6
	
	else {#7
		print NEW_FASTA ">$subject could not be found \n";

	}#7



}#1




print ".....Done.\n$count lines processed";
close FASTAFILE;
close BLASTFILE;
close NEW_FASTA;

exit;
#****************************************************************************************************************
# subroutine to search header matching subject name and chop on the basis of given cordinates. return choped substring.

sub search_seq {#1

	my ($subject,$sstart,$extend,$revcom)=@_;
	$subject=~s/^\s+//;
	my @subject=split(/\s+/,$subject);

	#my @subject=split(/\|/,$subject);

	my $accession=$subject[0];

	chomp($accession);

	foreach $key(keys %seq){#2
		$key=~s/^\s+//;
		my@key_N=split(/\s+/,$key);
		$key=$key_N[0];
		
		if ($key=~/$accession/){#3
			my $new_seq=substr($seq {$key},$sstart,$extend);
			if($revcom==1){
				my $new_seq_rev= reverse $new_seq;
				$new_seq_rev=~ tr/ATGCatgc/TACGtacg/;
				$new_seq=$new_seq_rev;
			}
		return ($key,$new_seq);
		}#3

	}#2

}#1


