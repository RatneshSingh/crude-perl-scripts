#!/usr/bin/perl
use strict;
use warnings;

print "*****USAGE*****\n\nScript\tseqfile\tspidey_output_file\tgene_outputfile\tsummary_output\n\n";

my (%seq_hash);

my$seqfile=$ARGV[0];
my$spidey_output_file=$ARGV[1];
#my$gene_output="$seqfile.gene.out";
#my$summary_output="$seqfile.summary.out";
my$gene_output=$ARGV[2];
my$summary_output=$ARGV[3];

open GENE,">$gene_output" or die "gene_output not defined\n";
open FASTA,"$seqfile" or die "Seqfile not defined\n";
open SUMMARY,">$summary_output" or die "summary_output not defined\n";

####################################################################################################
# Read database sequence in to a hash %seq_hash{}, remove ">" from header. 
# remove any white spaces and newline from sequence.

print "reading Sequences from input file.....Plz wait...\n";
$/="\n>";

while(<FASTA>){#
    chomp;
    my($header,@sequence)=split("\n",$_);
    $header=~s/>//;
    $header=~s/\s//g;

    my$sequence= join("",@sequence);
    $sequence=~ s/\s//g;
    $sequence=~s/\n//g;
    #print "$header\n";
    #feed headers and sequences in hash.
    $seq_hash{$header}=$sequence;

    #print "$seq_hash{$header}\n\n";
}#
my @seq_count=keys (%seq_hash);
my $seq_count=@seq_count;

print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
@seq_count=();
############################################################################################

#########################################################################################################
# Parsing the spidey output file
print "reading Spidey out put file.........\n";
spidey_parse($spidey_output_file);











###########################################################################################
# Subroutine for parsing spidey output file for the gene name, co-ordinates of exon and 
# introns.Output table with all the co-ordinates arranged in a nice format.

sub spidey_parse{
open SPIDYOUT,(my$spidey_output= shift(@_));
$/= "\n--SPIDEY";
#print "Spidey output file is:$spidey_output\n";

	while(<SPIDYOUT>){
		my ($genomic,$mRNA,$strand,$missing_mRNA_ends,%gene,@mRNA_info)=();
		my ($exon_N,$splice_site,$mRNA_coverage,%exon_start,%exon_end,%mRNA_start,%mRNA_end,%exon_percent_identity,%exon_gap,$percent_identity,%exon_mismatch)=0;
		my$chunk=$_;
		$chunk=~s/(--SPIDEY\s*version\s[\d\W]*?--)//g;
		$chunk=~s/(--SPIDEY)//g;

#print "\n\nChunk read is : \n$chunk\n";

		if($chunk=~/(Genomic:[\s\w\d\W\S\D]*?Missing mRNA ends:\s*neither)/){
		$chunk=$1;
		#print"\nprinting Dollar 1:$1\n";
		}  
#print "\n\nSelected chunk is: \n$1\nxxxxxx\n";

		@mRNA_info= split(/\n/,$chunk);  
#print "First line of the chunk is : $mRNA_info[0]\n";			
			
			foreach my$mRNA_info(@mRNA_info){
				
				chomp($mRNA_info);
				
				if($mRNA_info=~/^Genomic:[\s\w]*\|([\.\w\W\|\_\-]*?)\s+([\w\s]*)\,\s*([\d]*)\s*bp/){
				$genomic=$1;
				$genomic=~s/\s//g;

				#print SUMMARY"Genomic DNA: $genomic\n";
				print SUMMARY"$genomic\t";

				}
 		
 				if($mRNA_info=~/^mRNA:[\s\w]*\|([\w\W]*?)\s+([\s\w]*)\s*,\s*([\d]*)\s*bp/){
 				$mRNA=$1;
 				#print "mRNA: $mRNA\n";
 				print SUMMARY"$mRNA\n";

 				}
 		
 				if($mRNA_info=~/^Strand:\s*([minuspl]*)\s*$/){
 				$strand=$1;
 				#print "Strand: $strand\n";
 				}
 		
 				if($mRNA_info=~/^Number of exons:\s*([\d]*)\s*$/){
 				$exon_N=$1;
 				#print "Number of exons: $exon_N\n";
 				}

				if($mRNA_info=~/^Number of splicing sites:\s*([\d]*)\s*$/){
 				$splice_site=$1;
 				#print "Number of splice sites: $splice_site\n";
 				}
				
				if($mRNA_info=~/^mRNA coverage:\s*([\d]*)%\s*$/){
				$mRNA_coverage=$1;
 				#print "mRNA coverage: $mRNA_coverage\n";
 				}

				if($mRNA_info=~/^overall percent identity:\s*([\d\W]*)%\s*$/){
				$percent_identity=$1;
 				#print "percent_identity: $percent_identity\n";
 				}

				if($mRNA_info=~/^Missing mRNA ends:\s*([a-zA-Z]*)\s*$/){
				$missing_mRNA_ends=$1;
 				#print "Missing mRNA ends: $missing_mRNA_ends\n";
 				}

				#else{print"Could not parse line\n $mRNA_info\n Dont have any thing to match\n";}
				if($mRNA_info=~/^Exon/){
				for(my$i=1;$i<$exon_N+1;$i++){
						
						if($mRNA_info=~/Exon\s*$i+([\s\w\d\W]*?)splice\s*site/){
       						#my$exon_in=$1;
        					my@exon_info=split(/\s+/,$1);
        					#print "@exon_info\n";
        					($exon_start{$i},$exon_end{$i})=split(/-/,$exon_info[1])if defined $exon_info[1];
       						($mRNA_start{$i},$mRNA_end{$i})=split(/-/,$exon_info[3])if defined $exon_info[3];
       						$exon_percent_identity{$i}=$exon_info[6]if defined $exon_info[6];
       						$exon_mismatch{$i}=$exon_info[8]if defined $exon_info[8];
 							$exon_gap{$i}=$exon_info[10]if defined $exon_info[10];
						#print "Exon $i Start:$exon_start{$i}\tExon $i End:$exon_end{$i}\n"
						print SUMMARY"\tExon $i\t$exon_start{$i}\-$exon_end{$i}\n"	
        				}# if loop.				
				
				}# For loop.

			}
			
			
			next;
			}# foreach loop.
	
	
	my@a= sort values(%exon_start);
	my@b= sort values(%exon_end);
	my@c= sort(@a,@b);
	my$l=@c;
	#print "\n$a[$l-1]\n";
	#print "@a";
	
	
	
	#print "\ngene_start=$exon_start{1}\tgene_end=$exon_end{$exon_N}";
	
	($gene{$mRNA})= extract_gene($c[0],$c[$l-1],$genomic,$strand);
	print GENE">Gene_$mRNA\n$gene{$mRNA}\n" if defined $gene{$mRNA} ;
	#print"\n\n";
	}# while loop.

print "Done....\n";
close(GENE);
close(SUMMARY);
close(FASTA);
#close(GENE);

# information collected: Genomic DNA=>$genomic,mRNA Name=>$mRNA,Strand=>$strand,Number of exons=>$exon_N,
# splice sites=>$splice_site,mRNA coverage=>$mRNA_coverage,Missing mRNA ends: $missing_mRNA_ends
# Exon start=>$exon_start{$i},Exon_end=>$exon_end{$i},mRNA start=>$mRNA_start{$i},
# mRNA end=>$mRNA_end{Si},$exon_percent_identity{$i},$exon_gap{$i}




}# End of subroutine spidey_parse

###########################################################################################
# Subroutine for the extraction of gene from the genomic DNA based on co-ordinates obtained
# from parsing spidey summary file. Save gene sequences in a file. Will try to align cDNA 
# with genomic DNA.

sub extract_gene{

	my($exon_start,$exon_end,$genomic,$strand)=@_;
	chomp($exon_start,$exon_end,$genomic,$strand);
	print "Exon:$exon_start,$exon_end,\tGenomic:$genomic,\tStrand:$strand\n\n";

	#print "infor received in subroutine extract_gene:\n$exon_start,$exon_end,$genomic,$strand\n";
	my($gene,$length)=();
	if($exon_start<$exon_end){
		$length=$exon_end-$exon_start+1;
		$exon_start=$exon_start-1;
		#print "\nFetching substring with co-ordinates:Start:$exon_start\tEnd:\t$exon_end\tLength:\t$length\n"; 

		$gene=substr($seq_hash{$genomic},$exon_start,$length);
	}

	elsif($exon_start>$exon_end){
		$length=$exon_start-$exon_end+1;
		$exon_start=$exon_end-1;
		#my $exon_end_N=$exon_start;
		#print "\nFetching substring with co-ordinates:Start:$exon_start_N\tEnd:\t$exon_start\tLength:\t$length\n"; 
		$gene=substr($seq_hash{$genomic},$exon_start,$length);
	}

	else{print"\nExon Start and end figures either are equal or not a numericl value\n";}

	if($strand eq 'minus'){
		$gene= reverse($gene);
		$gene=~tr/ATGCNatgcn/TACGNtacgn/;
		}
	return($gene);
	
	($exon_start,$exon_end,$length)=0;
	($gene,$genomic,$strand)=();

}
