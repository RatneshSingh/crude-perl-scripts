#!/usr/bin/perl
use strict;
use warnings;

print "\n\n\n*****USAGE*****\nScript\tspidey_output_file\n\n";

my$spidey_output_file=$ARGV[0] if defined $ARGV[0];
#my$gene_output=$ARGV[2] if defined $ARGV[2];
#my$summary_output=$ARGV[3] if defined $ARGV[3];
open SUMMARY,">$spidey_output_file.summary" or die "spidey_output_file not defined\n";

#########################################################################################################
# Parsing the spidey output file
print "reading Spidey out put file.........\n";
spidey_parse($spidey_output_file);


###########################################################################################
# Subroutine for parsing spidey output file for the gene name, co-ordinates of exon and 
# introns.Output table with all the co-ordinates arranged in a nice format.

sub spidey_parse{
	open SPIDYOUT,(my$spidey_output= shift(@_)) or die "Cannot find Spidey output file\n";
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
				
				if($mRNA_info=~/^Genomic:[\s\w]*\|([\w\d\|_-]*)\s*([\w\d\s]*)\,\s*([\d]*)\s*bp/){
				$genomic=$1;
				#print SUMMARY"Genomic DNA: $genomic\n";
				print SUMMARY"$genomic\t";
				print "$genomic\t";


				}
 		
 				if($mRNA_info=~/^mRNA:[\s\w]*\|([\w\W]*?)\s+([\s\w]*)\s*,\s*([\d]*)\s*bp/){
 				$mRNA=$1;
 				#print "mRNA: $mRNA\n";
 				print SUMMARY"$mRNA\n";
 				print "$mRNA\n";


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
						print SUMMARY"\tExon $i\t$exon_start{$i}\-$exon_end{$i}\n";
						print "\tExon $i\t$exon_start{$i}\-$exon_end{$i}\n"	;
        				}# if loop.				
				
					}# For loop.

				}#if loop
			
			
			next;
			}# foreach loop.
	
	
	my@a= sort values(%exon_start);
	my@b= sort values(%exon_end);
	my@c= sort(@a,@b);
	my$l=@c;
	
	
	}# while loop.

print "Done....\n";
close(SUMMARY);

# information collected: Genomic DNA=>$genomic,mRNA Name=>$mRNA,Strand=>$strand,Number of exons=>$exon_N,
# splice sites=>$splice_site,mRNA coverage=>$mRNA_coverage,Missing mRNA ends: $missing_mRNA_ends
# Exon start=>$exon_start{$i},Exon_end=>$exon_end{$i},mRNA start=>$mRNA_start{$i},
# mRNA end=>$mRNA_end{Si},$exon_percent_identity{$i},$exon_gap{$i}




}# End of subroutine spidey_parse






