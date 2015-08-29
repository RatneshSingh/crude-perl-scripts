#!/usr/bin/perl
use strict;
use warnings;
use Bio::SearchIO;
use Getopt::Long;


# Usage to be printed out for user.
my$usage="

This script reads standard blast result and convert it into blast table format.
The blast results can be filtered using several parameter. 

Usage: Perl Script -i blast_result_file [options]

-l	Alignment length to filter. Results above this threshold will be reported [All]
-p	Percent value to filter. Results above this threshold will be reported [50].
-e	Evalue value to filter. Results below this threshold will be reported [All].
-c	Query Coverage value to filter. Results above this threshold will be reported.
	Coverage is calculated for each hsp individually [50]. 
-d	Description in hit. Hits matching description will be reported [All].
-o	Output file to save result [STDOUT]
-s	Sequence file used as database for blast.
-a	Report coverage in increment of[5] 
-b	Report percent in increment of[5]
-m	Minimum value for coverage to start counting from[60]
-x	Maximum value for coverage to stop counting on[100]
-n	Minimum value for percent to start counting from[60]
-y	Maximum value for percent to stop counting on[100]

";




our(%sequence,%sequencelength,%seqnamelist,%matrix,%matrixcalc,@blasttable);
our ($format,$infile,$outfile,$evalue,$coverage,$percent,$description,$alnlength,$OUT,$OUT2,$sequence,$binC,$binP,$minC,$maxC,$minP,$maxP) = ('blast');

GetOptions(
	   'f|format:s'   => \$format,
	   'i|infile:s'    => \$infile,
	   'o|outfile:s'   => \$outfile,
	   'e|evalue:f'	  => \$evalue,
	   'c|coverage:f' => \$coverage,
	   'p|percent:f'  => \$percent,
	   'd|description:s' => \$description,
	   'l|alnlength:i' => \$alnlength,
	   's|sequence:s' => \$sequence,
	   'a|inccoverage:i'=>\$binC,
	   'b|incpercent:i'=>\$binP,
	   'm|minC:i'=>\$minC,
	   'x|maxC:i'=>\$maxC,
	   'n|minP:i'=>\$minP,
	   'y|maxP:i'=>\$maxP 
	   );
	   

# Pick the file name from commandline if option flags are are not used.
if(!$infile){die "\nFile contining blast result in required for the script to run\n\n$usage";} 



my$sequencelength=0;
my$sequencenumber=0;

if($sequence){
	open SEQ,"$sequence" or die "Cannot open sequence file: $sequence\n$usage\n";
	while(<SEQ>){
	
		if(/>/){
			$sequencenumber++;
			next;
		} 
	
		else{
			my$N=($_=~tr/ATGCatgc/ATGCatgc/); 
			$sequencelength+=$N;
		}
	}
	
	
	close SEQ;
}
else {die "\nSequence file used as query for blast is required for the script to run\n\n$usage";}

print "Total number of sequences in sequence file:$sequencenumber\n";
print "Total number of sequences in sequence file:$sequencelength\n";





# print Usage instructions when blast file name is not provided.
print "$usage" if !defined $infile;


#Set parameters to default values for filtering when not defined.
if(!$coverage){$coverage=50;}
if(!$percent){$percent=50;}
if(!$evalue){$evalue=10;}
if(!$description){$description="";}
if(!$alnlength){$alnlength=1;}	   
if(!$binC){$binC=5;}
if(!$binP){$binP=5;}
if(!$minC){$minC=60;}
if(!$maxC){$maxC=100;}
if(!$minP){$minP=60;}
if(!$maxP){$maxP=100;}

print "Reading blast file\n";  
my $searchIO = Bio::SearchIO->new(-format => $format, -file   => $infile);


# Process blast results.
while( my $result = $searchIO->next_result ) {
    while( my $hit = $result->next_hit ) {
		while( my $hsp = $hit->next_hsp ) {
	    
	    	# calculate mismatches. Not provided by SearchIO methods.
		    #my $mismatch = $hsp->length('total') - $hsp->num_conserved - $hsp->gaps;
		   	
		   	# Calculate query coverage using Alignment length. it includes gaps.
		   	#my$current_coverage=($hsp->length('total')/$result->query_length)*100; 
	
		   	# Calculate query coverage using Q start and Q end. Does not includes gaps.
		   	my$current_coverage= ($hsp->length('query')/$result->query_length)*100; 	   	
		   	#print $hsp->query->strand < 0 ? ( $hsp->query->start - $hsp->query->end ):( $hsp->query->end - $hsp->query->start )."\n";
		   	
		   	
		   	#print join ("\t",$hsp->length('total'),$result->query_length, $current_coverage, "\n" );
		   	
		   	   if(	$hsp->hsp_length >= $alnlength && 
		   	   		$hsp->percent_identity >= $percent && 
		   	   		$hsp->evalue <= $evalue && 
		   	   		$current_coverage >= $coverage && 
		   	   		$hit->description =~ m/$description/i
		   	   	)
		   	   		{
		    				#print $hit->description."\n";
		    				push(@blasttable, join("\t", ($result->query_name,
									    sprintf("%.2f",$hsp->percent_identity),
									    $current_coverage,
									    $result->query_length)));
							
				
							$sequence{$result->query_name}{'length'}=$result->query_length;
												
#							if(!$sequence{$result->query_name}{'coverage'}){ 
#								$sequence{$result->query_name}{'coverage'}=$current_coverage;
#							}
#							
#							elsif($sequence{$result->query_name}{'coverage'}<$current_coverage){
#								$sequence{$result->query_name}{'coverage'}=$current_coverage;	
#							}
#							
#							
#							
#							if(!$sequence{$result->query_name}{'percent'}){ 
#								$sequence{$result->query_name}{'percent'}=$hsp->percent_identity;
#							}
#							
#							elsif($sequence{$result->query_name}{'percent'}<$hsp->percent_identity){
#								$sequence{$result->query_name}{'percent'}=$hsp->percent_identity;	
#							}
							

							
							
							
							
				
				}	
				
				}



	}
}





# store the names of sequences in seperate hash for every sombination of coverage and percent identity

for(my$c=$minC;$c<=$maxC;$c=$c+$binC){
	for(my$p=$minP;$p<=$maxP;$p=$p+$binP){
		foreach my$blastline(@blasttable){
			my@blastelements=split(/\t/,$blastline);
				
			if($blastelements[1] >= $p && $blastelements[2] >= $c){
				$matrix{$p}{$c}{$blastelements[0]}=();
 
			}
		}
	}		
}			

# calculate the number percent 			
for(my$c=$minC;$c<=$maxC;$c=$c+$binC){
	for(my$p=$minP;$p<=$maxP;$p=$p+$binP){
		
		$matrixcalc{$p}{$c}{'length'}=0;
		$matrixcalc{$p}{$c}{'reads'}=0;
		
		foreach my$seqname(keys %{$matrix{$p}{$c}}){
			$matrixcalc{$p}{$c}{'reads'}=$matrixcalc{$p}{$c}{'reads'}+1; 
			$matrixcalc{$p}{$c}{'length'}=$matrixcalc{$p}{$c}{'length'}+$sequence{$seqname}{'length'}; 
				
		}		
		# calculate the percent contribution of read number and length to total
		$matrixcalc{$p}{$c}{'perlength'}=($matrixcalc{$p}{$c}{'length'}*100)/$sequencelength;
		$matrixcalc{$p}{$c}{'perreads'}=($matrixcalc{$p}{$c}{'reads'}*100)/$sequencenumber;

		
	}
	
}	

	
		
		#print "\nCalculated read number: $matrix{$p}{$c}{'reads'}\n";
		






if( $outfile ) {
	open($OUT,'>',"length.$outfile") || die "Unable to open $outfile for writing\n $usage";
	open($OUT2,'>',"number.$outfile") || die "Unable to open $outfile for writing\n $usage";
}
else {
	$OUT = \*STDOUT; 
	$OUT2 = \*STDOUT;
}




print $OUT "Percent contribution of filtered reads to total reads\n";
# print matrix in %

print $OUT 'c\p';

for(my$p=$minP;$p<=$maxP;$p=$p+$binP){
	print $OUT "\t$p";
	#print $OUT2 "\t$p";
}


for(my$c=$minC;$c<=$maxC;$c=$c+$binC){
	
	print $OUT "\n$c";
	#print $OUT2 "\n$c";


	for(my$p=$minP;$p<=$maxP;$p=$p+$binP){
		#my$perlength=($matrix{$p}{$c}{'length'}*100)/$sequencelength;
		#my$pernumber=($matrix{$p}{$c}{'reads'}*100)/$sequencenumber;

		print $OUT "\t",sprintf("%.5f",$matrixcalc{$p}{$c}{'perreads'});
		#print $OUT2 "\t",sprintf("%.5f",$pernumber);
		
		#print "length:$matrix{$p}{$c}{'length'}\nnumber:$matrix{$p}{$c}{'reads'}\n";
		#print "Perlength:$perlength\tPernumber:$pernumber\n";
		
	}
}

print "\n\n\n";



# print matrix in %
print $OUT2 "Percent contribution of sum of the filtered read length to the sum of total read length\n";
print $OUT2 'c\p';

for(my$p=$minP;$p<=$maxP;$p=$p+$binP){
	#print $OUT "\t$p";
	print $OUT2 "\t$p";
}



for(my$c=$minC;$c<=$maxC;$c=$c+$binC){
	
	#print $OUT "\n$c";
	print $OUT2 "\n$c";


	for(my$p=$minP;$p<=$maxP;$p=$p+$binP){
		#my$perlength=($matrix{$p}{$c}{'length'}*100)/$sequencelength;
		#my$pernumber=($matrix{$p}{$c}{'reads'}*100)/$sequencenumber;

		#print $OUT "\t",sprintf("%.5f",$perlength);
		print $OUT2 "\t",sprintf("%.5f",$matrixcalc{$p}{$c}{'perlength'});
		
		#print "length:$matrix{$p}{$c}{'length'}\nnumber:$matrix{$p}{$c}{'reads'}\n";
		#print "Perlength:$perlength\tPernumber:$pernumber\n";
		
	}
}

print "\n";