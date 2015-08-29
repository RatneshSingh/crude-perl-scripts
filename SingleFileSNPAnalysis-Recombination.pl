#/usr/perl/bin -w
use strict;
use Time::localtime;
use Getopt::Std;

open LOG, ">>$0.Log.txt";
print LOG ctime(), "\n@ARGV\n";
print ctime(), "\n$0 @ARGV\n";

my ($SNP_equals_Ref2,
				$SNP_equals_Ref1_and_Ref2,
            	$SNP_equals_Ref2_and_other,
            	$SNP_equals_Ref1_and_Ref2_and_other,
            	$ref1,
            	$ref2,
				$other,
            	$non_processed_SNPs
);

open SNP1, "$ARGV[0]" or die "cannot find $ARGV[0]";

my $totalSNP1 = 0;
my%recomb;
my%mappingList;


# Read file 1 for processing;
print "Reading file: $ARGV[0]";
while (<SNP1>) {
	next if $_=~/^\s*$/;
    $totalSNP1++;
    my $blank                  = ();
    my $Mapping                = ();
    my $ReferencePosition      = ();
    my $ConsensusPosition      = ();
    my $VariationType          = ();
    my $Length                 = ();
    my $Reference              = ();
    my $Variants               = ();
    my $AlleleVariations1      = ();
    my $Frequencies            = ();
    my $Counts                 = ();
    my $Coverage               = ();
    my $Variant_1              = ();
    my $Frequencyof_1          = ();
    my $Countof_1              = ();
    my $Variant_2              = ();
    my $Frequencyof_2          = ();
    my $Countof_2              = ();
    my $OverlappingAnnotations = ();
    my $AminoAcidChange        = ();
    s/[\,]//g;
    s/[\"]+/\"/g;
    (
        $blank,                  $Mapping,       $ReferencePosition,
        $ConsensusPosition,      $VariationType, $Length,
        $Reference,              $Variants,      $AlleleVariations1,
        $Frequencies,            $Counts,        $Coverage,
        $Variant_1,              $Frequencyof_1, $Countof_1,
        $Variant_2,              $Frequencyof_2, $Countof_2,
        $OverlappingAnnotations, $AminoAcidChange
    ) = split( /\"/, $_ );

#print "\n--->Mapping: $Mapping,Refposition: $ReferencePosition,ConcensusPosition: $ConsensusPosition,Variation: $VariationType, Length: $Length, Reference: $Reference, Variants:$Variants, AlleleVariations:$AlleleVariations1, Frequencies:$Frequencies, Counts:$Counts, Coverage:$Coverage,Var1: $Variant_1, Freq1:$Frequencyof_1,Count1: $Countof_1, Var2:$Variant_2, Freq2:$Frequencyof_2, Count2:$Countof_2, OverlapAnnot:$OverlappingAnnotations, AaChange:$AminoAcidChange       \n";
    $Mapping           =~ s/\s//g;
    
	$mappingList{$Mapping}=();
	$recomb{$Mapping}{'len'}=1 if !defined $recomb{$Mapping}{'len'};
	$recomb{$Mapping}{'len'}=$ConsensusPosition if $recomb{$Mapping}{'len'} < $ConsensusPosition ;
	

        #clean and sort allele variation for comparison
        $AlleleVariations1=~ s/\W//;

        my $Allelevarlen1 = length($AlleleVariations1);

        #print "\nBefore: $AlleleVariations1\t$AlleleVariations1\n";

        #comparison of SNPs begin from here


            if ( $Allelevarlen1 == 1 ) {
				$SNP_equals_Ref2++;
				$ref2++;

				$recomb{$Mapping}{'EE'}++; 
				

			}

            
            elsif ( $Allelevarlen1 == 2 ) {
            	if(shareref($AlleleVariations1,$Reference)>0){
            		$SNP_equals_Ref1_and_Ref2++;
            		$ref1++;
            		$ref2++;
            		
            	$recomb{$Mapping}{'CE'}++;

            	}
            	elsif(shareref($AlleleVariations1,$Reference)==0){
            		$SNP_equals_Ref2_and_other++;
            		$ref2++;
            		$other++;   
            	$recomb{$Mapping}{'EE'}++;
				        	
            	}
            	
            }

                
            elsif ( $Allelevarlen1 > 2 ) {
            	if(shareref($AlleleVariations1,$Reference)>0){
            		$SNP_equals_Ref1_and_Ref2_and_other++;
            		$ref1++;
            		$ref2++;
					$other++; 
            	
				$recomb{$Mapping}{'CE'}++;
            	}
            	
            	elsif(shareref($AlleleVariations1,$Reference)==0){
            		$SNP_equals_Ref2_and_other++;
					$ref2++;
            		$other++;
            	
				$recomb{$Mapping}{'CE'}++;
            	
            	}


            
            }
            
            else { $non_processed_SNPs++;}

 }   

###########################################################
# Read file 2 for processing;
open SNP2, "$ARGV[1]" or die "cannot find $ARGV[1]";

my $totalSNP2 = 0;


my ($SNP2_equals_Ref2,
				$SNP2_equals_Ref1_and_Ref2,
            	$SNP2_equals_Ref2_and_other,
            	$SNP2_equals_Ref1_and_Ref2_and_other,
            	$ref21,
            	$ref22,
				$other2,
            	$non_processed_SNPs2
);

print ".........Done\n";

print "Reading file: $ARGV[1]...";

while (<SNP2>) {
	next if $_=~/^\s*$/;
    $totalSNP2++;
    my $blank                  = ();
    my $Mapping                = ();
    my $ReferencePosition      = ();
    my $ConsensusPosition      = ();
    my $VariationType          = ();
    my $Length                 = ();
    my $Reference              = ();
    my $Variants               = ();
    my $AlleleVariations2      = ();
    my $Frequencies            = ();
    my $Counts                 = ();
    my $Coverage               = ();
    my $Variant_1              = ();
    my $Frequencyof_1          = ();
    my $Countof_1              = ();
    my $Variant_2              = ();
    my $Frequencyof_2          = ();
    my $Countof_2              = ();
    my $OverlappingAnnotations = ();
    my $AminoAcidChange        = ();
    s/[\,]//g;
    s/[\"]+/\"/g;
    (
        $blank,                  $Mapping,       $ReferencePosition,
        $ConsensusPosition,      $VariationType, $Length,
        $Reference,              $Variants,      $AlleleVariations2,
        $Frequencies,            $Counts,        $Coverage,
        $Variant_1,              $Frequencyof_1, $Countof_1,
        $Variant_2,              $Frequencyof_2, $Countof_2,
        $OverlappingAnnotations, $AminoAcidChange
    ) = split( /\"/, $_ );

#print "\n--->Mapping: $Mapping,Refposition: $ReferencePosition,ConcensusPosition: $ConsensusPosition,Variation: $VariationType, Length: $Length, Reference: $Reference, Variants:$Variants, AlleleVariations:$AlleleVariations1, Frequencies:$Frequencies, Counts:$Counts, Coverage:$Coverage,Var1: $Variant_1, Freq1:$Frequencyof_1,Count1: $Countof_1, Var2:$Variant_2, Freq2:$Frequencyof_2, Count2:$Countof_2, OverlapAnnot:$OverlappingAnnotations, AaChange:$AminoAcidChange       \n";
    $Mapping           =~ s/\s//g;
    
	$mappingList{$Mapping}=();
	$recomb{$Mapping}{'len'}=1 if !defined $recomb{$Mapping}{'len'};
	$recomb{$Mapping}{'len'}=$ConsensusPosition if $recomb{$Mapping}{'len'} < $ConsensusPosition ;
	

	
	

        #clean and sort allele variation for comparison
        $AlleleVariations2=~ s/\W//;

        my $Allelevarlen2 = length($AlleleVariations2);

        #print "\nBefore: $AlleleVariations1\t$AlleleVariations1\n";

        #comparison of SNPs begin from here


            if ( $Allelevarlen2 == 1 ) {
				$SNP2_equals_Ref2++;
				$ref22++;

				$recomb{$Mapping}{'2EE'}++; 
				

			}

            
            elsif ( $Allelevarlen2 == 2 ) {
            	if(shareref($AlleleVariations2,$Reference)>0){
            		$SNP2_equals_Ref1_and_Ref2++;
            		$ref21++;
            		$ref22++;
            		
            	$recomb{$Mapping}{'2CE'}++;

            	}
            	elsif(shareref($AlleleVariations2,$Reference)==0){
            		$SNP2_equals_Ref2_and_other++;
            		$ref22++;
            		$other2++;   
            	$recomb{$Mapping}{'2EE'}++;
				        	
            	}
            	
            }

                
            elsif ( $Allelevarlen2 > 2 ) {
            	if(shareref($AlleleVariations2,$Reference)>0){
            		$SNP2_equals_Ref1_and_Ref2_and_other++;
            		$ref21++;
            		$ref22++;
					$other++; 
            	
				$recomb{$Mapping}{'2CE'}++;
            	}
            	
            	elsif(shareref($AlleleVariations2,$Reference)==0){
            		$SNP2_equals_Ref2_and_other++;
					$ref22++;
            		$other2++;
            	
				$recomb{$Mapping}{'2CE'}++;
            	
            	}


            
            }
            
            else { $non_processed_SNPs2++;}

 }   



print ".....Done\n";


open OUT2,">$ARGV[0].CC-EElists.txt";
 
print OUT2"RefContig\tMaxConcensusLength\t1CE\tCEper100Bp\t1EE\tEEper100Bp\t2CE\t2CEper100Bp\t2EE\t2EEper100Bp\n"; 


# set values for binning of frequency per 100bp values
#my$j=0;
my  $incr=1;
my  $max=100;
my  $init=0;
my%histogram;
print "\nCalculating frequency of EE and CE type of SNPs per 100bp\n";
foreach my$mappingNames(keys %mappingList){
		
		#$j++;
		#print "Processing $j: $mappingNames\n";

		
		if(!defined $recomb{$mappingNames}{'CE'}){$recomb{$mappingNames}{'CE'}=0;}
		if(!defined $recomb{$mappingNames}{'EE'}){$recomb{$mappingNames}{'EE'}=0;}
		if(!defined $recomb{$mappingNames}{'2CE'}){$recomb{$mappingNames}{'2CE'}=0;}
		if(!defined $recomb{$mappingNames}{'2EE'}){$recomb{$mappingNames}{'2EE'}=0;}

		$recomb{$mappingNames}{'EEper100bp'}= my$EEper100bp=$recomb{$mappingNames}{'EE'}*100/$recomb{$mappingNames}{'len'};
		$recomb{$mappingNames}{'CEper100bp'}= my$CEper100bp=$recomb{$mappingNames}{'CE'}*100/$recomb{$mappingNames}{'len'};
		
		$recomb{$mappingNames}{'EE2per100bp'}= my$EE2per100bp=$recomb{$mappingNames}{'2EE'}*100/$recomb{$mappingNames}{'len'};
		$recomb{$mappingNames}{'CE2per100bp'}= my$CE2per100bp=$recomb{$mappingNames}{'2CE'}*100/$recomb{$mappingNames}{'len'};
		
		$recomb{$mappingNames}{'EEtoEE2ratio'}=($recomb{$mappingNames}{'EE'})/($recomb{$mappingNames}{'2EE'}+1);
		$recomb{$mappingNames}{'CEtoCE2ratio'}=($recomb{$mappingNames}{'CE'})/($recomb{$mappingNames}{'2CE'}+1);

		
		
	print OUT2"$mappingNames\t$recomb{$mappingNames}{'len'}\t$recomb{$mappingNames}{'CE'}\t$CEper100bp\t$recomb{$mappingNames}{'EE'}\t$EEper100bp\t$recomb{$mappingNames}{'2CE'}\t$CE2per100bp\t$recomb{$mappingNames}{'2EE'}\t$EE2per100bp\n";

	my$k=0;
	for(my$i=$init;$i < $max; $i=$i+$incr){
		#$k++;
		#print"$mappingNames:\t$k\n";
		if($recomb{$mappingNames}{'EEper100bp'}>$i && $recomb{$mappingNames}{'EEper100bp'}<$i+$incr ){
			$histogram{$i}{'EEper100bp'}++;
		}
		
		if($recomb{$mappingNames}{'CEper100bp'}>$i && $recomb{$mappingNames}{'CEper100bp'}<$i+$incr ){
			$histogram{$i}{'CEper100bp'}++;
		}

		if($recomb{$mappingNames}{'EE2per100bp'}>$i && $recomb{$mappingNames}{'EE2per100bp'}<$i+$incr ){
			$histogram{$i}{'EE2per100bp'}++;
		}
		
		if($recomb{$mappingNames}{'CE2per100bp'}>$i && $recomb{$mappingNames}{'CE2per100bp'}<$i+$incr ){
			$histogram{$i}{'CE2per100bp'}++;
		}
		
		
		if($recomb{$mappingNames}{'EEtoEE2ratio'}>$i && $recomb{$mappingNames}{'EEtoEE2ratio'}<$i+$incr ){
			$histogram{$i}{'EEtoEE2ratio'}++;
		}
		
		if($recomb{$mappingNames}{'CEtoCE2ratio'}>$i && $recomb{$mappingNames}{'CEtoCE2ratio'}<$i+$incr ){
			$histogram{$i}{'CEtoCE2ratio'}++;
		}

	
	}
	
	
}

 
#my$incr=1;
#my$max=1;
#my$init=0;
#my%histogram;
#for(my$i=$init;$i<=$max;$i+$incr){
#
#	foreach my$contigNames(keys %mappingList){
#		if($recomb{$contigNames}{'EEper100bp'}>$i && $recomb{$contigNames}{'EEper100bp'}<$i+$incr ){
#			$histogram{$i}{'EEper100bp'}++;
#		}
#		
#		if($recomb{$contigNames}{'CEper100bp'}>$i && $recomb{$contigNames}{'CEper100bp'}<$i+$incr ){
#			$histogram{$i}{'CEper100bp'}++;
#		}
#
#		if($recomb{$contigNames}{'EE2per100bp'}>$i && $recomb{$contigNames}{'EE2per100bp'}<$i+$incr ){
#			$histogram{$i}{'EE2per100bp'}++;
#		}
#		
#		if($recomb{$contigNames}{'CE2per100bp'}>$i && $recomb{$contigNames}{'CE2per100bp'}<$i+$incr ){
#			$histogram{$i}{'CE2per100bp'}++;
#		}
#
#	}
#}

open OUT3,">$ARGV[0].frequencyHistogram.txt" or die "Cannot open output file for Frequecny table\n";
open OUT4,">$ARGV[0].EEandCEratios.txt" or die "Cannot open output file for Frequecny table\n";

print OUT3"SNPPer100Bp\tFreq_EE\tFreq_CE\tFreq_2EE\tFreq_2CE\n";
print OUT4"ratios\tFreq_EEto2EE\tFreq_CEto2CE\n";

for(my$i=$init;$i<=$max;$i=$i+$incr){
	$histogram{$i}{'EEper100bp'}=0 if !defined $histogram{$i}{'EEper100bp'};
	$histogram{$i}{'CEper100bp'}=0 if !defined $histogram{$i}{'CEper100bp'};
	$histogram{$i}{'EE2per100bp'}=0 if !defined $histogram{$i}{'EE2per100bp'};
	$histogram{$i}{'CE2per100bp'}=0 if !defined $histogram{$i}{'CE2per100bp'};
	print OUT3 round($i)."\t$histogram{$i}{'EEper100bp'}\t$histogram{$i}{'CEper100bp'}\t$histogram{$i}{'EE2per100bp'}\t$histogram{$i}{'CE2per100bp'}\n";
	print OUT4 round($i)."\t$histogram{$i}{'EEtoEE2ratio'}\t $histogram{$i}{'CEtoCE2ratio'}\n";
}











 
# 
# print"
#TotalSNP1:$totalSNP1\t".round($totalSNP1*100/$totalSNP1)."
#SNP_equals_Ref2:$SNP_equals_Ref2\t".round($SNP_equals_Ref2*100/$totalSNP1)."
#SNP_equals_Ref1_and_Ref2:$SNP_equals_Ref1_and_Ref2\t".round($SNP_equals_Ref1_and_Ref2*100/$totalSNP1)."
#SNP_equals_Ref2_and_other:$SNP_equals_Ref2_and_other\t".round($SNP_equals_Ref2_and_other*100/$totalSNP1)."
#SNP_equals_Ref1_and_Ref2_and_other:$SNP_equals_Ref1_and_Ref2_and_other\t".round($SNP_equals_Ref1_and_Ref2_and_other*100/$totalSNP1)."
#
#\n\n\nBased on total number fo nucleotides compared. Note that Each SNP loci can have more than one nucleotides\n
#
#Ref1:$ref1\t".round($ref1*100/($ref1+$ref2+$other))."
#Ref2:$ref2\t".round($ref2*100/($ref1+$ref2+$other))."
#Others:$other\t".round($other*100/($ref1+$ref2+$other))."
#non-processed_SNPs:\t".round($non_processed_SNPs*100/($ref1+$ref2+$other))."
# ";
#     

####################################################################
# subroutines

sub shareref {
    my $alleles1 = my $ref = ();
    ( $alleles1, $ref ) = @_;
    my $shareref = 0;
    $alleles1 =~ s/\W//g;
    $ref      =~ s/\W//g;

    my @alleles1 = split( undef, $alleles1 );
    my @ref      = split( undef, $ref );

    foreach my $al1 (@alleles1) {
        foreach my $al2 (@ref) {
            if ( $al1 eq $al2 ) {
                $shareref++;
            }
        }
    }

    return ($shareref);
}

sub sharent {
    my $alleles1 = my $alleles2 = ();
    ( $alleles1, $alleles2 ) = @_;
    my $sharent = 0;
    $alleles1 =~ s/\W//g;
    $alleles2 =~ s/\W//g;

    my @alleles1 = split( undef, $alleles1 );
    my @alleles2 = split( undef, $alleles2 );

    foreach my $al1 (@alleles1) {
        foreach my $al2 (@alleles2) {
            if ( $al1 eq $al2 ) {
                $sharent++;
            }
        }
    }

    return ($sharent)

}



sub round{
	my $decimals=3;
    my( $number) = @_;
    substr( $number + ( '0.' . '0' x $decimals . '5' ), 0, $decimals + length(int($number)) + 1 );
}

