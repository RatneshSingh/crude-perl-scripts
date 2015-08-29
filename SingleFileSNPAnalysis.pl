#/usr/perl/bin -w
use strict;
use Time::localtime;

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



while (<SNP1>) {

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
    


        #clean and sort allele variation for comparison
        $AlleleVariations1=~ s/\W//;

        my $Allelevarlen1 = length($AlleleVariations1);

        #print "\nBefore: $AlleleVariations1\t$AlleleVariations1\n";

        #comparison of SNPs begin from here


            if ( $Allelevarlen1 == 1 ) {
				$SNP_equals_Ref2++;
				$ref2++;
            }

            
            elsif ( $Allelevarlen1 == 2 ) {
            	if(shareref($AlleleVariations1,$Reference)>0){
            		$SNP_equals_Ref1_and_Ref2++;
            		$ref1++;
            		$ref2++;
            	}
            	elsif(shareref($AlleleVariations1,$Reference)==0){
            		$SNP_equals_Ref2_and_other++;
            		$ref2++;
            		$other++;            	
            	}
            	
            }

                
            elsif ( $Allelevarlen1 > 2 ) {
            	if(shareref($AlleleVariations1,$Reference)>0){
            		$SNP_equals_Ref1_and_Ref2_and_other++;
            		$ref1++;
            		$ref2++;
					$other++; 
            	}
            	
            	elsif(shareref($AlleleVariations1,$Reference)==0){
            		$SNP_equals_Ref2_and_other++;
					$ref2++;
            		$other++;
            	
            	}


            
            }
            
            else { $non_processed_SNPs++;}

 }   
 
 
 
 print"
TotalSNP1:$totalSNP1\t".round($totalSNP1*100/$totalSNP1)."
SNP_equals_Ref2:$SNP_equals_Ref2\t".round($SNP_equals_Ref2*100/$totalSNP1)."
SNP_equals_Ref1_and_Ref2:$SNP_equals_Ref1_and_Ref2\t".round($SNP_equals_Ref1_and_Ref2*100/$totalSNP1)."
SNP_equals_Ref2_and_other:$SNP_equals_Ref2_and_other\t".round($SNP_equals_Ref2_and_other*100/$totalSNP1)."
SNP_equals_Ref1_and_Ref2_and_other:$SNP_equals_Ref1_and_Ref2_and_other\t".round($SNP_equals_Ref1_and_Ref2_and_other*100/$totalSNP1)."

\n\n\nBased on total number fo nucleotides compared. Note that Each SNP loci can have more than one nucleotides\n

Ref1:$ref1\t".round($ref1*100/($ref1+$ref2+$other))."
Ref2:$ref2\t".round($ref2*100/($ref1+$ref2+$other))."
Others:$other\t".round($other*100/($ref1+$ref2+$other))."
non-processed_SNPs:\t".round($non_processed_SNPs*100/($ref1+$ref2+$other))."
 ";
     

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

