#/usr/perl/bin -w
use strict;
use Time::localtime;

open LOG, ">>$0.Log.txt";
print LOG ctime(), "\n@ARGV\n";
print ctime(), "\n$0 @ARGV\n";

#my($Mapping, $ReferencePosition, $ConsensusPosition, $VariationType, $Length, $Reference, $Variants, $AlleleVariations2, $Frequencies, $Counts, $Coverage, $Variant_1, $Frequencyof_1, $Countof_1, $Variant_2, $Frequencyof_2, $Countof_2, $OverlappingAnnotations, $AminoAcidChange);
my (
    %SNP,         %SNP1,             %SNP2,
    %snp_mapping, %RefPosition,      $SNP1count,
    $SNP2count,   $SNP2_1_bothcount, $SNP1_2_bothcount
);

open SNP1, "$ARGV[0]" or die "cannot find $ARGV[0] 1st SNP file";
open SNP2, "$ARGV[1]" or die "cannot find $ARGV[1] 2nd SNP file";
my $totalSNP1 = 0;

#my@SNP1_list;
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
	s/^\s+//g;
    (
        $Mapping,       $ReferencePosition,
        $ConsensusPosition,      $VariationType, $Length,
        $Reference,              $Variants,      $AlleleVariations1,
        $Frequencies,            $Counts,        $Coverage,
        $Variant_1,              $Frequencyof_1, $Countof_1,
        $Variant_2,              $Frequencyof_2, $Countof_2,
        $OverlappingAnnotations, $AminoAcidChange
    ) = split( /\t/, $_ );

 ##Test print
 #print "\nblank :$blank Mapping :$Mapping ReferencePosition :$ReferencePosition ConsensusPosition :$ConsensusPosition VariationType :$VariationType Length :$Length Reference :$Reference Variants :$Variants AlleleVariations1 :$AlleleVariations1 Frequencies :$Frequencies Counts :$Counts Coverage :$Coverage Variant_1 :$Variant_1 Frequencyof_1 :$Frequencyof_1 Countof_1 :$Countof_1 Variant_2 :$Variant_2 Frequencyof_2 :$Frequencyof_2 Countof_2 :$Countof_2 OverlappingAnnotations :$OverlappingAnnotations AminoAcidChange :$AminoAcidChange\n";
 #print "\n--->Mapping: $Mapping,Refposition: $ReferencePosition,ConcensusPosition: $ConsensusPosition,Variation: $VariationType, Length: $Length, Reference: $Reference, Variants:$Variants, AlleleVariations:$AlleleVariations2, Frequencies:$Frequencies, Counts:$Counts, Coverage:$Coverage,Var1: $Variant_1, Freq1:$Frequencyof_1,Count1: $Countof_1, Var2:$Variant_2, Freq2:$Frequencyof_2, Count2:$Countof_2, OverlapAnnot:$OverlappingAnnotations, AaChange:$AminoAcidChange       \n";
    $Mapping           =~ s/\s+//g;
    $ReferencePosition =~ s/\D+//g;

    $SNP1{$Mapping}                  = ();
    $RefPosition{$ReferencePosition} = ();

    # my$snpline=$Mapping1.$ReferencePosition1;
    #$snpline=~s/\s//g;

#push(@SNP1_list,$snpline);
#$SNP1{$snpline}=();
#
#$SNP{$Mapping}=$Mapping; # used $SNP{} instead of $SNP1 or2 to reduce the number of hash elements. They will be shared in both the files
#$SNP1{$ReferencePosition}=$ReferencePosition;
#$SNP1{$Mapping.$ReferencePosition}{'ConsensusPosition'}=$ConsensusPosition;
#$SNP1{$Mapping.$ReferencePosition}{'VariationType'}=$VariationType;
#$SNP1{$Mapping.$ReferencePosition}{'Length'}=$Length;
    $SNP1{ $Mapping . $ReferencePosition }{'Reference'} = $Reference;

    #$SNP1{$Mapping.$ReferencePosition}{'Variants'}=$Variants;
    $SNP1{ $Mapping . $ReferencePosition }{'AlleleVariations'} =
      $AlleleVariations1;
    $SNP1{ $Mapping . $ReferencePosition }{'Frequencies'} = $Frequencies;

#$SNP1{$Mapping.$ReferencePosition}{'Counts'}=$Counts;
#$SNP1{$Mapping.$ReferencePosition}{'Coverage'}=$Coverage;
#$SNP1{$Mapping.$ReferencePosition}{'Variant_1'}=$Variant_1;
#$SNP1{$Mapping.$ReferencePosition}{'Frequencyof_1'}=$Frequencyof_1;
#$SNP1{$Mapping.$ReferencePosition}{'Countof_1'}=$Countof_1;
#$SNP1{$Mapping.$ReferencePosition}{'Variant_2'}=$Variant_2;
#$SNP1{$Mapping.$ReferencePosition}{'Frequencyof_2'}=$Frequencyof_2;
#$SNP1{$Mapping.$ReferencePosition}{'Countof_2'}=$Countof_2;
#$SNP1{$Mapping.$ReferencePosition}{'OverlappingAnnotations'}=$OverlappingAnnotations;
#$SNP1{$Mapping.$ReferencePosition}{'AminoAcidChange'}=$AminoAcidChange;

}

close SNP1;
open OUT, ">Canephora-Mokka-Catimor_SNPs.txt";
print OUT
"Mapping\tReference Position\tReference\t$ARGV[0] Allele Variation\tFrequencies\t$ARGV[1] Allele Variation\tFrequencies\n";

my (
    $All1_eq_All2_and_len1_ofcourseNoRef1_gt1, $CaneEugenSNP1,
    $CaneEugenSNP2,                            $SNP_DiffIn1n2nRef,
	$totalSNP2,


);

$totalSNP2 = 0;

#my@SNP2_list;

my (
    $totalCommonSNP,

    $All1_eq_All2_and_len1_ofcourseNoRef1,

    $A1eqA2_len2_noRef,

    $A1eqA2_len2_yesRef1,

    $A1eqA2_len2up_noRef1,

    $A1eqA2_len2up_yesRef1,

    $A1neqA2_len2_yesRef1,

    $A1neqA2_len2_sharent_noRef1,

    $A1neqA2_len2_sharent_OneShareRef1,

    $A1neqA2_len2plus_yesRef1,

    $A1neqA2_len2plus_sharent_noRef1,

    $A1neqA2_len2plus_sharent_OneShareRef1,

    $A1neqA2_oneslengthis1_sharent_BothShareRef1,

	$A1neqA2_oneslengthis1_sharent_OnlyOneShareRef1,

    $A1neqA2_oneslengthis1_sharent_NoneShareRef1,

	$A1neqA2_All1_lengthis1_sharent_All2_ShareRef1,

	$A1neqA2_All2_lengthis1_sharent_All1_ShareRef1,

	$A1neqA2_oneslengthis1_sharent_unprocessed,

    $A1neqA2_len1,

    $A1neqA2_len2orup_Nosharent_ShareRef1,

    $A1neqA2_len2orup_Nosharent_OnlyOneShareRef1,

    $nonprocessedCommonSNP,

    $SNP_inOnlyfile2_Ref1_infile1,

    $SNP_inOnlyfile2_ShareRef1_Ref1_infile1,

    $SNP_inOnlyfile2_DontShareRef1_Ref1_infile1,

    $nonprocessedSNP,

    $A1eqA2_len2_unprocessed,
);

while (<SNP2>) {

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
	s/^\s+//g;
    (
        $Mapping,       $ReferencePosition,
        $ConsensusPosition,      $VariationType, $Length,
        $Reference,              $Variants,      $AlleleVariations2,
        $Frequencies,            $Counts,        $Coverage,
        $Variant_1,              $Frequencyof_1, $Countof_1,
        $Variant_2,              $Frequencyof_2, $Countof_2,
        $OverlappingAnnotations, $AminoAcidChange
    ) = split( /\t/, $_ );

#print "\n--->Mapping: $Mapping,Refposition: $ReferencePosition,ConcensusPosition: $ConsensusPosition,Variation: $VariationType, Length: $Length, Reference: $Reference, Variants:$Variants, AlleleVariations:$AlleleVariations2, Frequencies:$Frequencies, Counts:$Counts, Coverage:$Coverage,Var1: $Variant_1, Freq1:$Frequencyof_1,Count1: $Countof_1, Var2:$Variant_2, Freq2:$Frequencyof_2, Count2:$Countof_2, OverlapAnnot:$OverlappingAnnotations, AaChange:$AminoAcidChange       \n";
    $Mapping           =~ s/\s+//g;
    $ReferencePosition =~ s/\D+//g;

    $SNP2{$Mapping}                  = ();
    $RefPosition{$ReferencePosition} = ();

    # my$snpline=$Mapping2.$ReferencePosition2;
    # $snpline=~s/\s//g;
    #push(@SNP2_list,$snpline);

    # $SNP2{$snpline}=();

    #hash eats up too much memory. use other method instead
    #$SNP{$Mapping}=$Mapping;
    #$SNP2{$ReferencePosition}=$ReferencePosition;
    #$SNP2{$Mapping.$ReferencePosition}{'ConsensusPosition'}=$ConsensusPosition;
    #$SNP2{$Mapping.$ReferencePosition}{'VariationType'}=$VariationType;
    #$SNP2{$Mapping.$ReferencePosition}{'Length'}=$Length;
    $SNP2{ $Mapping . $ReferencePosition }{'Reference'} = $Reference;

    #$SNP2{$Mapping.$ReferencePosition}{'Variants'}=$Variants;
    $SNP2{ $Mapping . $ReferencePosition }{'AlleleVariations'} =
      $AlleleVariations2;
    $SNP2{ $Mapping . $ReferencePosition }{'Frequencies'} = $Frequencies;

#$SNP2{$Mapping.$ReferencePosition}{'Counts'}=$Counts;
#$SNP2{$Mapping.$ReferencePosition}{'Coverage'}=$Coverage;
#$SNP2{$Mapping.$ReferencePosition}{'Variant_1'}=$Variant_1;
#$SNP2{$Mapping.$ReferencePosition}{'Frequencyof_1'}=$Frequencyof_1;
#$SNP2{$Mapping.$ReferencePosition}{'Countof_1'}=$Countof_1;
#$SNP2{$Mapping.$ReferencePosition}{'Variant_2'}=$Variant_2;
#$SNP2{$Mapping.$ReferencePosition}{'Frequencyof_2'}=$Frequencyof_2;
#$SNP2{$Mapping.$ReferencePosition}{'Countof_2'}=$Countof_2;
#$SNP2{$Mapping.$ReferencePosition}{'OverlappingAnnotations'}=$OverlappingAnnotations;
#$SNP2{$Mapping.$ReferencePosition}{'AminoAcidChange'}=$AminoAcidChange;

    if ( exists $SNP1{ $Mapping . $ReferencePosition }{'Reference'} ) {

        #count total SNP compared
        $totalCommonSNP++;

        #clean and sort allele variation for comparison
        #$AlleleVariations2=~ s/\W//;
		$AlleleVariations2=~ s/\///;
        my $AlleleVariations1 =
          $SNP1{ $Mapping . $ReferencePosition }{'AlleleVariations'};
        $AlleleVariations1 =~ s/\///;


        my $Allelevarlen1 = length($AlleleVariations1);
        my $Allelevarlen2  = length($AlleleVariations2);

        #print "\nBefore: $AlleleVariations2\t$AlleleVariations1\n";

        # if any of position has more than one allele variation, sort them alphabetically for easy comparison.
        if ( $Allelevarlen2 > 1 || $Allelevarlen1 > 1 ) {
            my @allvar2u = split( //, $AlleleVariations2);
            my @allvar1u = split( //, $AlleleVariations1 );

            my @allvar2s = sort { lc($a) cmp lc($b) } @allvar2u;
            my @allvar1s = sort { lc($a) cmp lc($b) } @allvar1u;

            $AlleleVariations2 = join( "", @allvar2s );
            $AlleleVariations1 = join( "", @allvar1s );

            #print "After :$AlleleVariations2\t$AlleleVariations1\n";
        }
		################################################################################################################################################################
        #comparison of SNPs begin from here
		################################################################################################################################################################
        #
		# if both genomes have same alleles at all1=all2
        if ( $AlleleVariations1 eq $AlleleVariations2) {

            # all1=all2 and length-All1-N-All2 = 1.
            if ( $Allelevarlen1 == 1 && $Allelevarlen2 == 1 && shareref( $AlleleVariations1, $Reference ) == 0 && shareref( $AlleleVariations2, $Reference ) == 0) {
                $All1_eq_All2_and_len1_ofcourseNoRef1++;
            }

            # all1=all2 and length-All1-N-All2 = 2.
            elsif ( $Allelevarlen1 == 2 && $Allelevarlen2 == 2 ) {

                # all1=all2 and length-All1-N-All2 = 2. (all1=all2)!=~ ref1
                if (   shareref( $AlleleVariations1, $Reference ) == 0
                    && shareref( $AlleleVariations2, $Reference ) == 0 )
                {

# all1 and all2 are equal and do not share ref1 but may share ref2 alongwith other genome.
                    $A1eqA2_len2_noRef++;
                }

                # all1=all2 and length-All1-N-All2 = 2. (all1=all2)=~ ref1
                elsif (shareref( $AlleleVariations1, $Reference ) >= 1
                    && shareref( $AlleleVariations2, $Reference ) >= 1 )
                {

                    # allele1 and allele2 share ref1 and other is ref 2
                    $A1eqA2_len2_yesRef1++;
                }

                else { $A1eqA2_len2_unprocessed++; }
            }

            # all1=all2 and length-All1-N-All2 > 2.
            elsif ( $Allelevarlen1 > 2 && $Allelevarlen2 > 2 ) {

                # all1=all2 and length-All1-N-All2 >2. (all1=all2)!=~ ref1
                if (   shareref( $AlleleVariations1, $Reference ) == 0
                    && shareref( $AlleleVariations2, $Reference ) == 0 )
                {

                   # Do not share ref1 and may share ref2 but have other genomes
                    $A1eqA2_len2up_noRef1++;
                }

                # all1=all2 and length-All1-N-All2 >2. (all1=all2)=~ ref1
                elsif (shareref( $AlleleVariations1, $Reference ) >= 1
                    && shareref( $AlleleVariations2, $Reference ) >= 1 )
                {

                    # share ref1 and may be ref2 but have other genome
                    $A1eqA2_len2up_yesRef1;
                }

            }

            else {
                print
"\nError.$AlleleVariations1 And $AlleleVariations2 are not considered equal\n ";
            }

        }

        # all1 != all2
        elsif ( $AlleleVariations1 ne $AlleleVariations2) {

            #print "- $AlleleVariations1 \t $AlleleVariations2\n";

            # all1 != all2, All1 =~sharent  All2
            if ( sharent( $AlleleVariations1, $AlleleVariations2) >= 1 ) {

                # all1 != all2, All1 =~sharent  All2, length(all1)=lengthAll2
                if ( $Allelevarlen1 == $Allelevarlen2 ) {

                    if ( $Allelevarlen1 == 2 && $Allelevarlen2 == 2 ) {

  # all1 != all2, All1 =~sharent  All2, length(all1)=lengthAll2=2, both shareref
                        if (   shareref( $AlleleVariations1, $Reference ) >= 1
                            && shareref( $AlleleVariations2, $Reference ) >= 1 )
                        {

                            # share ref1 and may be ref2
                            $A1neqA2_len2_yesRef1++;
                        }

  # all1 != all2, All1 =~sharent  All2, length(all1)=lengthAll2, both ! shareref
                        elsif (shareref( $AlleleVariations1, $Reference ) == 0
                            && shareref( $AlleleVariations2, $Reference ) == 0 )
                        {

      #  do not share ref1 and one might share ref2 while other has other genome
                            $A1neqA2_len2_sharent_noRef1++;
                        }

# all1 != all2, All1 =~sharent  All2, length(all1)=lengthAll2=2, atleast one shareref and other dont
                        elsif (
                            (
                                shareref( $AlleleVariations1, $Reference ) == 0
                                && shareref( $AlleleVariations2, $Reference ) > 0
                            )
                            || ( shareref( $AlleleVariations1, $Reference ) > 0
                                && shareref( $AlleleVariations2, $Reference ) ==
                                0 )
                          )
                        {

       #  one shares ref1 and other do not. atleast one has part of other genome
                            $A1neqA2_len2_sharent_OneShareRef1++;
                        }

                    }

                    elsif ( $Allelevarlen1 > 2 && $Allelevarlen2 > 2 ) {

 # all1 != all2, All1 =~sharent  All2, length(all1)=lengthAll2>=2, both shareref
                        if (   shareref( $AlleleVariations1, $Reference ) >= 1
                            && shareref( $AlleleVariations2, $Reference ) >= 1 )
                        {

                            # share ref1 and may be ref2
                            $A1neqA2_len2plus_yesRef1++;
                        }

  # all1 != all2, All1 =~sharent  All2, length(all1)=lengthAll2, both ! shareref
                        elsif (shareref( $AlleleVariations1, $Reference ) == 0
                            && shareref( $AlleleVariations2, $Reference ) == 0 )
                        {

      #  do not share ref1 and one might share ref2 while other has other genome
                            $A1neqA2_len2plus_sharent_noRef1++;
                        }

# all1 != all2, All1 =~sharent  All2, length(all1)=lengthAll2=2, atleast one shareref and other dont
                        elsif (
                            (
                                shareref( $AlleleVariations1, $Reference ) == 0
                                && shareref( $AlleleVariations2, $Reference ) > 0
                            )
                            || ( shareref( $AlleleVariations1, $Reference ) > 0
                                && shareref( $AlleleVariations2, $Reference ) ==
                                0 )
                          )
                        {

       #  one shares ref1 and other do not. atleast one has part of other genome
                            $A1neqA2_len2plus_sharent_OneShareRef1++;
                        }

                    }
                }

                # all1 != all2, All1 =~sharent  All2, length(all1) !=lengthAll2
                elsif ( $Allelevarlen1 != $Allelevarlen2 ) {

                    if (   shareref( $AlleleVariations2, $Reference ) > 0
                        && shareref( $AlleleVariations1, $Reference ) > 0 )
                    {
                        # all1 != all2, All1 =~sharent  All2, length(all1) !=lengthAll2, one length is 1

                        if ( $Allelevarlen1 == 1 || $Allelevarlen2 == 1 ) {

                            $A1neqA2_oneslengthis1_sharent_BothShareRef1++;

                        }
                    }


				   elsif (shareref( $AlleleVariations2, $Reference ) > 0
                        || shareref( $AlleleVariations1, $Reference ) > 0 )
                    {
						## all1 != all2, All1 =~sharent  All2, length(all1) !=lengthAll2, one length is 1, Only one share ref1
                        if ( $Allelevarlen1 == 1 || $Allelevarlen2 == 1 ) {

                            $A1neqA2_oneslengthis1_sharent_OnlyOneShareRef1++;

                        	if($Allelevarlen1 == 1 && shareref( $AlleleVariations2, $Reference ) > 0 && shareref( $AlleleVariations1, $Reference ) == 0){

								#Allele2 share Ref1, and Allele1 is 1nt long.
								$A1neqA2_All1_lengthis1_sharent_All2_ShareRef1++

							}

							elsif($Allelevarlen2 == 1 && shareref( $AlleleVariations1, $Reference ) > 0 && shareref( $AlleleVariations2, $Reference ) == 0){

								#Allele1 share Ref1, and Allele2 is 1 nt long.
								$A1neqA2_All2_lengthis1_sharent_All1_ShareRef1++


							}

						}







					}

                    elsif (shareref( $AlleleVariations2, $Reference ) == 0
                        && shareref( $AlleleVariations1, $Reference ) == 0 )
                    {

                        $A1neqA2_oneslengthis1_sharent_NoneShareRef1++;

                    }

                    else {
                        $A1neqA2_oneslengthis1_sharent_unprocessed++;
                    }

                }

            }

            # all1 != all2, All1 !=~sharent  All2
            elsif ( sharent( $AlleleVariations1, $AlleleVariations2) == 0 ) {

                # all1 != all2, All1 !=~sharent  All2, length(all1)=lengthAll2=1
                if ( $Allelevarlen1 == 1 && $Allelevarlen2 == 1 ) {

#all1 and all2 both are different and do not share ref1 at least one is other genome
                    $A1neqA2_len1++;
                }

               # all1 != all2, All1 !=~sharent  All2, length(all1)=lengthAll2>=2
                elsif ( $Allelevarlen2 >= 2 && $Allelevarlen1 >= 2 ) {

# all1 != all2, All1 !=~sharent  All2, length(all1)=lengthAll2>=2, All1-All2 =~shareref
                    if (   shareref( $AlleleVariations2, $Reference ) >= 1
                        && shareref( $AlleleVariations1, $Reference ) >= 1 )
                    {

                        #share ref2 and have other genomes
                        $A1neqA2_len2orup_Nosharent_ShareRef1++;
                    }

# all1 != all2, All1 !=~sharent  All2, length(all1)=lengthAll2>=2, All1 !~sharent  All2
                    elsif (shareref( $AlleleVariations2, $Reference ) >= 1
                        || shareref( $AlleleVariations1, $Reference ) )
                    {

                        #share ref2 and have other genomes
                        $A1neqA2_len2orup_Nosharent_OnlyOneShareRef1++;
                    }
                }

            }

        }

        # non processed SNPs
        else {
            $nonprocessedCommonSNP++;

        }

#print OUT"$Mapping\t$ReferencePosition\t$SNP1{$Mapping.$ReferencePosition}{'Reference'}\t$SNP1{$Mapping.$ReferencePosition}{'AlleleVariations'}\t$SNP1{$Mapping.$ReferencePosition}{'Frequencies'}\t$AlleleVariations2\t$Frequencies\n"
    }

    elsif ( exists $SNP1{$Mapping} && exists $SNP2{$Mapping} ) {

        #SNPs in flile2 but ref1 in file1
        $SNP_inOnlyfile2_Ref1_infile1++;

        if ( shareref( $AlleleVariations2, $Reference ) > 0 ) {
            $SNP_inOnlyfile2_ShareRef1_Ref1_infile1++;
        }

        if ( shareref( $AlleleVariations2, $Reference ) == 0 ) {
            $SNP_inOnlyfile2_DontShareRef1_Ref1_infile1++;
        }

    }

    else { $nonprocessedSNP++; }

}

print "\n
Total SNPs in $ARGV[0]:$totalSNP1

Total SNPs in $ARGV[1]:$totalSNP2

totalCommonSNPcompared:$totalCommonSNP\n\t"
. eval{round($totalCommonSNP*100/$totalSNP1)}."   of File1\n\t"
. eval{round($totalCommonSNP*100/$totalSNP2)}."   of File2\n\n"
.
"****Allele1 is equal to Allele2****\n
All1_eq_All2_and_len1_OfcourseNoRef1 mean Ref2:$All1_eq_All2_and_len1_ofcourseNoRef1  "
  .  eval{round($All1_eq_All2_and_len1_ofcourseNoRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($All1_eq_All2_and_len1_ofcourseNoRef1 * 100 /$totalSNP2)}.

  "\nA1eqA2_len2_noRef1 means have ref2 and other genome:$A1eqA2_len2_noRef  "
  .  eval{round($A1eqA2_len2_noRef * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1eqA2_len2_noRef * 100 /$totalSNP2)}.

  "\nA1eqA2_len2_yesRef1:$A1eqA2_len2_yesRef1  "
  . eval{ round($A1eqA2_len2_yesRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1eqA2_len2_yesRef1 * 100 /$totalSNP2)}.

  "\n	A1eqA2_len2_unprocessed:$A1eqA2_len2_unprocessed  "
  .  eval{round($A1eqA2_len2_unprocessed * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1eqA2_len2_unprocessed * 100 / $totalSNP2)}.

  "\nA1eqA2_len2up_noRef1:$A1eqA2_len2up_noRef1  "
  .  eval{round($A1eqA2_len2up_noRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1eqA2_len2up_noRef1 * 100 / $totalSNP2)}.

  "\nA1eqA2_len2up_yesRef1 Mean yesRef2 too:$A1eqA2_len2up_yesRef1  "
  .  eval{round($A1eqA2_len2up_yesRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1eqA2_len2up_yesRef1 * 100 / $totalSNP2)}.

  "\n\n\n***Allele is not equal to Allele2****\n
A1neqA2_len2_yesRef1:$A1neqA2_len2_yesRef1  "
  .  eval{round($A1neqA2_len2_yesRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_len2_yesRef1 * 100 / $totalSNP2)}.

  "\nA1neqA2_len2_sharent_noRef1:$A1neqA2_len2_sharent_noRef1  "
  .  eval{round($A1neqA2_len2_sharent_noRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_len2_sharent_noRef1 * 100 / $totalSNP2)}.

  "\nA1neqA2_len2_sharent_OneShareRef1:$A1neqA2_len2_sharent_OneShareRef1  "
  .  eval{round($A1neqA2_len2_sharent_OneShareRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_len2_sharent_OneShareRef1 * 100 /$totalSNP2)}.

  "\nA1neqA2_len2plus_yesRef1:$A1neqA2_len2plus_yesRef1  "
  .  eval{round($A1neqA2_len2plus_yesRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_len2plus_yesRef1 * 100 / $totalSNP2)}.

  "\nA1neqA2_len2plus_sharent_noRef1:$A1neqA2_len2plus_sharent_noRef1  "
  .  eval{round($A1neqA2_len2plus_sharent_noRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_len2plus_sharent_noRef1 * 100 / $totalSNP2)}.

"\nA1neqA2_len2plus_sharent_OneShareRef1:$A1neqA2_len2plus_sharent_OneShareRef1  "
  .  eval{round($A1neqA2_len2plus_sharent_OneShareRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_len2plus_sharent_OneShareRef1 * 100 / $totalSNP2)}.

"\nA1neqA2_oneslengthis1_sharent_BothShareRef1:$A1neqA2_oneslengthis1_sharent_BothShareRef1  "
  .  eval{round($A1neqA2_oneslengthis1_sharent_BothShareRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_oneslengthis1_sharent_BothShareRef1 * 100 /$totalSNP2)}.

"\n\nA1neqA2_oneslengthis1_sharent_OnlyOneShareRef1:$A1neqA2_oneslengthis1_sharent_OnlyOneShareRef1  "
  .  eval{round($A1neqA2_oneslengthis1_sharent_OnlyOneShareRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_oneslengthis1_sharent_OnlyOneShareRef1 * 100 /$totalSNP2)}.

  "\n\tA1neqA2_All1_lengthis1_sharent_All2_ShareRef1:$A1neqA2_All1_lengthis1_sharent_All2_ShareRef1  "
  . eval{round($A1neqA2_All1_lengthis1_sharent_All2_ShareRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_All1_lengthis1_sharent_All2_ShareRef1 * 100 / $totalSNP2)}.

  "\n\tA1neqA2_All2_lengthis1_sharent_All1_ShareRef1:$A1neqA2_All2_lengthis1_sharent_All1_ShareRef1  "
  . eval{round($A1neqA2_All2_lengthis1_sharent_All1_ShareRef1 * 100/$totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_All2_lengthis1_sharent_All1_ShareRef1 * 100/$totalSNP2)}.

  "\n\nA1neqA2_oneslengthis1_sharent_NoneShareRef1:$A1neqA2_oneslengthis1_sharent_NoneShareRef1  "
  .  eval{round($A1neqA2_oneslengthis1_sharent_NoneShareRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_oneslengthis1_sharent_NoneShareRef1 * 100 / $totalSNP2)}.

"\n	A1neqA2_oneslengthis1_sharent_unprocessed:$A1neqA2_oneslengthis1_sharent_unprocessed	"
  .  eval{round($A1neqA2_oneslengthis1_sharent_unprocessed * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_oneslengthis1_sharent_unprocessed * 100 /$totalSNP2)}.

"\nA1neqA2_len1 ofcourse don't sharent, don't share ref1,one is other genome:$A1neqA2_len1  "
  .  eval{round($A1neqA2_len1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_len1 * 100 /$totalSNP2)}.

"\nA1neqA2_len2orup_Nosharent_ShareRef1:$A1neqA2_len2orup_Nosharent_ShareRef1  "
  .  eval{round($A1neqA2_len2orup_Nosharent_ShareRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_len2orup_Nosharent_ShareRef1 * 100 / $totalSNP2)}.

"\nA1neqA2_len2orup_Nosharent_OnlyOneShareRef1:$A1neqA2_len2orup_Nosharent_OnlyOneShareRef1  "
  .  eval{round($A1neqA2_len2orup_Nosharent_OnlyOneShareRef1 * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($A1neqA2_len2orup_Nosharent_OnlyOneShareRef1 * 100 / $totalSNP2)}.

  "\nNonprocessedCommonSNP:$nonprocessedCommonSNP  "
  .  eval{round($nonprocessedCommonSNP * 100 / $totalCommonSNP)}." . \% total SNPs-->  "
  . eval{round($nonprocessedCommonSNP * 100 /$totalSNP2)}.

"\n\n\nSNP(Ref2/and other)_inOnlyfile2___Means_Ref1_infile1:$SNP_inOnlyfile2_Ref1_infile1  "
  .  eval{round($SNP_inOnlyfile2_Ref1_infile1 * 100 / $totalSNP2)}." . \% total SNPs-->  "
  . eval{round($SNP_inOnlyfile2_Ref1_infile1 * 100 /$totalSNP2)}.

"\n \tSNP_inOnlyfile2_ShareRef1_Ref1_infile1:$SNP_inOnlyfile2_ShareRef1_Ref1_infile1   "
  .  eval{round($SNP_inOnlyfile2_ShareRef1_Ref1_infile1 * 100 / $totalSNP2)}
  .

"\n \tSNP_inOnlyfile2_DontShareRef1_Ref1_infile1:$SNP_inOnlyfile2_DontShareRef1_Ref1_infile1   "
  . eval{round($SNP_inOnlyfile2_DontShareRef1_Ref1_infile1 * 100 / $totalSNP2)}
  .

  "\n\nNonprocessedSNP(. \% total SNPs--> SNPs in file):$nonprocessedSNP  " .  eval{round($nonprocessedSNP * 100 / $totalSNP2)}
;


my $percentRef2 =
  ( $All1_eq_All2_and_len1_ofcourseNoRef1 +
      $A1eqA2_len2_noRef +
      $A1eqA2_len2_yesRef1 +
      $A1eqA2_len2up_noRef1 +
      $A1eqA2_len2up_yesRef1 +
      $A1neqA2_len2_yesRef1 +
      $A1neqA2_len2_sharent_noRef1 +
      $A1neqA2_len2_sharent_OneShareRef1 +
      $A1neqA2_len2plus_yesRef1 +
      $A1neqA2_len2plus_sharent_noRef1 ) * 100 / $totalCommonSNP;

print "\n\n\n\nPercent contribution of Ref2 in $ARGV[0]: $percentRef2\n";
#print "Percent contribution of Ref1 in $ARGV[0]: $percentRef2\n";
close SNP2;

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
