#/usr/perl/bin -w
use strict;
use Time::localtime;
use Getopt::Long;

our ($ref_seq,@sam_seq,$col_refname,$col_refpos,$col_AlleleVar,	$col_freq,$col_counts,$col_coverage,$format,$help);
$format='clcbio';
GetOptions(
		'r|reference:s' => \$ref_seq,
		's|samples:s'=>\@sam_seq,
		'f|format:s'=>\$format, 			#Format of SNP files. Program which produced SNP files. e.g. clcbio,soapdenovo,unknown
		'col_refname:i'=>\$col_refname,		#if unknown, column number conatining name of reference sequence.
		'col_refpos:i'=>\$col_refpos,		#if unknown, column number containing position of SNP on reference sequence.
		'col_AlleleVar:i'=>\$col_AlleleVar,	#if unknown, column number containing allele variations iformation.
		'col_freq:i'=>\$col_freq,			#if unknown, column number conataining allele frequencies information.
		'col_counts:i'=>\$col_counts,		#if unknown, column number containing read counts per variations.
		'col_coverage:i'=>\$col_coverage,	#if unknown, column number conatining total number of reads at a loci.
		'h|help'=>\$help);

die help() if $help;
die help() if !$ref_seq;
open LOG, ">>$0.Log.txt";
print LOG ctime(), "\n@ARGV\n";
print ctime(), "\n$0 @ARGV\n";



#my($Mapping, $ReferencePosition, $ConsensusPosition, $VariationType, $Length, $Reference, $Variants, $AlleleVariations2, $Frequencies, $Counts, $Coverage, $Variant_1, $Frequencyof_1, $Countof_1, $Variant_2, $Frequencyof_2, $Countof_2, $OverlappingAnnotations, $AminoAcidChange);
my ( %SNP, %SNP_summary);

#open SNP2, "$ARGV[1]" or die "cannot find $ARGV[1] 2nd SNP file";
my $totalSNP1 = 0;

##Create hash of SNPs found in reference sequence
($SNP{$ref_seq},$SNP_summary{$ref_seq}) = Get_SNPs($ref_seq);
if(@sam_seq>0){
	foreach my $sam_seq(@sam_seq){
		($SNP{$sam_seq},$SNP_summary{$sam_seq}) = Get_SNPs($sam_seq);
		####test print
		print "Read SNP from :$sam_seq\n";
	}
}

# compare SNPs in samples with reference.
foreach my$sam_names(@sam_seq){

	###test print
	#print "\ncurrent file:$sam_names\n";

	foreach my$loc_snp(keys %{$SNP{$sam_names}}){

		###test print
		#print "\ncurrent location:$loc_snp\n";

		#### if sample SNP also exits in reference sequence.
		if(exists $SNP{$ref_seq}{$loc_snp}){

			if($SNP{$ref_seq}{$loc_snp}{'AlleleVariations'} eq $SNP{$sam_names}{$loc_snp}{'AlleleVariations'}){
				$SNP_summary{$sam_names}{'2.Loci with NO DIFFERENCE in SNPs'}++;
			}

			elsif($SNP{$ref_seq}{$loc_snp}{'AlleleVariations'} ne $SNP{$sam_names}{$loc_snp}{'AlleleVariations'}){
				$SNP_summary{$sam_names}{'3.0.Loci with DIFFERENCE in SNPs'}++;

				if(sharent($SNP{$ref_seq}{$loc_snp}{'AlleleVariations'} ,$SNP{$sam_names}{$loc_snp}{'AlleleVariations'})>0){
					$SNP_summary{$sam_names}{"3.1.\tLoci with shared variations"}++;
				}

				if(sharent($SNP{$ref_seq}{$loc_snp}{'AlleleVariations'} ,$SNP{$sam_names}{$loc_snp}{'AlleleVariations'})==0){
					$SNP_summary{$sam_names}{"3.2.\tLoci with no common variations"}++;
				}

			}

		}

		#### if sample SNP does not exits in reference sequence.
		else{
			###test print
			#print "$loc_snp does not exists\n"
			$SNP_summary{$sam_names}{'4.Loci unique in sample sequence'}++;


		}

	}
	$SNP_summary{$sam_names}{'5.Loci unique in Reference sequence'}=$SNP_summary{$ref_seq}{'1.Total number of Loci'}-$SNP_summary{$sam_names}{'2.Loci with NO DIFFERENCE in SNPs'}-$SNP_summary{$sam_names}{'3.0.Loci with DIFFERENCE in SNPs'};
}



foreach my$samseq(keys %SNP_summary){
	print "\n\n\nBetween $ref_seq and  $samseq\n";
	foreach my$summary(sort keys %{$SNP_summary{$samseq}}){
		print "\t$summary: $SNP_summary{$samseq}{$summary} \(".round(eval($SNP_summary{$samseq}{$summary}*100/$SNP_summary{$samseq}{'1.Total number of Loci'}),2)."\%\)\n" if $summary ne '5.Loci unique in Reference sequence';
		print "\t$summary: $SNP_summary{$samseq}{$summary} \(".round(eval($SNP_summary{$samseq}{$summary}*100/$SNP_summary{$ref_seq}{'1.Total number of Loci'}),2)."\%\)\n" if $summary eq '5.Loci unique in Reference sequence';
	}
}




























####################################################################
# subroutines

sub Get_SNPs {

    my $snp_file = shift;
    open SNP1, "$snp_file";    # or die "cannot find Reference SNP file $ref_seq.";

	my %totalSNP1;
	my %SNP_NEW;
	   $totalSNP1{'1.Total number of Loci'}=0;
    while (<SNP1>) {
        $totalSNP1{'1.Total number of Loci'}++;
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
        s/^\s+//g;

		# split lines if format is clcbio.
		(
            $Mapping,           $ReferencePosition, $ConsensusPosition, $VariationType,          $Length,    $Reference,     $Variants,
            $AlleleVariations1, $Frequencies,       $Counts,            $Coverage,               $Variant_1, $Frequencyof_1, $Countof_1,
            $Variant_2,         $Frequencyof_2,     $Countof_2,         $OverlappingAnnotations, $AminoAcidChange
        ) = split( /\t/, $_ ) if $format =~/'clcbio'/gi;




        ##Test print
#print "\nblank :$blank Mapping :$Mapping ReferencePosition :$ReferencePosition ConsensusPosition :$ConsensusPosition VariationType :$VariationType Length :$Length Reference :$Reference Variants :$Variants AlleleVariations1 :$AlleleVariations1 Frequencies :$Frequencies Counts :$Counts Coverage :$Coverage Variant_1 :$Variant_1 Frequencyof_1 :$Frequencyof_1 Countof_1 :$Countof_1 Variant_2 :$Variant_2 Frequencyof_2 :$Frequencyof_2 Countof_2 :$Countof_2 OverlappingAnnotations :$OverlappingAnnotations AminoAcidChange :$AminoAcidChange\n";
#print "\n--->Mapping: $Mapping,Refposition: $ReferencePosition,ConcensusPosition: $ConsensusPosition,Variation: $VariationType, Length: $Length, Reference: $Reference, Variants:$Variants, AlleleVariations:$AlleleVariations2, Frequencies:$Frequencies, Counts:$Counts, Coverage:$Coverage,Var1: $Variant_1, Freq1:$Frequencyof_1,Count1: $Countof_1, Var2:$Variant_2, Freq2:$Frequencyof_2, Count2:$Countof_2, OverlapAnnot:$OverlappingAnnotations, AaChange:$AminoAcidChange       \n";
        $Mapping           =~ s/\s+//g;
        $ReferencePosition =~ s/\D+//g;

		#sort allele variations in ascending order for easier comparison later.frequencies and counts need to be sorte accordingly.
		($AlleleVariations1,$Frequencies,$Counts)=sort_al($AlleleVariations1, $Frequencies,$Counts) if $Variants>1;

        $SNP_NEW {"$Mapping.$ReferencePosition"} ={
                'Mapping'                => $Mapping,
                'ReferencePosition'      => $ReferencePosition,
                'ConsensusPosition'      => $ConsensusPosition,
                'VariationType'          => $VariationType,
                'Length'                 => $Length,
                'Reference'              => $Reference,
                'Variants'               => $Variants,
                'AlleleVariations'       => $AlleleVariations1,
                'Frequencies'            => $Frequencies,
                'Counts'                 => $Counts,
                'Coverage'               => $Coverage
                #'Variant_1'              => $Variant_1,
                #'Frequencyof_1'          => $Frequencyof_1,
                #'Countof_1'              => $Countof_1,
                #'Variant_2'              => $Variant_2,
                #'Frequencyof_2'          => $Frequencyof_2,
                #'Countof_2'              => $Countof_2,
                #'OverlappingAnnotations' => $OverlappingAnnotations,
                #'AminoAcidChange'        => $AminoAcidChange
            }
        ;
    }
	close SNP1;
	return ( \%SNP_NEW,\%totalSNP1);
}

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

sub round {
        my $number=shift;
		my $decimals = shift;
		$decimals=$decimals?$decimals:3;
        return(substr( $number + ( '0.' . '0' x $decimals . '5' ), 0, $decimals + length( int($number) ) + 1 ));

    }

sub sort_al{
	my($AlleleVariations, $Frequencies,$Counts)=@_;
	my(%hash_al,@AlleleVariations, @Frequencies,@Counts);


	my@alvar=split(/\//,$AlleleVariations);
	## return same thing back if they are already sorted.
	if(join "",@alvar eq join "",sort@alvar){
		return ($AlleleVariations, $Frequencies,$Counts)
	}
	else{
		my@alfreq=split(/\//,$Frequencies);
		my@alcount=split(/\//,$Counts);
		for(my$i=0;$i<@alvar;$i++){
			$hash_al{$alvar[$i]}{'freq'}=$alfreq[$i];
			$hash_al{$alvar[$i]}{'count'}=$alcount[$i];
		}
		foreach(sort @alvar){
			push(@AlleleVariations,$_);
			push(@Frequencies,$hash_al{$_}{'freq'});
			push(@Counts,$hash_al{$_}{'count'});
		}

		return(join('/',@AlleleVariations),join('/',@Frequencies),join('/',@Counts));
	}
}

sub remove_nonword{
	my $stringn=shift;
	$stringn=~s/\W//g;
	return $stringn;

}

sub help{

"This script summarizes the SNP data obtained from different SNP analysis softwares
usage: perl Script	-f input_format -r reference -s samples1 -s sample2 ......

-r|-reference	File containing SNP infomation by aligning reference reads on reference sequence.
-s|-sample		Files containing sample SNP information by aligning sample reads on reference sequence.
-f|format		Format of SNP files. Program which produced SNP files.
				clcbio,soapdenovo,unknown
-col_refname	if unknown, column number conatining name of reference sequence.
-col_refpos		if unknown, column number containing position of SNP on reference sequence.
-col_AlleleVar	if unknown, column number containing allele variations iformation.
-col_freq		if unknown, column number conataining allele frequencies information.
-col_counts		if unknown, column number containing read counts per variations.
-col_coverage	if unknown, column number conatining total number of reads at a loci.
";




}