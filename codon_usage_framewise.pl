#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long;
our($opt_s,$opt_o,$opt_p,$place,$opt_trim_at_stop,$help);

$place=2;
$opt_p="freq";
$opt_trim_at_stop="yes";
my$options=GetOptions(
    "sequence=s"=>\$opt_s,
    "out=s"=>\$opt_o,
    "print=s"=>\$opt_p,
    "roundoff=i"=>$place,
    "trim_at_stop=s"=>\$opt_trim_at_stop,
    "help"=>\$help

);


my$usage="\n
Program to count codon usage frequency in all the three frames.
The frequency in other framesis wrt to observed frequency in frame 1.
perl script options...
Options:
    -sequence Sequence file containing CDS sequencres in fasta format
    -out    Outputfile to print the results
    -print  Freq: Total number of codons.
            Rel_freq: Frequency of codon/Total num of synonymous codons.
            Both: Freq and Rel_freq.
            Obs_exp-ratio: Observed vs expected ratio of codon.
            Log_obs_exp_ratio : Log of \"Rel_freq\"/Expected frequency.

            stop_codon_freq: Total number of codons.
            stop_codon_rel_freq: Frequency of codon/Total num of synonymous codons.
            stop_codon_obs_exp_ratio: Observed vs expected ratio of codon.
            stop_codon_both
            stop_codon_log_obs_exp_ratio: Log of \"Rel_freq\"/Expected frequency.

    -trim_at_stop yes|no


    roundoff    Roundoff the decimal values to this many places[2].
\n";
die "\n$usage\n" if $help;
die "\n$usage\n" if !$opt_s;


my $Ref_seq_hash   = ReadFasta( $opt_s);
my $Ref_codon_hash = codon_hash();

my %codon_count;
my %count_synonymous;
foreach my $seqname (keys %{$Ref_seq_hash} ) {
    $$Ref_seq_hash{$seqname}=~s/\s+//g;
    next if $$Ref_seq_hash{$seqname} !~ m/^ATG/i;
    if($opt_p !~ /stop_codon/i){
        count_codon( $$Ref_seq_hash{$seqname}, 1, \%codon_count,$opt_trim_at_stop );
        count_codon( $$Ref_seq_hash{$seqname}, 2, \%codon_count,$opt_trim_at_stop );
        count_codon( $$Ref_seq_hash{$seqname}, 3, \%codon_count,$opt_trim_at_stop );
    }
    elsif($opt_p =~ /stop_codon/i){
        count_stopcodon($$Ref_seq_hash{$seqname}, 1, \%codon_count,$opt_trim_at_stop );
    }
}

open OUT,">$opt_o" or die "\nCannot create output file\n";

my($col1,$col2,$col3);
if(lc$opt_p =~ /stop_codon/i){$col1='T'; $col2='TA'; $col3='TG';}
else{$col1=1; $col2=2; $col3=3;}

foreach($col1,$col2,$col3){$_=~s/\'//g}

if (lc$opt_p =~ /stop_codon/i){print OUT " \tcoded_amino\tExpected_frequency\tAfter_\"T\"\tAfter_\"TA\"\tAfter_\"TG\"\n" ;}
else{print OUT "\tcoded_amino\tExpected_frequency\tframe1\tframe2\tframe3\n";}
foreach ( keys %codon_count ) {
    #print "Keys: $_\n";
    $count_synonymous{ $$Ref_codon_hash{$_}{amino} } += $codon_count{$_}{1}{count};

}

foreach ( keys %codon_count ) {

    print OUT "$_\t$$Ref_codon_hash{$_}{amino}\t$$Ref_codon_hash{$_}{exp_freq}\t$codon_count{$_}{$col1}{count}\t$codon_count{$_}{$col2}{count}\t$codon_count{$_}{$col3}{count}\n" if (lc$opt_p eq "freq" || lc$opt_p eq "both"||lc$opt_p eq "stop_codon_freq"||lc$opt_p eq "stop_codon_both");
    print OUT $_. "\t"
      . $$Ref_codon_hash{$_}{amino} . "\t".$$Ref_codon_hash{$_}{exp_freq}."\t"
      . roundoff($codon_count{$_}{$col1}{count} / $count_synonymous{ $$Ref_codon_hash{$_}{amino} },2 ). "\t"
      . roundoff($codon_count{$_}{$col2}{count} / $count_synonymous{ $$Ref_codon_hash{$_}{amino} },2 ). "\t"
      . roundoff($codon_count{$_}{$col3}{count} / $count_synonymous{ $$Ref_codon_hash{$_}{amino} },2 ). "\n" if (lc$opt_p eq "rel_freq"|| lc$opt_p eq "both"||lc$opt_p eq "stop_codon_rel_freq"||lc$opt_p eq "stop_codon_both");


    print OUT $_. "\t"
      . $$Ref_codon_hash{$_}{amino} . "\t".$$Ref_codon_hash{$_}{exp_freq}."\t"
      . roundoff(($codon_count{$_}{$col1}{count} / $count_synonymous{ $$Ref_codon_hash{$_}{amino} })/$$Ref_codon_hash{$_}{exp_freq},2 ). "\t"
      . roundoff(($codon_count{$_}{$col2}{count} / $count_synonymous{ $$Ref_codon_hash{$_}{amino} })/$$Ref_codon_hash{$_}{exp_freq},2 ). "\t"
      . roundoff(($codon_count{$_}{$col3}{count} / $count_synonymous{ $$Ref_codon_hash{$_}{amino} })/$$Ref_codon_hash{$_}{exp_freq},2 ). "\n" if (lc$opt_p eq "obs_exp_ratio"||lc$opt_p eq "stop_codon_obs_exp_ratio");


    print OUT $_. "\t"
      . $$Ref_codon_hash{$_}{amino} . "\t".$$Ref_codon_hash{$_}{exp_freq}."\t"
      . log_base(roundoff(($codon_count{$_}{$col1}{count} / $count_synonymous{ $$Ref_codon_hash{$_}{amino} })/$$Ref_codon_hash{$_}{exp_freq},2 ),2). "\t"
      . log_base(roundoff(($codon_count{$_}{$col2}{count} / $count_synonymous{ $$Ref_codon_hash{$_}{amino} })/$$Ref_codon_hash{$_}{exp_freq},2 ),2). "\t"
      . log_base(roundoff(($codon_count{$_}{$col3}{count} / $count_synonymous{ $$Ref_codon_hash{$_}{amino} })/$$Ref_codon_hash{$_}{exp_freq},2 ),2). "\n" if (lc$opt_p eq "log_obs_exp_ratio"||lc$opt_p eq "stop_codon_log_obs_exp_ratio");
}

########################################################################################################################################
# subroutines

sub ReadFasta {    # to read fasta format files into hash. returns hash.

    my $seqfile = shift;

    my ( $header, @sequence );
    chomp $seqfile;
    open FASTA, "$seqfile";
    print "reading Sequences from input file.....Plz wait...\n";
    my %seq_hash = ();

    #$seq_hash{'RS_Concatenated'}="";

    $/ = "\n>";    # Change record seperator to read Fasta
    my $last_N = 1;
    while (<FASTA>) {
        chomp;
        ( $header, @sequence ) = split( "\n", $_ );

        $header =~ s/>//;       # Remove Leading > from Header
        $header =~ s/\s*$//;    # Remove trailing spaces from header
        $header =~ s/^\s*//;    # Remove Leading spaces from Header

        my $sequence = join( "", @sequence );
        $sequence =~ s/\s//g;
        $sequence =~ s/\n//g;

        if ( $header =~ /^\s*$/ ) { next; }

        # Make sure no duplicate header exists else they will be overwritten and lost during hashing due to same key.
        if ( !exists $seq_hash{$header} ) {
            $seq_hash{$header} = $sequence;    #feed headers and sequences in hash.
                                               #$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
        }
        else {

            # find a uniq header name by adding a number at the end. If header still exists, increase the number by one
            while ( exists $seq_hash{$header} ) { $header = $header . $last_N; $last_N++; }

            $seq_hash{$header} = $sequence;

            #$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;

        }
    }

    my @seq_count = keys(%seq_hash);
    my $seq_count = @seq_count;

    print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
    @seq_count = ();
    $/         = "\n";    # Record seperator set back to default newline.

    return ( \%seq_hash );

}

#-------------------------------------End ReadFasta---------------------------------------+

sub min {
    @_ = sort { $a <=> $b } @_;
    return $_[0];

}

sub max {

    @_ = sort { $a <=> $b } @_;
    return $_[-1];

}

sub roundoff{
    my$float=shift;
    my$place=shift;
    my$before_decimal=index($float, ".");
    my$temp_float=substr($float,0,$place+$before_decimal+2);
    my$add="0.".("0"x$place)."5";
    $add*=1;
    $temp_float+=$add;
    my$new_float=substr($temp_float,0,$place+$before_decimal+1);
    $new_float=~s/\.// if $place==0;
    return($new_float);
}

sub log_base {
my $n = shift;
my $base= shift;
return log($n)/log($base);
}

sub revcomp {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq;
}

sub GC_content {
    my $sequence   = shift;
    my $GCnumber   = $sequence =~ tr/GCgc/GCgc/;
    my $length_seq = $sequence =~ tr/ATGCatgc/ATGCatgc/;

    my $GCcontent = $GCnumber / $length_seq;
    return ($GCcontent);
}

sub GC_content_percent {
    my $sequence   = shift;
    my $GCnumber   = $sequence =~ tr/GCgc/GCgc/;
    my $length_seq = $sequence =~ tr/ATGCatgc/ATGCatgc/;

    my $GCcontent = $GCnumber / $length_seq;
    return ( $GCcontent * 100 );
}

sub remove_extension {
    my $filename = shift;
    my @names = split( /\./, $filename );
    $names[-1] = "";
    return join(@names);
}

sub getfilename {
    my $rawfilename = shift;

    #remove folder names from rawfilename names.
    my @names = split( /[\/\\]/, $rawfilename );
    return $names[-1];

}


######################################
sub codon_hash {
    my $codon_hash;
    my %g;
    $g{TCA}{amino} = 'S';
    $g{TCC}{amino} = 'S';
    $g{TCG}{amino} = 'S';
    $g{TCT}{amino} = 'S';
    $g{TTC}{amino} = 'F';
    $g{TTT}{amino} = 'F';
    $g{TTA}{amino} = 'L';
    $g{TTG}{amino} = 'L';
    $g{TAC}{amino} = 'Y';
    $g{TAT}{amino} = 'Y';
    $g{TAA}{amino} = '*';
    $g{TAG}{amino} = '*';
    $g{TGC}{amino} = 'C';
    $g{TGT}{amino} = 'C';
    $g{TGA}{amino} = '*';
    $g{TGG}{amino} = 'W';
    $g{CTA}{amino} = 'L';
    $g{CTC}{amino} = 'L';
    $g{CTG}{amino} = 'L';
    $g{CTT}{amino} = 'L';
    $g{CCA}{amino} = 'P';
    $g{CCC}{amino} = 'P';
    $g{CCG}{amino} = 'P';
    $g{CCT}{amino} = 'P';
    $g{CAC}{amino} = 'H';
    $g{CAT}{amino} = 'H';
    $g{CAA}{amino} = 'Q';
    $g{CAG}{amino} = 'Q';
    $g{CGA}{amino} = 'R';
    $g{CGC}{amino} = 'R';
    $g{CGG}{amino} = 'R';
    $g{CGT}{amino} = 'R';
    $g{ATA}{amino} = 'I';
    $g{ATC}{amino} = 'I';
    $g{ATT}{amino} = 'I';
    $g{ATG}{amino} = 'M';
    $g{ACA}{amino} = 'T';
    $g{ACC}{amino} = 'T';
    $g{ACG}{amino} = 'T';
    $g{ACT}{amino} = 'T';
    $g{AAC}{amino} = 'N';
    $g{AAT}{amino} = 'N';
    $g{AAA}{amino} = 'K';
    $g{AAG}{amino} = 'K';
    $g{AGC}{amino} = 'S';
    $g{AGT}{amino} = 'S';
    $g{AGA}{amino} = 'R';
    $g{AGG}{amino} = 'R';
    $g{GTA}{amino} = 'V';
    $g{GTC}{amino} = 'V';
    $g{GTG}{amino} = 'V';
    $g{GTT}{amino} = 'V';
    $g{GCA}{amino} = 'A';
    $g{GCC}{amino} = 'A';
    $g{GCG}{amino} = 'A';
    $g{GCT}{amino} = 'A';
    $g{GAC}{amino} = 'D';
    $g{GAT}{amino} = 'D';
    $g{GAA}{amino} = 'E';
    $g{GAG}{amino} = 'E';
    $g{GGA}{amino} = 'G';
    $g{GGC}{amino} = 'G';
    $g{GGG}{amino} = 'G';
    $g{GGT}{amino} = 'G';

    # expested frequence
    $g{GCT}{exp_freq} = 0.25;
    $g{GCA}{exp_freq} = 0.25;
    $g{GCC}{exp_freq} = 0.25;
    $g{GCG}{exp_freq} = 0.25;
    $g{AGA}{exp_freq} = 0.17;
    $g{AGG}{exp_freq} = 0.17;
    $g{CGA}{exp_freq} = 0.17;
    $g{CGT}{exp_freq} = 0.17;
    $g{CGG}{exp_freq} = 0.17;
    $g{CGC}{exp_freq} = 0.17;
    $g{AAT}{exp_freq} = 0.5;
    $g{AAC}{exp_freq} = 0.5;
    $g{GAT}{exp_freq} = 0.5;
    $g{GAC}{exp_freq} = 0.5;
    $g{TGT}{exp_freq} = 0.5;
    $g{TGC}{exp_freq} = 0.5;
    $g{TGA}{exp_freq} = 0.33;
    $g{TAA}{exp_freq} = 0.33;
    $g{TAG}{exp_freq} = 0.33;
    $g{CAA}{exp_freq} = 0.5;
    $g{CAG}{exp_freq} = 0.5;
    $g{GAA}{exp_freq} = 0.5;
    $g{GAG}{exp_freq} = 0.5;
    $g{GGA}{exp_freq} = 0.25;
    $g{GGT}{exp_freq} = 0.25;
    $g{GGG}{exp_freq} = 0.25;
    $g{GGC}{exp_freq} = 0.25;
    $g{CAT}{exp_freq} = 0.5;
    $g{CAC}{exp_freq} = 0.5;
    $g{ATT}{exp_freq} = 0.33;
    $g{ATC}{exp_freq} = 0.33;
    $g{ATA}{exp_freq} = 0.33;
    $g{TTG}{exp_freq} = 0.17;
    $g{CTT}{exp_freq} = 0.17;
    $g{CTG}{exp_freq} = 0.17;
    $g{CTC}{exp_freq} = 0.17;
    $g{TTA}{exp_freq} = 0.17;
    $g{CTA}{exp_freq} = 0.17;
    $g{AAG}{exp_freq} = 0.5;
    $g{AAA}{exp_freq} = 0.5;
    $g{ATG}{exp_freq} = 1;
    $g{TTT}{exp_freq} = 0.5;
    $g{TTC}{exp_freq} = 0.5;
    $g{CCA}{exp_freq} = 0.25;
    $g{CCT}{exp_freq} = 0.25;
    $g{CCC}{exp_freq} = 0.25;
    $g{CCG}{exp_freq} = 0.25;
    $g{TCT}{exp_freq} = 0.17;
    $g{TCA}{exp_freq} = 0.17;
    $g{AGT}{exp_freq} = 0.17;
    $g{TCC}{exp_freq} = 0.17;
    $g{AGC}{exp_freq} = 0.17;
    $g{TCG}{exp_freq} = 0.17;
    $g{ACA}{exp_freq} = 0.25;
    $g{ACT}{exp_freq} = 0.25;
    $g{ACC}{exp_freq} = 0.25;
    $g{ACG}{exp_freq} = 0.25;
    $g{TGG}{exp_freq} = 1;
    $g{TAT}{exp_freq} = 0.5;
    $g{TAC}{exp_freq} = 0.5;
    $g{GTT}{exp_freq} = 0.25;
    $g{GTG}{exp_freq} = 0.25;
    $g{GTC}{exp_freq} = 0.25;
    $g{GTA}{exp_freq} = 0.25;

    return ( \%g );

}

sub count_codon {
    my $sequence = shift;
    my $frame    = shift;
    my $Ref_hash = shift;
    my$trim_at_stop=shift;
    $sequence=uc$sequence;
    my$stop_codons=0;
    #print "\n\nSequence:$sequence\nLength:".length($sequence)."\n" if($stop_codons >1);
    for ( my $i = $frame-1; $i <= length($sequence); $i = $i + 3 ) {

        my $codon = substr( $sequence, $i, 3 );
        $codon =~ s/\s+//g;
        next if $codon =~ m/[^ATGC]/i;
        next if length($codon) < 3;

        if ($stop_codons ==0 && lc$trim_at_stop eq "yes"){$$Ref_hash{$codon}{$frame}{'count'}++}
        elsif(lc$trim_at_stop eq "no"){$$Ref_hash{$codon}{$frame}{'count'}++}


        $stop_codons++ if ($codon eq "TAA"||$codon eq "TAG"||$codon eq "TGA");

        #if($stop_codons >1){if ($codon eq "TAA"||$codon eq "TAG"||$codon eq "TGA"){print "*$codon"}else{print "$codon"}}

    }

    if(lc$trim_at_stop eq "yes" && $stop_codons > 1){print "\nWarning: Sequence has $stop_codons stop codon\nTruncating sequence to first stop codon\n\n"}

}


sub count_stopcodon{
my $sequence = shift;
    my $frame    = shift;
    my $Ref_hash = shift;
    my$trim_at_stop=shift;
    $sequence=uc$sequence;
    my$stop_codons=0;

    for ( my $i = $frame-1; $i <= length($sequence)-3; $i = $i + 3 ) {

       my $codon = substr( $sequence, $i, 3 );
        $codon =~ s/\s+//g;

        $stop_codons++ if ($codon eq "TAA"||$codon eq "TAG"||$codon eq "TGA");
        if ($stop_codons ==0 && lc$trim_at_stop eq "yes"){
            if(substr($codon,2,1) eq "T"){
                my$next_codon=substr($sequence,$i+3,3);
                next if $next_codon =~ m/[^ATGC]/i;
                next if length($next_codon) < 3;
                $$Ref_hash{$next_codon}{$frame}{'count'}++;
                $$Ref_hash{$next_codon}{T}{'count'}++;
             }

            if(substr($codon,1,2) eq "TA"){
                my$next_codon=substr($sequence,$i+3,3);
                next if $next_codon =~ m/[^ATGC]/i;
                next if length($next_codon) < 3;
                $$Ref_hash{$next_codon}{$frame}{'count'}++;
                $$Ref_hash{$next_codon}{TA}{'count'}++;
             }

            if(substr($codon,1,2) eq "TG"){
                my$next_codon=substr($sequence,$i+3,3);
                next if $next_codon =~ m/[^ATGC]/i;
                next if length($next_codon) < 3;
                $$Ref_hash{$next_codon}{$frame}{'count'}++;
                $$Ref_hash{$next_codon}{TG}{'count'}++;
             }

        }

        elsif(lc$trim_at_stop eq "no"){
            if(substr($codon,2,1) eq "T"){
                my$next_codon=substr($sequence,$i+3,3);
                next if $next_codon =~ m/[^ATGC]/i;
                next if length($next_codon) < 3;
                $$Ref_hash{$next_codon}{$frame}{'count'}++;
                $$Ref_hash{$next_codon}{T}{'count'}++;
             }

            if(substr($codon,1,2) eq "TA"){
                my$next_codon=substr($sequence,$i+3,3);
                next if $next_codon =~ m/[^ATGC]/i;
                next if length($next_codon) < 3;
                $$Ref_hash{$next_codon}{$frame}{'count'}++;
                $$Ref_hash{$next_codon}{TA}{'count'}++;
             }

            if(substr($codon,1,2) eq "TG"){
                my$next_codon=substr($sequence,$i+3,3);
                next if $next_codon =~ m/[^ATGC]/i;
                next if length($next_codon) < 3;
                $$Ref_hash{$next_codon}{$frame}{'count'}++;
                $$Ref_hash{$next_codon}{TG}{'count'}++;
             }

        }
    }

}