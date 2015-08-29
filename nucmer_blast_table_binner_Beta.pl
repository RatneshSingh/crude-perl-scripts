#!/usr/bin/perl -w

# This script is created to parse blast or show-coords result (in table format) to filter
# hits based on user defined criteria. It groups the hits based on the %identity and produces
# data which can be plottted directly using excel or graphing software. This script can also extract the
# hit regions form the query and subject sequences based on the blast co-ordinates for further analysis.
# it can also pick the query and subject sequences (full seq, not extracted) for selected hits if sequence file is provided.

use strict;
use Getopt::Std;
use Getopt::Long;
our ( $opt_f, $opt_m, $opt_x, $opt_l, $opt_q, $opt_r, $opt_b, $opt_o, $opt_t, $opt_s, $opt_u, $opt_p, $opt_h, $opt_tail, $opt_head, $opt_tableout, $opt_refseq, $opt_queryseq,
    $opt_getpair, $opt_getpairaligned, $opt_extractseqaligned, $opt_d, $opt_c, $opt_extractseq, $firstBestHit );
our ( @selected, %uniq_pair );
our ( $ref_qseq, $ref_rseq, $opt_escape_spchar );

#getopt('fmxlqrbotsuph');

# initialize default values for some options

$opt_m    = 1;
$opt_x    = 100;
$opt_l    = 1;
$opt_q    = 0;         # dont make it 1. it might filter out lot of sequences if query or subject is very long and coverage is in decimal points.
$opt_r    = 0;         # dont make it 1. it might filter out lot of sequences if query or subject is very long and coverage is in decimal points.
$opt_b    = 1;
$opt_s    = 'no';
$opt_u    = 'no';
$opt_p    = 'blast';
$opt_c    = 1;
$opt_d    = " ";
$opt_head = 0;
$opt_tail = 0;

# use getopt long to accomodate too many options. Getopt::std can only take single alphabet and i am running out of appropriate flags.
my $result = GetOptions(
    "file|f=s"            => \$opt_f,    # string [required]
    "program|p=s"         => \$opt_p,    # Program name e.g blast or show-coords [blast]
    "maxp|maxpercent|x=f" => \$opt_x,    # max limit of the percent for bin [100]
    "len|length|l=i"      => \$opt_l,    # min alignment length to be included [all]

    # min % of query/subject length to be included in the alignment in order to be counted as hit.
    #works for show-coords. Use Blast_parser_BIOPERL program to parse regular blast result
    #into table for this options to work. Default table from blast does not have coverage included.
    "qcoverage|qcov|q=f"                   => \$opt_q,
    "rcoverage|referencecoverage|rcov|r=f" => \$opt_r,

    # size of each bin. [1]. eg Count all the hits 70% to 80% in one bin use bin size of 10.
    "bin|binsize|b=f" => \$opt_b,

    "out|output|o=s" => \$opt_o,         # output file to save results

    # if you want to print the resulting table along with binning result. use this options.
    #will output results which has passed the filtering criteria.
    "table|t"         => \$opt_t,
    "tout|tableout=s" => \$opt_tableout,    # save table in a file

    #  use this options to exclude self hits in case of selfblast. if query and reference(subject)
    # have same name it will not be processed and will be rejected.
    "selfhits|s=s" => \$opt_s,

    #To use only one result (with the longest alignment length) amongst multiple hits between two sequences.
    "uniqpair|u|uniqp|uniq=s" => \$opt_u,

    # Use only first best hit basedon length/percent for calculation. Ignore other hits
    "FirstBestHit|fbh=s" => \$firstBestHit,

    "minp|minpercent|m=f" => \$opt_m,       # max limit of the percent for bin [100]
    "help|h"              => \$opt_h,

    # for help file to be printed. program will not run and help file will be printed if this flag is present.
    "ref_seq|subject_seq=s" => \$opt_refseq,      # Subject sequence in fasta format
    "query_seq=s"           => \$opt_queryseq,    # query sequence file in fasta format

    # get the sequence for the query and reference pair in each hit and save into a file.
    "get_pair|gp=s" => \$opt_getpair,

    "get_pair_aligned|gpa=s"    => \$opt_getpairaligned,       # save aligned pairs in this folder
    "extract_seq|es=s"          => \$opt_extractseq,           # Extract fragments for selected based on coordinates in table.
    "extract_seq_aligned|esa=s" => \$opt_extractseqaligned,    # align extracted sequences in this folder
    "tail=i"                    => \$opt_tail,                 #add extra length in the extracted sequence on tail.
    "head=i"                    => \$opt_head,                 #add extra length in the extracted sequence on head.

    # delimiter to break the name of sequence in fasta file. Often needed when the names in fasta
    #format is longer than oe different than the one used in blast or show-coords file.
    "delimiter|d=s" => \$opt_d,

    # used in conjunction with delimiter option. Tells which column should be used as sequence name for comparision purpose.
    "column|col|c=i" => $opt_c,

    "escape_spchar|esc" => \$opt_escape_spchar                 # to escape special characters in the name if they are causing problems.

);

my $usage = "
Please run show-coords utility with options; -c -l -q   on nucmer output (.delta) before feeding to this script. or
Please run Blast_parser_BIOPERL.pl program on regular blast output to create table format.

usage: perl $0 -options....

-f,--file		Name of inputfile.
			This file will be output of blast result in table format or
			from running show-coords utility with options; -c -l -q   on nucmer output (.delta).
-p,--program		The name of program from which input is obtained. blast,show-coords[show-coords]
-m,--minp		Minimum percent to start binning[1]
-x,--maxp		maximum percent to stop binning at[100]
-l,--len		Alignment length cutoff for filtered sequences[1]
-q,--qcov		Query coverage cutoff for filtering[1]
-r,--rcov		Reference coverage cutoff for filtering[1]
-b,--bin		bin size[1]
-o,--out		outputfile. Save output in file
-t,--table		TRUE|FALSE. Print filtered table on STDOUT [FALSE]
--tout			File name. Print fileterd table in this file.
-s,--selfhits		yes/no/solo. Include self hits?
			yes: Count hits to self.e.g.  Include alignmens of 'A' to 'A'.[no]
			no : Do  not count all the hits to self.
			solo: Count the hits to self only for those showing non-self hits also [Not implemented yet].
-u,--uniqpair		yes: Use only uniq pairs for calculations. Removes reverse pairs from counting
						e.g. Between alignments of 'A' to 'B' and 'B' to 'A' only one alignment
						( with the largest Aln length) is included in the calculations.
						Also, Multiple hits between 'A' and 'B' will be reduced to one hit.
					no: Include all the pairs in the calculation. This is default [no].
-FirstBestHit|fbh	PerId|AlnLen|Evalue|BitScore.
					Only use first best hit for each query for calculations. First best hit will be decided on user defined parameter.
					PerId: Hit with the highest percent identity is best hit .
					AlnLen: Hit with the  largest alignment length is best hit.
					Evalue: Hit with the lowest evalue is best hit.
					BitScore: Hit with the largest bitscore is the best hit.


Following options are in beta mode:
--ref_seq		Sequence used as reference in nucmer.
--subject_seq		Sequence used as database for blast.
--query_seq		Sequence used as query file in nucmer or blast.
--get_pair|gp		Save query and subject sequence pair for each selected hit.
--get_pair_aligned|gpa	Save the alignment of the sequence pair in this folder
--extract_seq|es		query|subject|both. Extract hit regions using blast-coordinates for each selected hit.
--extract_seq_aligned|esa	Align extracted sequences and save in this folder.
--tail			Extract and save query and reference regions using blast-coordinates for each selected hit.
--head			Extract and save query and reference regions using blast-coordinates for each selected hit.
-h,--help		Print help and exit.
-escape_spchar|esc	Convert all the special character to '_' if they are causing problems
";

if ($opt_h) { die $usage }

# check the values of $firstBestHit
if(defined $firstBestHit){
	if (lc$firstBestHit eq 'perid' || lc$firstBestHit eq 'alnlen' || lc$firstBestHit eq 'evalue' || lc$firstBestHit eq 'bitscore'){}
	else{die "Please use appropriate value for option -FirstBestHit|fbh.\n Valid options are: PerId, AlnLen,Evalue,BitScore\n\n"}
}



#$opt_p='blast' if !defined $opt_p;
if   ( $opt_p eq 'show-coords' || $opt_p eq 'blast' || $opt_p eq 'blastn' ) { }
else                                                                        { die "Please provide proper value for -p\nAvailable options include: show-coords,blast, blastn\n" }

#$opt_f="test_blastFile.blastn";
# Acquire input to parse
if    ($opt_f)     { open INFILE, $opt_f   or die "\nCannot find file $opt_f\n$!" }
elsif ( $ARGV[0] ) { open INFILE, $ARGV[0] or die "Cannot find file $ARGV[0]\n\$!" }
else               { die "Provide a input file (table format from blast or show-coords) to parse\n\n$usage\n"; }

# get filenames from path. slashes in  path mess up the file names for creating file saving results
my $qfilename = getfilename($opt_queryseq) if $opt_queryseq;
my $rfilename = getfilename($opt_refseq)   if $opt_refseq;

# Read sequence file based on what need to done (query,subject,both) and what filenames are provided
if ( lc $opt_extractseq eq 'query' ) {
    undef $opt_extractseq if !$opt_queryseq;
    print "\n Warning: The query sequences will not be extracted as as query sequence file is not provided.\n"
      if !$opt_queryseq;
    $ref_qseq = ReadFasta($opt_queryseq) if $opt_queryseq;
    open EXTRACTED, ">EXTRACTED.$qfilename" or die "Unable to create output file for extracted sequence\n$!";
}

elsif ( lc $opt_extractseq eq 'subject' || lc $opt_extractseq eq 'reference' ) {
    undef $opt_extractseq if !$opt_refseq;
    print "\n Warning: The Reference sequences will not be extracted as as reference sequence file is not provided.\n"
      if !$opt_refseq;
    $ref_rseq = ReadFasta($opt_refseq) if $opt_refseq;
    open EXTRACTED, ">EXTRACTED.$rfilename" or die "Unable to create output file for extracted sequence\n$!";
}

if ( lc $opt_extractseq eq 'both' || $opt_extractseqaligned ) {
    undef $opt_extractseq if ( !$opt_refseq || !$opt_queryseq );
    undef $opt_extractseqaligned
      if ( !$opt_refseq || !$opt_queryseq );    # if no sequence files are provided, how sequences can be extracted and aligned.
    print "\n Warning: The Reference and query sequences will not be extracted as as one or both of the reference or query sequences files are not provided.\n"
      if ( !$opt_refseq || !$opt_queryseq );
    $ref_qseq = ReadFasta($opt_queryseq) if ( $opt_refseq || $opt_queryseq );
    $ref_rseq = ReadFasta($opt_refseq)   if ( $opt_refseq || $opt_queryseq );
    open EXTRACTED, ">EXTRACTED.$qfilename._.$rfilename"
      or die "Unable to create output file for extracted sequence\n$!"
      if $opt_extractseq;
    mkdir $opt_extractseqaligned if $opt_extractseqaligned;
}

# to extract pairs, create file to save sequeces
if ( $opt_getpair || $opt_getpairaligned ) {
    undef $opt_getpair        if ( !$opt_refseq || !$opt_queryseq );
    undef $opt_getpairaligned if ( !$opt_refseq || !$opt_queryseq );
    print "\n Warning: The Reference and query sequences will not be extracted as as one or both of the reference or query sequences files are not provided.\n"
      if ( !$opt_refseq || !$opt_queryseq );
    $ref_qseq = ReadFasta($opt_queryseq) if ( $opt_refseq || $opt_queryseq );
    $ref_rseq = ReadFasta($opt_refseq)   if ( $opt_refseq || $opt_queryseq );
}

# open output file to save results with the name user provided or else create automatically.
my $tablename = getfilename($opt_f);
if   ( !$opt_o ) { open OUT, ">Alignment_binned_.$tablename"; }
else             { open OUT, ">$opt_o"; }

# print program and option information inthe output file for records.
print OUT
  "These results are obtained from following command:\n$0   -file $opt_f    -minp $opt_m   -maxp $opt_x    -len $opt_l    -qcov $opt_q    -rcov $opt_r    -bin $opt_b    -selfhits $opt_s    -uniqpair $opt_u    -p $opt_p\n\n";

# print the headers
my ( %freq, %length );
for ( my $i = $opt_m + $opt_b; sprintf( "%0.2f", $i ) <= $opt_x; $i = $i + $opt_b ) {

    print OUT "\t", sprintf( "%0.2f", $i );
    print "\t", sprintf( "%0.2f", $i );

}

print OUT"\n";
print "\n";

# initialize bins and set value to zero;
my $total_length_compared = 0;
for ( my $i = $opt_m + $opt_b; sprintf( "%0.2f", $i ) <= $opt_x; $i = $i + $opt_b ) {
    $freq{ sprintf( "%0.2f", $i ) } = 0;
    $length{ sprintf( "%0.2f", $i ) } = 0;
}
$. = 0;

#skip initial 5 lines and then parse
if ( $opt_p eq 'show-coords' ) { do <INFILE> while $. <= 5 }

while (<INFILE>) {

    # Test if the file type is actually the file type user provided.
    #my@line_check=split(/\s+/,$_);
    # next if empty line.
    if ( $_ =~ /^\s*$/ ) { next }
    my $line = $_;
    my ( $LENR, $LENQ, $NAMEQ, $NAMER, $PER_IDY, $ALNLEN, $LEN1, $LEN2, $MISMATCH, $GAP, $S1, $E1, $S2, $E2, $EVALUE, $BITSCORE, $COVR, $COVQ );

    # parse lines based on what program they come from. e.g. show-coords, blast
    if ( lc $opt_p eq 'show-coords' ) {
        ( $S1, $E1, $S2, $E2, $LEN1, $LEN2, $PER_IDY, $LENR, $LENQ, $COVR, $COVQ, $NAMER, $NAMEQ ) = parse_showcoords($_);
    }
    elsif ( lc $opt_p eq 'blastn' || lc $opt_p eq 'blast' ) {
        ( $NAMEQ, $NAMER, $PER_IDY, $ALNLEN, $LEN1, $LEN2, $MISMATCH, $GAP, $S1, $E1, $S2, $E2, $EVALUE, $BITSCORE, $COVQ, $COVR, $LENQ, $LENR ) = parse_blast($_);
    }

    # next if filtering options are not met
    if ( !( ( $LEN1 >= $opt_l || $LEN2 >= $opt_l ) && $COVR >= $opt_r && $COVQ >= $opt_q && $PER_IDY >= $opt_m && $PER_IDY <= $opt_x ) ) {
        next;
    }

    # Dont proceed if self hit are not allowed and query and subject are same
    if ( $opt_s eq 'no' && $NAMER eq $NAMEQ ) { next }

    # if no uniq pairs and no self hits needed
    if ( lc $opt_u eq 'no' ) {
        push( @selected, $_ );    # if lc$opt_t eq 'yes';
        for ( my $i = $opt_m + $opt_b; $i <= $opt_x; $i = $i + $opt_b ) {

            if (   $PER_IDY > sprintf( "%0.2f", $i ) - $opt_b
                && $PER_IDY <= sprintf( "%0.2f", $i ) )
            {
                $freq{ sprintf( "%0.2f", $i ) }++;
                $length{ sprintf( "%0.2f", $i ) } += $LEN1;
                $total_length_compared += $LEN1;
            }                     #print "$PER_IDY is included in group". sprintf("%0.2f",$i)."\n"}
        }
        if ( $opt_extractseq && ( lc $opt_extractseq eq 'query' || lc $opt_extractseq eq 'both' ) ) {
            my ( $header, $sequence ) = extract_seq( $NAMEQ, $S1, $E1, $ref_qseq, $opt_tail, $opt_head );
            print EXTRACTED">$header\n$sequence\n";
        }
        if ( $opt_extractseq && ( lc $opt_extractseq eq 'subject' || lc $opt_extractseq eq 'both' ) ) {
            my ( $header1, $sequence1 ) = extract_seq( $NAMER, $S2, $E2, $ref_rseq, $opt_tail, $opt_head );
            print EXTRACTED">$header1\n$sequence1\n";
        }
        elsif ($opt_extractseqaligned) {

            my ( $header,  $sequence )  = extract_seq( $NAMEQ, $S1, $E1, $ref_qseq, $opt_tail, $opt_head );
            my ( $header1, $sequence1 ) = extract_seq( $NAMER, $S2, $E2, $ref_rseq, $opt_tail, $opt_head );

            #$header=~s/\s+/_/g;
            $header = special_char_to_underscore($header) if $opt_escape_spchar;

            #$header1=~s/\s+/_/g;
            $header1 = special_char_to_underscore($header1) if $opt_escape_spchar;

            open EXTRACTED, ">$opt_extractseqaligned/$header.$header1.fasta";
            print EXTRACTED">$header\n$sequence\n";
            print EXTRACTED">$header1\n$sequence1\n";

            system("muscle -in $opt_extractseqaligned/$header.$header1.fasta -out $opt_extractseqaligned/$header.$header1.muscle.aln");
            system("rm $opt_extractseqaligned/$header.$header1.fasta");

        }

    }

    # if uniq pairs needed, the Extraction of sequences will be done later as best one need to be calculated among all.
    elsif ( lc $opt_u eq 'yes' && !defined $firstBestHit ) {    # collect lines for uniq hits and process later

        # Arange the namer and nameq in alphabetical order to create the uniq name.
        my $uniq_name = $NAMER lt $NAMEQ ? ( $NAMER . $NAMEQ ) : ( $NAMEQ . $NAMER );

        # if this is first time, create hash key and value pair for aln_length, per_idy and whole line.
        if ( !exists $uniq_pair{'aln_length'}{$uniq_name} ) {
            $uniq_pair{'aln_length'}{$uniq_name} = ( $LEN1 + $LEN2 ) / 2;    # take average of query and subject lengths in the alignment
            $uniq_pair{'per_idy'}{$uniq_name}    = $PER_IDY;
            $uniq_pair{'whole_line'}{$uniq_name} = $_;
        }

        # if already exists, compare the new value to the existing values and use longest one.
        elsif ( exists $uniq_pair{'aln_length'}{$uniq_name} ) {

            # Select the Alignment length and Per_IDY fo rthe longest alignment if more than one alignment for same pair are detected.
            $uniq_pair{'aln_length'}{$uniq_name} = ( $LEN1 + $LEN2 ) / 2
              if $uniq_pair{'aln_length'}{$uniq_name} < ( $LEN1 + $LEN2 ) / 2;
            $uniq_pair{'per_idy'}{$uniq_name} = $PER_IDY
              if $uniq_pair{'aln_length'}{$uniq_name} < ( $LEN1 + $LEN2 ) / 2;
            $uniq_pair{'whole_line'}{$uniq_name} = $_ if $uniq_pair{'aln_length'}{$uniq_name} < ( $LEN1 + $LEN2 ) / 2;
        }

    }

    # elsif only FirstBestHit is requested select only one best hit per query and ignore others.
    elsif ( lc $opt_u eq 'yes' && $firstBestHit ) {    # collect lines for uniq hits and process later

        # Arange the namer and nameq in alphabetical order to create the uniq name.
        my $uniq_name = $NAMEQ;

        # if this is first time, create hash key and value pair for aln_length, per_idy and whole line.
        if ( !exists $uniq_pair{'aln_length'}{$uniq_name} ) {
            $uniq_pair{'aln_length'}{$uniq_name} = ( $LEN1 + $LEN2 ) / 2;    # take average of query and subject lengths in the alignment
            $uniq_pair{'per_idy'}{$uniq_name}    = $PER_IDY;
            $uniq_pair{'evalue'}{$uniq_name}     = $EVALUE;
            $uniq_pair{'bitscore'}{$uniq_name}   = $BITSCORE;
            $uniq_pair{'whole_line'}{$uniq_name} = $_;
        }

        # if already exists, compare the new value to the existing values and use longest one.
        elsif ( exists $uniq_pair{'aln_length'}{$uniq_name} ) {
            my $update_values = 'false';

            if ( lc $firstBestHit eq 'perid'    && $uniq_pair{'per_idy'}{$uniq_name} < $PER_IDY )   { $update_values = 'true' }
            if ( lc $firstBestHit eq 'evalue'   && $uniq_pair{'evalue'}{$uniq_name} > $EVALUE )     { $update_values = 'true' }
            if ( lc $firstBestHit eq 'bitscore' && $uniq_pair{'bitscore'}{$uniq_name} < $BITSCORE ) { $update_values = 'true' }
            if ( lc $firstBestHit eq 'alnlen' && $uniq_pair{'aln_length'}{$uniq_name} < ( $LEN1 + $LEN2 ) / 2 ) { $update_values = 'true' }

            if ( $update_values eq 'true' ) {

                # Select the Alignment length and Per_IDY for the longest alignment if more than one alignment for same pair are detected.
                $uniq_pair{'aln_length'}{$uniq_name} = ( $LEN1 + $LEN2 ) / 2;    # if $uniq_pair{'aln_length'}{$uniq_name} < ( $LEN1 + $LEN2 ) / 2;
                $uniq_pair{'per_idy'}{$uniq_name}    = $PER_IDY;                 # if $uniq_pair{'aln_length'}{$uniq_name} < ( $LEN1 + $LEN2 ) / 2;
                $uniq_pair{'whole_line'}{$uniq_name} = $_;                       # if $uniq_pair{'aln_length'}{$uniq_name} < ( $LEN1 + $LEN2 ) / 2;
                $uniq_pair{'evalue'}{$uniq_name}     = $EVALUE;                  # if $uniq_pair{'aln_length'}{$uniq_name} < ( $LEN1 + $LEN2 ) / 2;
                $uniq_pair{'bitscore'}{$uniq_name}   = $BITSCORE;                # if $uniq_pair{'aln_length'}{$uniq_name} < ( $LEN1 + $LEN2 ) / 2;
            }
        }

    }

    else { print "This one passed filteration but is not processed:$_\n" }

}

# calculate frequencies of uniq pairs if $opt_u. will have to again split the selected lines and retrieve all the info again as we are out of the read file loop..
if ( lc $opt_u eq 'yes' ) {

    foreach my $blastline ( keys %{ $uniq_pair{'aln_length'} } ) {

        my ( $NAMEQ, $NAMER, $PER_IDY, $ALNLEN, $LEN1, $LEN2, $MISMATCH, $GAP, $S1, $E1, $S2, $E2, $EVALUE, $BITSCORE, $COVQ, $COVR, $LENQ, $LENR ) =
          parse_blast( $uniq_pair{'whole_line'}{$blastline} )
          if ( lc $opt_p eq 'blastn' || lc $opt_p eq 'blast' );
        ( $S1, $E1, $S2, $E2, $LEN1, $LEN2, $PER_IDY, $LENR, $LENQ, $COVR, $COVQ, $NAMER, $NAMEQ ) = parse_showcoords( $uniq_pair{'whole_line'}{$blastline} )
          if ( lc $opt_p eq 'show-coords' );

        #print "\nBlastline:$uniq_pair{'whole_line'}{$blastline}\nALNLEN:$ALNLEN\tLEN1:$LEN1\tLEN2:$LEN2"; # print for debugging.

        # calculate frequency and total alignment length for each bin.
        for ( my $i = $opt_m + $opt_b; $i <= $opt_x; $i = $i + $opt_b ) {

            if (   $uniq_pair{'per_idy'}{$blastline} > ( sprintf( "%0.2f", $i ) - $opt_b )
                && $uniq_pair{'per_idy'}{$blastline} <= sprintf( "%0.2f", $i ) )
            {
                $freq{ sprintf( "%0.2f", $i ) }++;
                $length{ sprintf( "%0.2f", $i ) } += $LEN1;
                $total_length_compared += $LEN1;
            }
        }

        # extract and write sequences in output file if write $opt_extractseq is true
        if ( $opt_extractseq && ( lc $opt_extractseq eq 'query' || lc $opt_extractseq eq 'both' ) ) {
            my ( $header, $sequence ) = extract_seq( $NAMEQ, $S1, $E1, $ref_qseq, $opt_tail, $opt_head );
            print EXTRACTED">$header\n$sequence\n";
        }
        if ( $opt_extractseq && ( lc $opt_extractseq eq 'subject' || lc $opt_extractseq eq 'both' ) ) {
            my ( $header1, $sequence1 ) = extract_seq( $NAMER, $S2, $E2, $ref_rseq, $opt_tail, $opt_head );
            print EXTRACTED">$header1\n$sequence1\n";
        }

        # aligne extracted sequence and write $opt_extractseqaligned
        if ($opt_extractseqaligned) {

            my ( $header,  $sequence )  = extract_seq( $NAMEQ, $S1, $E1, $ref_qseq, $opt_tail, $opt_head );
            my ( $header1, $sequence1 ) = extract_seq( $NAMER, $S2, $E2, $ref_rseq, $opt_tail, $opt_head );

            #$header=~s/\s+/_/g;
            $header = special_char_to_underscore($header) if $opt_escape_spchar;

            #$header1=~s/\s+/_/g;
            $header1 = special_char_to_underscore($header1) if $opt_escape_spchar;

            open EXTRACTED, ">$opt_extractseqaligned/$header.$header1.fasta";
            print EXTRACTED">$header\n$sequence\n";
            print EXTRACTED">$header1\n$sequence1\n";

            system("muscle -in $opt_extractseqaligned/$header.$header1.fasta -out $opt_extractseqaligned/$header.$header1.muscle.aln");
            system("rm $opt_extractseqaligned/$header.$header1.fasta");

        }

        # store selected blastline to print if user asks for table.
        push( @selected, $uniq_pair{'whole_line'}{$blastline} )

    }
}

# If No data passed through the filtering criteria, notify the user.
if ( scalar @selected <= 0 ) {
    die "No lines in the provided table could pass thorugh filter. Check your parameters or debug the script.\n\n";
}

print OUT"\%of total fragments";
print "\%of total fragments";
for ( my $i = $opt_m + $opt_b; sprintf( "%0.2f", $i ) <= $opt_x; $i = $i + $opt_b ) {

    #foreach percent bin print frequency of sequences (number)
    print OUT"\t" . sprintf( "%0.2f", $freq{ sprintf( "%0.2f", $i ) } * 100 / scalar @selected )
      if defined $freq{ sprintf( "%0.2f", $i ) };
    print "\t" . sprintf( "%0.2f", $freq{ sprintf( "%0.2f", $i ) } * 100 / scalar @selected )
      if defined $freq{ sprintf( "%0.2f", $i ) };
}
print OUT"\n\%of TotalLength";
print "\n\%of TotalLength";
for ( my $i = $opt_m + $opt_b; sprintf( "%0.2f", $i ) <= $opt_x; $i = $i + $opt_b ) {

    #foreach percent bin print frequency of sequences (total length)
    print OUT"\t" . sprintf( "%0.2f", $length{ sprintf( "%0.2f", $i ) } * 100 / $total_length_compared )
      if defined $length{ sprintf( "%0.2f", $i ) };
    print "\t" . sprintf( "%0.2f", $length{ sprintf( "%0.2f", $i ) } * 100 / $total_length_compared )
      if defined $length{ sprintf( "%0.2f", $i ) };
}

print "\n\n\n*********************************************************************************************************************
		Filtered table
*********************************************************************************************************************\n"
  if $opt_t;
print "NAMEQ\tNAMES\tPER_IDY\tALN_LEN\tMISMATCH\tGAP\tQStart\tQEnd\tSStart\tSEnd\tEVALUE\tBITSCORE\tCOVQ\tCOVR\tLENQ\tLENR\n"
  if $opt_t;
if ($opt_t) {
    foreach (@selected) { print "$_\n" }
}

# print table in fle if asked by the user ($opt_tableout)
open TABLE, ">$opt_tableout" if $opt_tableout;
if ($opt_tableout) {
    print TABLE "NAMEQ\tNAMES\tPER_IDY\tALN_LEN\tMISMATCH\tGAP\tQStart\tQEnd\tSStart\tSEnd\tEVALUE\tBITSCORE\tCOVQ\tCOVR\tLENQ\tLENR\n";
    foreach (@selected) { print TABLE "$_\n" }
}

###################################################################################################################################################
# print sequence pairs if asked by the user (getpair)
if ( $opt_getpair || $opt_getpairaligned ) {

    open PAIRSEQ, ">$opt_getpair"
      if $opt_getpair;    # if $opt_getpairaligned then open file later when directory is created
    mkdir $opt_getpairaligned if $opt_getpairaligned;
    my @pairs_not_aligned;
    my @delete_list;

    foreach (@selected) {
        my @blasttemp = split( /\s+/, $_ );
        foreach (@blasttemp) { s/^\s+//g; s/\s+$//g; }

        $blasttemp[0] = special_char_to_underscore( $blasttemp[0] ) if $opt_escape_spchar;    # escape special characters
        $blasttemp[1] = special_char_to_underscore( $blasttemp[1] ) if $opt_escape_spchar;    # escape special characters

        my $temp_query_seq   = $$ref_qseq{ $blasttemp[0] };
        my $temp_subject_seq = $$ref_rseq{ $blasttemp[1] };
        my $tempseq_name     = $blasttemp[1];
        if ( $blasttemp[8] > $blasttemp[9] ) {
            $temp_subject_seq = revcomp( $$ref_rseq{ $blasttemp[1] } );
            $tempseq_name     = $blasttemp[1] . " RevComp";
        }

        open PAIRSEQ, ">$opt_getpairaligned/$blasttemp[0].$blasttemp[1].fasta" if $opt_getpairaligned;
        my $both_pair = 'true';
        if ($temp_query_seq) { print PAIRSEQ ">$blasttemp[0]\n$temp_query_seq\n"; }

        #elsif($$ref_qseq{ $blasttemp[0] }){print PAIRSEQ ">$blasttemp[0]\n$$ref_qseq{ $blasttemp[0] }\n";}
        else {
            undef $both_pair;
            push( @pairs_not_aligned, $blasttemp[0] . " in pair " . $blasttemp[0] . "_" . $blasttemp[1] . " not found" );
            print "Following_sequence is not found:\n>$blasttemp[0]\n$$ref_qseq{ $blasttemp[0] }\nSeqAgain:$temp_query_seq\n";

        }
        if ($temp_subject_seq) { print PAIRSEQ ">$tempseq_name\n$temp_subject_seq\n"; }

        #elsif($$ref_rseq{ $blasttemp[1] }){print PAIRSEQ ">$tempseq_name\n",revcomp($$ref_rseq{ $blasttemp[1] }),"\n";}
        else {
            undef $both_pair;
            push( @pairs_not_aligned, $blasttemp[1] . " in pair " . $blasttemp[0] . "_" . $blasttemp[1] . " not found" );
            print "Following_sequence is not found:\n>$blasttemp[1]\n$$ref_rseq{$blasttemp[1]}\nSeqAgain:$temp_subject_seq\n";
        }

        close PAIRSEQ
          if $opt_getpairaligned;    # only close when aligned pair is needed, else keep it open to save next sequence pair.

        # align the pairs if requested to do so
        if ( $opt_getpairaligned && defined $both_pair ) {
            system("muscle -in $opt_getpairaligned/$blasttemp[0].$blasttemp[1].fasta -out $opt_getpairaligned/$blasttemp[0].$blasttemp[1].muscle.aln") == 0
              or print "***\nCould not run Muscle on seqs.\n";
            system("rm .\/$opt_getpairaligned/$blasttemp[0].$blasttemp[1].fasta") == 0
              or push( @delete_list, "$opt_getpairaligned/$blasttemp[0].$blasttemp[1].fasta" );
        }
        else { print "\nOne of the sequence in pair not found, So not aligning the sequences\n" }

    }
    print "\n\nSaved each pair of the selected blast hits in file:$opt_getpair\n\nThe number of pairs not found are.\n";
    print scalar @pairs_not_aligned, "\n";

    print join( "\n", @pairs_not_aligned );
    if ( scalar @delete_list > 0 ) {
        foreach (@delete_list) { system( "rm", $_ ) }
    }

}
############################################################################################################################################################

print "\n\n Total number of hits passed the filter criteria:" . scalar @selected . "\n";
print OUT "\n\n Total number of hits passed the filter criteria:" . scalar @selected . "\n";
print "\n Total Length of sequences passed the filter criteria:"
  . sprintf( "%0.2f", $total_length_compared / 1000 ) . " kb/"
  . sprintf( "%0.2f", $total_length_compared / 1000000 )
  . "Mb\n\n";
print OUT"\n Total Length of sequences passed the filter criteria:"
  . sprintf( "%0.2f", $total_length_compared / 1000 ) . " kb/"
  . sprintf( "%0.2f", $total_length_compared / 1000000 )
  . "Mb\n\n";

close TABLE;
close EXTRACTED;
close OUT;

###############################################################
#                          Subroutines                        #
###############################################################

sub ReadFasta {    # to read fasta format files into hash. returns hash.

    my $seqfile = shift;
    my $demo_header;

    my ( $header, @sequence );
    chomp $seqfile;
    open FASTA, "$seqfile";
    print "\n\n*****************************\nReading Sequences from input file $seqfile.....Plz wait...\n";

    #my%seq_hash=();
    my $seq_hash = {};    # create anonymous hash and put referece of the created hash in $seq_hash
                          #$seq_hash{'RS_Concatenated'}="";

    $/ = "\n>";           # Change record seperator to read Fasta
    my $last_N = 1;
    while (<FASTA>) {
        chomp;
        ( $header, @sequence ) = split( "\n", $_ );

        $header =~ s/>//;       # Remove Leading > from Header
        $header =~ s/\s*$//;    # Remove trailing spaces from header
        $header =~ s/^\s*//;    # Remove Leading spaces from Header
        my @scaffold_name = split( /$opt_d/, $header ) if ( $opt_c && $opt_d );
        $header = $scaffold_name[ $opt_c - 1 ] if ( $opt_c && $opt_d );
        my $sequence = join( "", @sequence );
        $sequence =~ s/\s//g;
        $sequence =~ s/\n//g;

        if ( $header =~ /^\s*$/ ) { next; }
        $demo_header = $header;

        # Make sure no duplicate header exists else they will be overwritten and lost during hashing due to same key.
        if ( !exists $$seq_hash{$header} ) {
            $header = special_char_to_underscore($header) if $opt_escape_spchar;
            $$seq_hash{$header} = $sequence;    #feed headers and sequences in hash.
                                                #$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
        }
        else {

            # find a uniq header name by adding a number at the end. If header still exists, increase the number by one
            while ( exists $$seq_hash{$header} ) { $header = $header . $last_N; $last_N++; }

            $header = special_char_to_underscore($header) if $opt_escape_spchar;
            $$seq_hash{$header} = $sequence;

            #$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;

        }
    }

    my @seq_count = keys(%$seq_hash);
    my $seq_count = @seq_count;

    print "Done....\nNumber of sequences read form $seqfile = $seq_count\n\nExample of one the headers in sequence\n$demo_header\n\n\n*****************************\n";
    @seq_count = ();
    $/         = "\n";    # Record seperator set back to default newline.

    return ($seq_hash);

}

#-------------------------------------End ReadFasta---------------------------------------+

sub extract_seq {
    my $query1   = shift;
    my $start_s1 = shift;
    my $end_s1   = shift;
    my $refseq   = shift;
    my $tail     = shift;
    my $head     = shift;
    my $len1     = $end_s1 - $start_s1 + 1 + $tail if ( $end_s1 > $start_s1 );
    my $strand   = 'plus';

    $query1 = special_char_to_underscore($query1) if $opt_escape_spchar;

    #$query1=special_char_to_underscore($query1) if $opt_escape_spchar;
    # if sstart is bigger than end
    $len1     = $start_s1 - $end_s1 + 1 if ( $end_s1 < $start_s1 );    #calculate length to extract
    $strand   = 'minus'                 if ( $end_s1 < $start_s1 );    # change the strand to 'minus'
    $start_s1 = $end_s1                 if ( $end_s1 < $start_s1 );    #make make end as start for extraction.

    # extract sequence
    my $new_sequence1 = substr( $$refseq{$query1}, $start_s1 - 1 - $head, $len1 ) if defined $$refseq{$query1};

    if ( $strand eq 'minus' ) {
        my $new_sequence2 = reverse $new_sequence1;
        $new_sequence2 =~ tr/atgcATGC/tacgTACG/;
        $new_sequence1 = uc $new_sequence2;

        #print "seq reversed\n\n";
    }

    my $returnedlength = length($new_sequence1);

    my $new_header1 = $query1 . '_' . 'start-' . $start_s1 . '_' . 'end-' . $end_s1;

    return ( $new_header1, $new_sequence1, $returnedlength ) if defined $new_sequence1;
    print "Cannot find $$query1\n" if !defined $$refseq{$query1};
}

sub getfilename {
    my $rawfilename = shift;

    #remove folder names from rawfilename names.
    my @names = split( /[\/\\]/, $rawfilename );
    return $names[-1];

}

sub parse_showcoords {
    my $line = shift;
    my ( $LENR, $LENQ, $NAMEQ, $NAMER, $PER_IDY, $LEN1, $LEN2, $MISMATCH, $GAP, $S1, $E1, $S2, $E2, $EVALUE, $BITSCORE, $COVR, $COVQ );
    if ( lc $opt_p eq 'show-coords' ) {

        # split line and catch values.
        # header is : [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS
        my ( $S1_E1, $S2_E2, $LEN1_LEN2, $PER_IDY_temp, $LENR_LENQ, $COVR_COVQ, $NAMER_NAMEQ ) = split( /\|/, $line );
        chomp( $S1_E1, $S2_E2, $LEN1_LEN2, $PER_IDY_temp, $LENR_LENQ, $COVR_COVQ, $NAMER_NAMEQ );
        $PER_IDY = $PER_IDY_temp;

        foreach ( $S1_E1, $S2_E2, $LEN1_LEN2, $PER_IDY, $LENR_LENQ, $COVR_COVQ, $NAMER_NAMEQ ) {
            s/^\s+//g;
        }    #Remove leading spaces at the start of each variable to assist in proper splitting usng space as delimiter

        ( $S1,   $E1 )   = split( /\s+/, $S1_E1 );
        ( $S2,   $E2 )   = split( /\s+/, $S2_E2 );
        ( $LEN1, $LEN2 ) = split( /\s+/, $LEN1_LEN2 );

        #($PER,$IDY)=split(/\s+/,$PER_IDY);
        ( $LENR,  $LENQ )  = split( /\s+/, $LENR_LENQ );
        ( $COVR,  $COVQ )  = split( /\s+/, $COVR_COVQ );
        ( $NAMER, $NAMEQ ) = split( /\s+/, $NAMER_NAMEQ );

        foreach ( $S1, $E1, $S2, $E2, $LEN1, $LEN2, $PER_IDY, $LENR, $LENQ, $COVR, $COVQ, $NAMER, $NAMEQ ) {
            s/\s+//g;
        }    # Remove all the white spaces from name and values.
        return ( $S1, $E1, $S2, $E2, $LEN1, $LEN2, $PER_IDY, $LENR, $LENQ, $COVR, $COVQ, $NAMER, $NAMEQ );

        #print "$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ\n";
    }
}

sub parse_blast {
    my $blastline = shift;
    my ( $LENR, $LENQ, $NAMEQ, $NAMER, $PER_IDY, $ALNLEN, $LEN1, $LEN2, $MISMATCH, $GAP, $S1, $E1, $S2, $E2, $EVALUE, $BITSCORE, $COVR, $COVQ ) = 1;

 # split line and catch values.
 # header is : $queryId, $subjectId, $percIdentity, $alnLength, $mismatchCount, $gapOpenCount, $queryStart, $queryEnd, $subjectStart, $subjectEnd, $eVal, $bitScore, $Query_coverage
    ( $NAMEQ, $NAMER, $PER_IDY, $ALNLEN, $MISMATCH, $GAP, $S1, $E1, $S2, $E2, $EVALUE, $BITSCORE, $COVQ, $COVR, $LENQ, $LENR ) = split( /\s+/, $blastline );
    chomp( $NAMEQ, $NAMER, $PER_IDY, $ALNLEN, $MISMATCH, $GAP, $S1, $E1, $S2, $E2, $EVALUE, $BITSCORE, $COVQ, $COVR, $LENQ, $LENR );

    #accurate query lengths and subject lengths in alignment can be calculated
    $LEN2 = abs( $E2 - $S2 ) + 1;
    $LEN1 = abs( $E1 - $S1 ) + 1;
    if ( $COVQ =~ /\d+/ ) { $COVQ = sprintf( "%.2f", $COVQ ) }
    else                  { $COVQ = 1 }    # if coverage is not provided make it =1;
    if ( $COVR =~ /\d+/ ) { $COVR = sprintf( "%.2f", $COVR ) }
    else                  { $COVR = 1 }    # if coverage is not provided make it =1;

    $LENQ = 'NA' if !$LENQ;
    $LENR = 'NA' if !$LENR;

    #$LENR=$LENQ='NA';
    #Remove trailing spaces at the start of each variable to assist in proper splitting usng space as delimiter
    foreach ( $NAMEQ, $NAMER, $PER_IDY, $ALNLEN, $LEN1, $LEN2, $MISMATCH, $GAP, $S1, $E1, $S2, $E2, $EVALUE, $BITSCORE, $COVQ, $COVR, $LENQ, $LENR ) {
        s/\s+//g;
    }
    return ( $NAMEQ, $NAMER, $PER_IDY, $ALNLEN, $LEN1, $LEN2, $MISMATCH, $GAP, $S1, $E1, $S2, $E2, $EVALUE, $BITSCORE, $COVQ, $COVR, $LENQ, $LENR );
}

sub revcomp {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq;
}

sub esc_chars {
    my $temp_string = shift;

    # will change, for example, a!!a to a\!\!a
    $temp_string =~ s/([^\w\d\.\_\-])/\\$1/g;
    return $temp_string;
}

sub special_char_to_underscore {
    my $temp_string = shift;

    # will change, for example, a!!a to a\!\!a
    $temp_string =~ s/([^\w\d\.\_\-])/_/g;
    return $temp_string;
}
