#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;

#####################################################################################################
# This script is made to extract sequences from sequence file provided by the user, based on the    #
# cordinates given in the blast file (tabular output). user can provide number of nuceotide to add  #
# on both sides This option is optional                                                             #
#                                                                                                   #
# Author : Ratnesh Singh                                                                            #
# version 1.10                                                                                       #
# contact for bugs: ratnesh@hawaii.edu                                                              #
#####################################################################################################

getopt('abcdefghijkmnoqrst');
## avialable characters:k l p u v w x y z
our ( $opt_s, $opt_b, $opt_h, $opt_t, $opt_o, $opt_m, $opt_q, $opt_g, $opt_d, $opt_c, $opt_e, $opt_i,$opt_n,$opt_r,$opt_a,$opt_f,$opt_j,$opt_k );
our (%seq,@coords);

# like the shell getopt, "d:" means d takes an argument
print "-sequence file: $opt_s\n"                      if defined $opt_s;
print "-blast file/sequence: $opt_b\n"                if defined $opt_b;
print "-add to head: $opt_h\n"                        if defined $opt_h;
print "-add to tail: $opt_t\n"                        if defined $opt_t;
print "-print output as: $opt_o\n"                    if defined $opt_o;
print "-input mode: $opt_m\n"                         if defined $opt_m;
print "-using names till first space while searching" if defined $opt_g;

#print "-reverse: $opt_r\n" if defined $opt_r;
print "Unprocessed by Getopt::Std:\n" if $ARGV[0];
foreach (@ARGV) {
    print "-- $_\n";
}

my $help = "\n\nThis script will read sequence name and cordinates from
blast file and extract sequences from sequence file provided

usage:\n use following options\n -m mode of input(auto or manual) [Defaulst is auto]
-s sequence file \n -b Blast out put file in table format
-h (optional) number of extra nt to add at head [Default 0]
-t (optional) number of extra nt to add at tail [Default 0]
-o output file [Default is: output_seq_extract]
-g use names till first space while searching. [use full length names while]
-d use this delimiter to split sequence name and use column[-c] to find seq.
-c column to match if using delimiter to split the sequence names.
\n\n\n usage for manual input mode:\n\n-s sequence file name containing all sequences
-n sequence/contig name to extract from.
-i start site to cut from
-e end site to cut till
-k multiple coordinates to extract sequence of CDS eg:\"start-end,start1-end1,start2-end2\"...
-a Strand of sequence. minus strand will be reverse complemented.
-h (optional) number of extra nt to add at head [Default 0]
-t (optional) number of extra nt to add at tail [Default 0]
-o output file [Default is: output_seq_extract]
-q query name for manual input[Default is: manual_query ]
-f seperater to use for attaching query name to coordinates for header[\"_\"].
-j do not add coordinates to header name. [0]
-r ask if user have any more coordinates to extract.
    Saves time for loading fasta sequences for more than one extractions.
\n";

die "\nThere is no sequence file specified with -s \n $help" if !defined($opt_s);

# Check if blast file is provided with 'Auto' option.
if ( !defined $opt_m ) {
    die "\nThere is no blast file specified with -b \n $help" if !defined($opt_b);
    print "\n option for -m flag is either missing or improper (not 'auto' or 'manual'). Program assuming default 'auto' option for file input\n";
    $opt_m = 'auto';
}

if(defined $opt_k){
    print "\n\nUsing multiple coordinates for extraction.\n";
    $opt_m='manual';
    die "-i and -e options are not valid with -k option." if ($opt_i || $opt_e);
    $opt_k=~s/\"//g;
    @coords=split /,|\s+|-/,$opt_k;
    die "\nuneven number of coordinates. Make sure to provide start and end for each exon.\n" if ((scalar@coords)%2>0);

}



elsif ( lc($opt_m) eq 'auto' ) { die "\nThere is no blast file specified with -b \n $help" if !defined($opt_b); }
elsif ( lc($opt_m) eq 'manual' ) { $opt_m = 'manual'; }
else                             { print "Have problems with -m option. Check the options\n"; }

$opt_h = $opt_h ? $opt_h : 0;
$opt_t = $opt_t ? $opt_t : 0;                       # if !defined ($opt_t);
$opt_o = $opt_o ? $opt_o : 'output_seq_extract';    # if !defined ($opt_o);

#$opt_m = 'auto' if !defined ($opt_m); # moved to above section
$opt_q = $opt_q ? $opt_q : 'manual_query';          # if !defined ($opt_q);
$opt_c = $opt_c ? $opt_c : 1;

print "saving output in:$opt_o\n";
if   ( $opt_m eq 'auto' ) { open OUT, ">$opt_o"  or die "cannot open Output file\n"; }
else                      { open OUT, ">>$opt_o" or die "cannot open Output file\n"; }

print "reading sequence file....Plz wait\n";
open FASTA, "$opt_s" or die "cannot open Sequence file:$opt_s\n";

$/ = "\n>";
while (<FASTA>) {
    chomp;

    my ( $header, @sequence ) = split( /\n/, $_ );
    $header =~ s/>//;

    if ( defined $opt_g ) { my @names = split( /\s/,       $header ); $header = $names[0]; }
    if ( defined $opt_d ) { my @names = split( /\Q$opt_d/, $header ); $header = $names[ $opt_c - 1 ]; }
    $header=~s/^\s+//g;
    $header=~s/\s+$//g;

    my $sequence = join( "", @sequence );
    $sequence =~ s/\s+//g;
    $sequence =~ s/\n+//g;

    $seq{$header}{'sequence'} = $sequence;
    $seq{$header}{'len'}      = length($sequence);

}
close(FASTA);
print ".............Done\n";
$/ = "\n";

###################################################################
#parse information from blast file and extract seq start and end  #
###################################################################
if   ( lc($opt_m) eq 'auto' ) { goto "AUTO"; }
else                          { goto "MANUAL"; }



AUTO: {

    open BLAST, "$opt_b" or die "cannot read blast file \n";
    my$count=0;
    while (<BLAST>) {
        my $line = $_;
        $count++;
        #	print "before line : $line\n";
        if ( $line =~ /^\s+$/ ) { print "\nBLAST line:$count is empty. Skipping\n"; next; }
        if ( $.==1 && $line =~ /query/i or /match/i or /score/ or /gap/ or /mismatch/ ) { print "\nBLAST line:$count is a header. Skipping\n"; next; }

        #	print "line: $line\n";
        my @line_info = split( /\t/, $line );
        my $query = $line_info[0];

        #	print "query:$query\n";
        my $subject = $line_info[1];
        if ( defined $opt_g ) { my @names = split( /\s/, $line_info[1] ); $line_info[1] = $names[0]; }

        #	print"subject:$subject\n";
        my $sstart = $line_info[8];

        #	print "Extracting --> $subject: sstart: $sstart\t";
        my $subend = $line_info[9];

        #	print " subend:$subend\n";

        #adjustment for increase in length at both ends

        my ( $start_s, $end_s, $strand );

        if   ( $sstart > $subend ) { $end_s   = $sstart + $opt_h - 1; $start_s = $subend - $opt_t - 1; $strand = 'minus'; }
        else                       { $start_s = $sstart - $opt_h - 1; $end_s   = $subend + $opt_t - 1; $strand = 'plus' }

        if ( $start_s < 0 ) { $start_s = 0; }

        # added to avoid negative values of start_s.

        my $len = $end_s - $start_s + 1;

        #print "Extracting --> $subject: sstart: $start_s\t end: $end_s \t length : $len \n";

        #my ($new_header,$new_sequence)=extract_seq($query,$subject,$start_s,$end_s,$len,$opt_h,$opt_t,$strand);
        ##print OUT">$new_header.$line_info[0].$line_info[1].$line_info[2].$line_info[3].$line_info[4].$line_info[5].$line_info[6].$line_info[7].$line_info[8].$line_info[9] \n$new_sequence\n" if defined $new_sequence;
        #print OUT">$new_header\n$new_sequence\n" if defined $new_sequence;

        if ( my ( $new_header, $new_sequence, $returnedlength ) = extract_seq( $query, $subject, $start_s, $end_s, $len, $opt_h, $opt_t, $strand ) ) {

#print OUT">$new_header.$line_info[0].$line_info[1].$line_info[2].$line_info[3].$line_info[4].$line_info[5].$line_info[6].$line_info[7].$line_info[8].$line_info[9] \n$new_sequence\n" if defined $new_sequence;
            print OUT">$new_header\n$new_sequence\n" if defined $new_sequence;

            #print "Header:$new_header\nSequence:$new_sequence\n" if defined $new_sequence;
            print "Extracted --> $subject: sstart: $start_s\t end: $end_s \t Expected length : $len\t Extracted length:$returnedlength \n";
        }
        else {

            print "\nError in extracting.\n";
        }

    }
    exit;
}

##################################################################################################
# for manual input method.
MANUAL: {
    #	my $subject=$opt_b if defined ($opt_b);
    #	if(!defined($opt_b)){print "type the name of contig to look for:";$subject=<STDIN>;	chomp($subject);};
    my ($sstart,$subend,$subject,$strand);

        if(!$opt_n){
        print "type the name of contig to look for:";
            $subject = <STDIN>;
            chomp($subject);
        }else{$subject=$opt_n}
        #### get sstart from commandline or ask user.
    if(!$opt_k){
        if(!$opt_i){
            print "print start site to cut:" ;
            $sstart = <STDIN>;
            if ( $sstart <= 0 ) { print "Start site cannot be less than 1\n\n"; exit; }
            chomp($sstart);
        }
        else{$sstart=$opt_i}

        #### get subend from commandline or ask user.
        if(!$opt_e){
            print "\nPrint subject end to cut till:";
            my $subend = <STDIN>;
            if ( $subend <= 0 ) { print "End site cannot be less than 1\n\n"; exit; }

        }else{$subend=$opt_e}

        if(!$opt_a){
            $strand="plus";
         }else{$strand=$opt_a}

        $strand="plus" if $strand eq '+';
        $strand="minus" if $strand eq '-';


        my ( $start_s, $end_s );
        if   ( $sstart > $subend ) { $end_s   = $sstart + $opt_h - 1; $start_s = $subend - $opt_t - 1; }
        else                       { $start_s = $sstart - $opt_h - 1; $end_s   = $subend + $opt_t - 1; }
        my $len = $end_s - $start_s + 1;
        print "Extracting --> $subject: sstart: $start_s\t end: $end_s \t length : $len \n";
        my $query = $opt_q;

        ####
        if ( my ( $new_header, $new_sequence, $returnedlength ) = extract_seq( $query, $subject, $start_s, $end_s, $len, $opt_h, $opt_t,$strand ) ) {

    #print OUT">$new_header.$line_info[0].$line_info[1].$line_info[2].$line_info[3].$line_info[4].$line_info[5].$line_info[6].$line_info[7].$line_info[8].$line_info[9] \n$new_sequence\n" if defined $new_sequence;
            print OUT">$new_header\n$new_sequence\n" if defined $new_sequence;

            #print "Header:$new_header\nSequence:$new_sequence\n" if defined $new_sequence;
            print "Extracted --> $subject: sstart: $start_s\t end: $end_s \t Expected length : $len\t Extracted length:$returnedlength \n";
        }
        else {

            print "\nError in extracting.\n";
        }

        ###
        if($opt_r){
        print "type q  to quit or m to extract more:";
        my $option = <STDIN>;
        chomp($option);
        if   ( lc($option) eq 'm' ) { goto "MANUAL"; }
        }
        else                        { last; }
    }
    elsif(defined$opt_k){
        my $query = $opt_q;
        $opt_j=1;
        my $part_seq;

        if(!$opt_a){$strand="plus";}else{$strand=$opt_a}

        $strand="plus" if $strand eq '+';
        $strand="minus" if $strand eq '-';
        ####

        @coords=sort{$a<=>$b}@coords;


        for(my$i=0;$i<scalar@coords;$i=$i+2){

            my ( $start_s, $end_s );
            if   ( $coords[$i] > $coords[$i+1] ) { $end_s   = $coords[$i] - 1; $start_s = $coords[$i+1] - 1; }
            else                                 { $start_s = $coords[$i] - 1; $end_s   = $coords[$i+1] - 1; }
            my $len = $end_s - $start_s + 1;
            my ( $new_header, $new_sequence, $returnedlength ) = extract_seq( $query, $subject, $start_s, $end_s, $len, 0, 0,'plus' ) ;

            $part_seq.=$new_sequence;
        }
            if ( $strand eq 'minus' ) {
            my $new_sequence2 = reverse $part_seq;
            $new_sequence2 =~ tr/atgcATGC/tacgTACG/;
            $part_seq = uc $new_sequence2;}
        print OUT">$query coord:$opt_k\n$part_seq\n" if length($part_seq)>0;

    }
}



close(OUT);
close(BLAST);
exit;

################################################################################
# subroutine for extraction of sequences
################################################################################

sub extract_seq {
    my ( $query1, $subject1, $start_s1, $end_s1, $len1, $opt1_h, $opt1_t, $strand ) = @_;

    my $new_sequence1 = substr( $seq{$subject1}{'sequence'}, $start_s1, $len1 ) if defined $seq{$subject1};

    if ( $strand eq 'minus' ) {
        my $new_sequence2 = reverse $new_sequence1;
        $new_sequence2 =~ tr/atgcATGC/tacgTACG/;
        $new_sequence1 = uc $new_sequence2;
        print "seq reversed\n\n";
    }

    #	my $length_seq = length($seq{$subject1}) if defined $seq{$subject1};
    my $length_seq = $seq{$subject1}{'len'} if defined $seq{$subject1};
    my $returnedlength = length($new_sequence1);
    print "$subject1 has : $length_seq nt \n" if defined $length_seq && !defined $opt_k;
    my$new_header1;
    if (!defined $opt_j){ $new_header1 = $query1 . $opt_f . $subject1 . '-' . 'start-' . $start_s1 . '_' . 'end-' . $end_s1;}else{$new_header1 = $query1}

    return ( $new_header1, $new_sequence1, $returnedlength ) if defined $new_sequence1;
    print "Cannot find $subject1\n" if !defined $seq{$subject1}{'sequence'};
}
