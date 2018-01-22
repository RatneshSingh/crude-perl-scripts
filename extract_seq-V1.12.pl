#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;
use List::Util qw/ max min /;
#####################################################################################################
# This script is made to extract sequences from sequence file provided by the user, based on the    #
# cordinates given in the blast file (tabular output). user can provide number of nuceotide to add  #
# on both sides This option is optional.This script will read sequence name and cordinates from     #
# blast file or user provided table and extract sequences from sequence file provided.              #
# Alternatively names and coordinates can be assigned manually                                      #                       #
#                                                                                                   #
# Author : Ratnesh Singh                                                                            #
# version 1.13                                                                                      #
# contact for bugs: ratnesh@hawaii.edu                                                              #
#####################################################################################################

getopts('a:b:c:d:e:f:gh:i:jk:l:m:n:o:pq:rs:t:u');
## avialable characters:v w x y z
our($opt_a,$opt_b,$opt_c,$opt_d,$opt_e,$opt_f,$opt_g,$opt_h,$opt_i,$opt_j,$opt_k,$opt_l,$opt_m,$opt_n,$opt_o,$opt_p,$opt_q,$opt_r,$opt_s,$opt_t,$opt_u);
our (%seq,@coords,%chunks);

# like the shell getopt, "d:" means d takes an argument
print "-sequence file: $opt_s\n"                      if defined $opt_s;
print "-blast file/sequence: $opt_b\n"                if defined $opt_b;
print "-Table file: $opt_l\n"                         if defined $opt_l;
print "-add to head: $opt_h\n"                        if defined $opt_h;
print "-add to tail: $opt_t\n"                        if defined $opt_t;
print "-print output as: $opt_o\n"                    if defined $opt_o;
print "-input mode: $opt_m\n"                         if defined $opt_m;
print "-using names till first space while searching" if defined $opt_g;
print "-saving extracted sequence in append mode\n"   if defined $opt_u;
#print "-reverse: $opt_r\n" if defined $opt_r;
print "Unprocessed by Getopt::Std:\n" if $ARGV[0];
foreach (@ARGV) {
    print "-- $_\n";
}

my $help = "\n\nThis script will read sequence name and cordinates from
blast file or user provided table and extract sequences from sequence file provided.
Alternatively names and coordinates can be assigned manually.

usage:\n use following options
-m mode of input(auto or manual) [Defaulst is auto]
-s sequence file

-b Blast out put file in table format OR
-l table of SeqName\tStart\tEnd\tStrand seperated by Tab/spaces
-p extract from query instead of default subject.
-h (optional) number of extra nt to add at head [Default 0]
-t (optional) number of extra nt to add at tail [Default 0]
-o output file [Default is: output_seq_extract]
-g use names till first space while searching. [use full length names]
-d use this delimiter to split sequence name and use column[-c] to find seq.
-c column to match if using delimiter to split the sequence names.
\n\n\n usage for manual input mode:\n\n
-s sequence file name containing all sequences
-n sequence/contig name to extract from.
-i start site to cut from
-e end site to cut till
-k multiple coordinates to extract sequence of CDS eg:\"seqname1:start-end:strand,seqname2:start1-end1:strand,seqname3:start2-end2:strand\"...
-a Strand of sequence. minus strand will be reverse complemented.
-h (optional) number of extra nt to add at head [Default 0]
-t (optional) number of extra nt to add at tail [Default 0]
-o output file [Default is: output_seq_extract]
-u append extracted sequence in output file [overwrite]
-q query name for manual input[Default is: manual_query ]
-f seperater to use for attaching query name to coordinates for header[\"_\"].
-j do not add coordinates to header name. [add_coordinates]
-r ask if user have any more coordinates to extract.
    Saves time for loading fasta sequences for more than one extractions.
\n";

die "\nThere is no sequence file specified with -s \n $help" if !defined($opt_s);

# Check if blast or table file is provided with 'Auto' option.

if ($opt_b && $opt_l) {
  print "\nblast file (-b) and table(-l) cannot be used together\nSwitching to blast file only\n";
  undef $opt_l;
}



if ( !defined $opt_m ) {
    if (!$opt_b && !$opt_l){die "\nThere is no blast file specified with -b or table with -l \n $help" ;}
    print "\n option for -m flag is either missing or improper (not 'auto' or 'manual'). Program assuming default 'auto' option for file input\n";
    $opt_m = 'auto';
}

if(defined $opt_k){
    print "\n\nUsing multiple coordinates for extraction.\n";
    $opt_m='manual';
    die "-i and -e options are not valid with -k option." if ($opt_i || $opt_e);
    $opt_k=~s/\"//g;
    @coords=split /,/,$opt_k;
    #die "\nuneven number of coordinates. Make sure to provide start and end for each exon.\n" if ((scalar@coords)%2>0);

}



elsif ( lc($opt_m) eq 'auto' ) { die "\nThere is no blast file specified with -b \n $help" if (!$opt_b && !$opt_l); }
elsif ( lc($opt_m) eq 'manual' ) { $opt_m = 'manual'; }
else                             { print "Having problems with -m option. Check the options\n"; }

$opt_h = $opt_h ? $opt_h : 0;
$opt_t = $opt_t ? $opt_t : 0;                       # if !defined ($opt_t);
$opt_o = $opt_o ? $opt_o : 'output_seq_extract';    # if !defined ($opt_o);
$opt_f = defined $opt_f ? $opt_f : "_";
#$opt_m = 'auto' if !defined ($opt_m); # moved to above section
#$opt_q = $opt_q ? $opt_q : 'manual_query';          # if !defined ($opt_q);
$opt_c = $opt_c ? $opt_c : 1;

print "saving output in:$opt_o\n";
if   ( $opt_u ) { open OUT, ">>$opt_o"  or die "cannot open Output file\n"; }
else                      { open OUT, ">$opt_o" or die "cannot open Output file\n"; }

my$ref_fhash=fasta_to_hash($opt_s);
#%seq=%$ref_fhash;
print "\n Read ".scalar(keys %seq). " sequences from $opt_s\n";
###################################################################
#parse information from blast file and extract seq start and end  #
###################################################################
if   ( lc($opt_m) eq 'auto' ) { goto "AUTO"; }
else                          { goto "MANUAL"; }



AUTO: {

    my$bfile=$opt_b if $opt_b;
      $bfile=$opt_l if $opt_l;
    open BLAST, "$bfile" or die "cannot read blast file \n" ;

    my$count=0;
    while (<BLAST>) {
        my $line = $_;
        $count++;
        #	print "before line : $line\n";
        if ( $line =~ /^\s+$/ ) { print "\nBLAST line:$count is empty. Skipping\n"; next; }
        if ( $.==1 && $line =~ /query/i or /match/i or /score/ or /gap/ or /mismatch/ ) { print "\nBLAST line:$count is a header. Skipping\n"; next; }
        
        my($query,$subject,$sstart,$subend,$qstart,$qend,$strand_s);
        if ($opt_b) {
            #	print "line: $line\n";
            my @line_info = split( /\t/, $line );
    

            #	print "query:$query\n";
            $subject = $line_info[1];
            if ( defined $opt_g ) { my @names = split( /\s/, $subject ); $subject = $names[0];}

            #	print"subject:$subject\n";
            $sstart = $line_info[8];

            #	print "Extracting --> $subject: sstart: $sstart\t";
            $subend = $line_info[9];

            ### collect info from query too for $opt_p
            $query = $line_info[0];
            $qstart = $line_info[6] ;
            $qend = $line_info[7] ; 
        }elsif($opt_l){
            #	print "line: $line\n";
            my @line_info = split( /\s+/, $line ) if $opt_l;
        

            #	print "query:$query\n";
            $subject = $line_info[0] if $opt_l;;
            if ( defined $opt_g ) { my @names = split( /\s/, $subject ); $subject = $names[0];}

            #	print"subject:$subject\n";
            $sstart = $line_info[1] if $opt_l;

            #	print "Extracting --> $subject: sstart: $sstart\t";
            $subend = $line_info[2] if $opt_l;

            $strand_s = $line_info[3] if ($line_info[3] && $opt_l);
            ### collect info from query too for $opt_p
            my $query = $subject;
            my $qstart = $sstart;
            my $qend = $subend ; 
        }else{
            die "\n Please provide either a blast file (-b)  or a table (-l) to get coordinate info for sequence extraction.\n"
            
        }
        
                      
    
    
    

        my $name_header=$query;
        my $extract_header=$subject;
        my $start_f=$sstart;
        my $end_f=$subend;
        my $strand_f=$strand_s;
        
        if ($opt_b && $opt_p) {
                $name_header=$subject;
                $extract_header=$query;
                $start_f=$qstart;
                $end_f=$qend;
                $strand_f=$strand_s;
        }
        
        
        
        
        my ( $new_header, $new_sequence, $returnedlength ) = get_region( $name_header, $extract_header, $start_f, $end_f, $strand_f );
        if ( $new_sequence ) {
            print OUT">$new_header\n$new_sequence\n";
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

    if(!$opt_k){
        my ($sstart,$subend,$seq_name,$strand);

        if(!$opt_n){
        print "type the name of contig to look for:";
            $seq_name = <STDIN>;
            chomp($seq_name);
        }else{$seq_name=$opt_n}
        #### get sstart from commandline or ask user.
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
            $subend = <STDIN>;
            
            if ( $subend <= 0 ) { print "End site cannot be less than 1\n\n"; exit; }

        }else{$subend=$opt_e}
    
        
        if(!$opt_a){
            $strand="plus";
         }else{$strand=$opt_a}

        $strand="plus" if $strand eq '+';
        $strand="minus" if $strand eq '-';
        my $header = $opt_q?$opt_q:$seq_name;
        chomp($header,$seq_name,$sstart,$subend,$strand);
        ####
        my ( $new_header, $new_sequence ) = get_region( $header, $seq_name, $sstart, $subend, $strand );
            
            
        if ( $new_sequence ) {
            print OUT">$new_header\n$new_sequence\n";
            print "Extracted --> $seq_name: sstart: $sstart\t end: $subend \t Expected length :". ($subend - $sstart + 1) ."\t Extracted length:".length($new_sequence)." \n";
        }
        else {

            print "\nError in extracting.\n";
        }

        ###
        if($opt_r){
            print "type q  to quit or m to extract more:";
            my $option = <STDIN>;
            chomp($option);
            if   ( lc($option) eq 'm' ) { $opt_n=$opt_i=$opt_e=$opt_a=undef; goto "MANUAL"; }
        }
        else                        { last; }
    }
    elsif($opt_k){
        my $query = $opt_q if $opt_q;
        ####
        foreach my$coord(@coords){
            (my$query_i,my$start_s,my$end_s,my$strand_s)=get_coords($coord);
            $query=$opt_q?$opt_q:$query_i;
            my ( $new_header, $new_sequence) = get_region( $query, $query_i, $start_s, $end_s, $strand_s );
            print OUT">$new_header\n$new_sequence\n" if length($new_sequence)>0;
        }
    }


}

close(OUT);
close(BLAST);
exit;

################################################################################
# subroutine for extraction of sequences
################################################################################
sub fasta_to_hash{
    my$fasta_file=shift;
    
    print "reading sequence file:  $fasta_file ....\n";
    open FASTA, "$fasta_file" or die "cannot open Sequence file:$fasta_file\n";

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
    return(\%seq);
    
}



sub extract_seq {
    my ( $query1, $subject1, $start_s1, $end_s1, $len1, $opt1_h, $opt1_t, $strand ) = @_;

    my $new_sequence1 = substr( $seq{$subject1}{'sequence'}, $start_s1, $len1 ) if defined $seq{$subject1};

    if ( lc$strand eq 'minus' ) {
        $new_sequence1 = reverse_seq($new_sequence1);
        print "seq reversed\n\n";
    }

    #	my $length_seq = length($seq{$subject1}) if defined $seq{$subject1};
    my $length_seq = $seq{$subject1}{'len'} if defined $seq{$subject1};
    my $returnedlength = length($new_sequence1);
    print "$subject1 has : $length_seq nt \n" if defined $length_seq && !defined $opt_k;
    my$new_header1;
    if (!defined $opt_j){ my$head=$query1 eq $subject1?$query1:$subject1.$opt_f.$query1;$new_header1 = $head  . '-' . 'start-' . $start_s1 . '_' . 'end-' . $end_s1;}else{$new_header1 = $query1}

    return ( $new_header1, $new_sequence1, $returnedlength ) if defined $new_sequence1;
    print "Cannot find $subject1\n" if !defined $seq{$subject1}{'sequence'};
}

sub get_coords{
    my$info=shift;
    my($q,$se,$s,$e,$str)=undef;
    if ($info =~ m/\:/) {($q,$se,$str)=split /:/,$info}
    ($s,$e)=split /-/,$se;
    if (!$str) {
        if ($opt_a) {$str=$opt_a;}
        elsif($s < $e){$str='plus'}
        elsif($s > $e){$str='minus'}
    }
    
    chomp($q,$s,$e,$str);
    return($q,$s,$e,$str);
}


sub reverse_seq{
    my$rseq=shift;
    $rseq=reverse $rseq;
    $rseq =~ tr/atgcATGC/tacgTACG/;
    return(uc $rseq);
       
}

sub get_region{
    my($hsname, $ssname, $sstart, $send, $strand)=@_;


    if (!$strand) {
        if ($opt_a){$strand=$opt_a}         
        elsif ($send < $sstart){$strand='minus' }
        else{$strand='plus'}
    }
    
        $strand="plus" if $strand eq '+';
        $strand="minus" if $strand eq '-';
    
    
    $sstart= min($sstart,$send);$sstart=$sstart  - $opt_h;
    $send  = max($sstart,$send);$send=$send  + $opt_t;
    
    
    
    
    
    if (!$seq{$ssname}{'sequence'}) {
        print "Cannot find $ssname\n";
        return(undef,undef,undef);
    }
    else{
    
        if ($send > $seq{$ssname}{'len'}){print "\nEnd coord ($send)is larger than sequence length ($seq{$ssname}{'len'}). Changing it to the seq length.\n";$send=$seq{$ssname}{'len'}}
        if ($sstart < 1){print "\nStart coord ($sstart) is less than 1. Changing it to 1\n";$sstart=1 }
        my$len1=$send-$sstart+1;
   
        my $new_sequence1 = substr( $seq{$ssname}{'sequence'}, $sstart-1, $len1 ) if defined $seq{$ssname};

        if ( lc$strand eq 'minus' ) {$new_sequence1 = reverse_seq($new_sequence1);}

        #	my $length_seq = length($seq{$subject1}) if defined $seq{$subject1};
        my $length_seq = $seq{$ssname}{'len'} if defined $seq{$ssname};
        print "$ssname has : $length_seq nt \n" if defined $length_seq && !defined $opt_k;
        my$new_header1;
        if (!defined $opt_j){ my$head=$hsname eq $ssname?$hsname.$opt_f:$hsname.$opt_f.$ssname."-";$new_header1 = $head .'start-' . $sstart . '_' . 'end-' . $send;}else{$new_header1 = $hsname}
        print "Extracted --> $ssname: sstart: $sstart\t end: $send \t strand: $strand\tExpected length :". ($send - $sstart + 1) ."\t Extracted length:".length($new_sequence1)." \n*****\n";
        return ( $new_header1, $new_sequence1);
    }
}
