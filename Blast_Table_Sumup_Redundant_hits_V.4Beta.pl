#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;

#####################################################################################################
# This script is made to remove overlapping regions in hits and find unique coordinates to          #
# identify hit region. Also pulls out the hit region or whole sequence fro subject to align later   #
#
# Author : Ratnesh Singh                                                                            #
# version 1.2                                                                                       #
# contact for bugs: ratnesh@hawaii.edu                                                              #
# last updated: 01/03/2018
#####################################################################################################

our ( $opt_b, $opt_o, $opt_s, $help, $opt_pd, $opt_pc, $opt_v, $opt_and, $opt_or, $opt_e, $opt_l, $opt_p, $opt_qm, $opt_qx, $opt_sm, $opt_sx, $opt_final_qcoverage_max,
    $opt_final_qcoverage_min, $opt_final_length, $vectrim,$Query_New,$opt_nofilter,$manual_col,$opt_col_qc,$opt_col_sc,$opt_col_ql,$opt_col_sl,$opt_fs);

#$opt_b="ALL_MADS-Box.SRF-TF.50OrMore_tblastn_MADSmaskedLotusGenome.tblastn.1e-5.table";
#my $opt_v                   = "No";
my $opt_a    = 0;
my $opt_d    = " ";
my $opt_c    = 1;
my $opt_h    = 0;
my $opt_t    = 0;
my $from_end = 100;

#my $opt_e                   = 10;
#my $opt_l                   = 1;
#my $opt_p                   = 1;
my $opt_k = 'subject';

#my $opt_qm                  = 0;
#my $opt_qx                  = 100;
#my $opt_sm                  = 0;
#my $opt_sx                  = 100;
#my $opt_final_qcoverage_min = 0;
#my $opt_final_qcoverage_max = 101;
#my $opt_final_length = 1;
$opt_and = 'true';

#my $opt_filter              = 'all'
my ( %blast, %uniq_hit, %sequence );

#getopt('bodcsvaehtlpkqj');

my $result = GetOptions(
    "blastfile|b=s"         => \$opt_b,
    "seq|s=s"               => \$opt_s,
    "fullseq|fs"            => \$opt_fs,
    "verbose|v"             => \$opt_v,
    "allowed|a=s"           => \$opt_a,
    "delimiter|d=s"         => \$opt_d,
    "column|c=i"            => \$opt_c,
    "pat_delim|pd=s"        => \$opt_pd,
    "pat_col|pc=i"          => \$opt_pc,
    "head|hd=i"             => \$opt_h,
    "tail|t=i"              => \$opt_t,
    "evalue|e=f"            => \$opt_e,
    "aln_len|l=i"           => \$opt_l,
    "perid|p=f"             => \$opt_p,
    "target|k=s"            => \$opt_k,
    "qcov_max|qcx=f"        => \$opt_qx,
    "qcov_min|qcm=f"        => \$opt_qm,
    "scov_min|scm=f"        => \$opt_sm,
    "scov_max|scx=f"        => \$opt_sx,
    "final_qcov_min|fqcm=f" => \$opt_final_qcoverage_min,
    "final_qcov_max|fqcx=f" => \$opt_final_qcoverage_max,
    "final_length|fl=i"     => \$opt_final_length,
    "vectrim|vt"            => \$vectrim,
    "from_end|fe=i"         => \$from_end,
    "o"                     => \$opt_o,
    "and"                   => \$opt_and,
    "or"                    => \$opt_or,
    "no_filter|nf"          =>\$opt_nofilter,
    "col_qc|cqc=i"          =>\$opt_col_qc,
    "col_sc|csc=i"          =>\$opt_col_sc,
    "col_ql|cql=i"          =>\$opt_col_ql,
    "col_sl|csl=i"          =>\$opt_col_sl,
    "help"                  => \$help
);

$manual_col='TRUE' if ($opt_col_qc||$opt_col_sc||$opt_col_ql||$opt_col_sl);


##################################################
# set filtering option scope (all or any)
if ($opt_or) {
    undef $opt_and;
    if ( !defined( $opt_e || $opt_l || $opt_p || $opt_qm || $opt_qx || $opt_sm || $opt_sx ) ) {
        $opt_p = 1;
        print "\nSetting percent identity to 1 as none of the options in e,l,p,qm,qx,sm and sx is provided\n" if $opt_v;
    }
    if ( !defined( $opt_final_qcoverage_max || $opt_final_qcoverage_min || $opt_final_length ) ) {

        $opt_final_length = 1;
    }
    $opt_e                   = -1          if !defined $opt_e;                      #&& $opt_or;
    $opt_l                   = 10000000000 if !defined $opt_l;                      #&& $opt_or;
    $opt_p                   = 101         if !defined $opt_p;                      #&& $opt_or;
    $opt_qm                  = 101         if !defined $opt_qm;                     #&& $opt_or;
    $opt_qx                  = -1          if !defined $opt_qx;                     #&& $opt_or;
    $opt_sm                  = 101         if !defined $opt_sm;                     #&& $opt_or;
    $opt_sx                  = -1          if !defined $opt_sx;                     #&& $opt_or;
    $opt_final_qcoverage_max = -1          if !defined $opt_final_qcoverage_max;    #&& $opt_or;
    $opt_final_qcoverage_min = 101         if !defined $opt_final_qcoverage_min;    #&& $opt_or;
    $opt_final_length        = 10000000000 if !defined $opt_final_length;           #&& $opt_or
}

# set some values for filtering due to different expectations for and and or situtaions
$opt_e                   = 10  if !defined $opt_e                   && $opt_and;
$opt_l                   = 1   if !defined $opt_l                   && $opt_and;
$opt_p                   = 1   if !defined $opt_p                   && $opt_and;
$opt_qm                  = 0   if !defined $opt_qm                  && $opt_and;
$opt_qx                  = 100 if !defined $opt_qx                  && $opt_and;
$opt_sm                  = 0   if !defined $opt_sm                  && $opt_and;
$opt_sx                  = 100 if !defined $opt_sx                  && $opt_and;
$opt_final_qcoverage_max = 100 if !defined $opt_final_qcoverage_max && $opt_and;
$opt_final_qcoverage_min = 0   if !defined $opt_final_qcoverage_min && $opt_and;
$opt_final_length        = 1   if !defined $opt_final_length        && $opt_and;

#$opt_e                   = -1          if !defined $opt_e                   && $opt_or;
#$opt_l                   = 10000000000 if !defined $opt_l                   && $opt_or;
#$opt_p                   = 101         if !defined $opt_p           		&& $opt_or;
#$opt_qm                  = 101         if !defined $opt_qm                  && $opt_or;
#$opt_qx                  = -1          if !defined $opt_qx                  && $opt_or;
#$opt_sm                  = 101         if !defined $opt_sm                  && $opt_or;
#$opt_sx                  = -1          if !defined $opt_sx                  && $opt_or;
#$opt_final_qcoverage_max = -1          if !defined $opt_final_qcoverage_max && $opt_or;
#$opt_final_qcoverage_min = 101         if !defined $opt_final_qcoverage_min && $opt_or;
##################################################

my $usage = "\n\nThis script will read the blast file in table format and
and combine multiple overlapping hits in to one hit and will extract the sequence from sequence file
if asked.

usage:perl script -options.
-b  blast file in table format.
-o  save the resulting hit table and sequences.
-a  allowed distance between two hit region on same sequence.
    hits occuring less that this distance apart will be considered
    overlapping and will be combined under one large hit [0].
-s  Sequence file to extract hit regions from.
-fs Pull the whole sequence out instead of extracting.
-col_qc|cqc  column number containing query coverage info[13]
-col_sc|csc  column number containing subject coverage info[14]
-col_ql|cql  column number containing query length info[15]
-col_sl|csl  column number containing subject length info[16]
Blast filtering options:
# These filters act on each hit during summarization.
-nf Do not use any filters.
-e  E value. Exclude results higher that the evalue [10]
-l  Aln_length. Exclude results with alignment length smaller than this[1]
-p  Percent_identity. Exclude results with per id lower than this value[0]
-qcov_max|qcx	Exclude hits with query coverage value greater than this value[100]
-qcov_min|qcm	Exclude hits with query coverage value lesser than this[0]
-scov_max|scx	Exclude hits with subject coverage value greater than this value[100]
-scov_min|scm	Exclude hits with subject coverage value lesser than this[0]

# these final filters act on summarized results
-final_qcov_min|fqcm	Remove hits with final qcoverage smaller than this[0]
-final_qcov_max|fqcx	Remove hits with final qcoverage larger than this[100]
-final_length|fl	Include only hits with final alignment length larger than this[1]

-and	include hit if all of the filtering conditions are true [true]
-or		include hit if any of the filtering condition is true

Other options:
-k  query|subject. Key to use for comparison [subject].
-d  delimiter to split the genomic header and use -c column as sequence names.
-c  Column number to be used as genomic sequence header after splitting by -d delimiter.
-v  Yes|No. Print detailed report [No]
-hd  Extra sequence to extract from upstream from hit region
-t  Extra sequence to extract from downstream from hit region
-vectrim	Trim vector ends based on blast results.
-from_end|fe	Maximum allowed distance of vector region from ends.
-help|h	Show usage help.
\n";

die "\n\n$usage\n\n" if $help;

if ( lc $opt_k ne 'subject' && lc $opt_k ne 'query' ) {
    die "*******Error: -k only accepts 'subject' or 'query' as options. \n";
}
die "\nThere is no blast file specified with -b \n $usage" if !defined($opt_b);
####################################################################
# set file handles for printing.
my $fh_tabout     = \*STDOUT;
my $fh_seqout     = \*STDOUT;
my $fh_vectrimout = \*STDOUT;
my $outfile;
if ($opt_o) {
    open TABOUT, ">$opt_b.$opt_k.summary.table";
    $fh_tabout = \*TABOUT;
}

if ( $opt_s && $opt_o ) {
    $outfile="$opt_b.extracted_seq.fasta";
    $outfile="$opt_b.full_seq.fasta" if $opt_fs;
    open SEQOUT, ">$outfile";
    $fh_seqout = \*SEQOUT;
}

if ($vectrim) {
    open VECTRIMOUT, ">$opt_b.vectrimmed.fasta";
    $fh_vectrimout = \*VECTRIMOUT;
}
####################################################################

###################################################################
# parse information from blast file                               #
###################################################################

open BLAST, "$opt_b" or die "cannot read blast file \n";

print STDOUT "Reading blast file....\n";
while (<BLAST>) {
    my $line = $_;
    chomp($line);
    if ( $line =~ /^\s*$/ ) { next }
    if ( $line =~ /query/i or /match/i or /score/ or /gap/ or /mismatch/ ) {
        next;
    }

    my ( $Query, $Subject_New, $Identity, $Aln_length, $Mismatch, $Gap, $Q_start, $Q_end, $New_S_start, $New_S_end, $E_value, $Bit_score,$qcov,$scov,$qlen,$slen,@ext) =
      split( /\s+/, $line );
    if ($manual_col) {
        my@extinfo=($qcov,$scov,$qlen,$slen,@ext);
        $qcov=$extinfo[$opt_col_qc - 13] if $opt_col_qc;
        $scov=$extinfo[$opt_col_sc - 13] if $opt_col_sc;
        $qlen=$extinfo[$opt_col_ql - 13] if $opt_col_ql;
        $slen=$extinfo[$opt_col_sl - 13] if $opt_col_sl;
    }
 #########################################################################################################
    # filter the blast results for query($opt_q) aln_length($opt_l),e-value($opt_e),percent_identity($opt_p)
    #### assign the max(start, end) to length values if table does not contain them
    $qlen||=max($Q_start, $Q_end);
    $slen||=max($New_S_start, $New_S_end);
    ### retain the longest value as query and subject length if not in the table.
    $sequence{$Query}{'len'}       = $qlen if (!$sequence{$Query}{'len'}||$sequence{$Query}{'len'} < $qlen);
    $sequence{$Subject_New}{'len'} = $slen if (!$sequence{$Subject_New}{'len'}||$sequence{$Subject_New}{'len'} < $slen);

    $qcov||=$Aln_length*100/$sequence{$Query}{'len'};
    $scov||=$Aln_length*100/$sequence{$Subject_New}{'len'};
    chomp( $Query, $Subject_New, $Identity, $Aln_length, $Mismatch, $Gap, $Q_start, $Q_end, $New_S_start, $New_S_end, $E_value, $Bit_score, $qcov, $scov, $qlen, $slen );

#if (!$opt_nofilter) {

    if (!$opt_nofilter &&
        $opt_or    # proceed if hit passed any of the filters
        && !( $Identity >= $opt_p || $Aln_length >= $opt_l || $E_value <= $opt_e || $qcov >= $opt_qm || $qcov <= $opt_qx || $scov <= $opt_sx || $scov >= $opt_sm )
      )
    {

        print "\nAnyTrue:Skipping this because none of these is true\n
    	  pid: $Identity >= $opt_p
		  alnlen: $Aln_length >= $opt_l
    	  eval:$E_value <= $opt_e
          qcovm:$qcov >= $opt_qm
		  qcovx: $qcov <= $opt_qx
          scovm: $scov >= $opt_sm
		  scovx: $scov <= $opt_sx
    	  " if $opt_v;

        next;
    }
    elsif (!$opt_nofilter &&
        $opt_and    # proceed only if hit passed all the filters
        && (   $Identity < $opt_p
            || $Aln_length < $opt_l
            || $E_value > $opt_e
            || $qcov < $opt_qm
            || $qcov > $opt_qx
            || $scov > $opt_sx
            || $scov < $opt_sm )
      )
    {

        print "\nAllTrue:Skipping this because at least one of these is  true
    	pid: $Identity < $opt_p
    	aln_len: $Aln_length < $opt_l
    	eval: $E_value > $opt_e
    	qcov_min: $qcov < $opt_qm
    	qcov_max: $qcov > $opt_qx
    	scov_min: $scov > $opt_sx
    	scov_max: $scov < $opt_sm " if $opt_v;

        next;
    }
    else {
        print "\nThis Hit passed all or atleast 1 of the filter.
		$Query\t$Subject_New
    	pid: $Identity >= $opt_p
    	aln_len: $Aln_length >= $opt_l
    	eval: $E_value <= $opt_e
    	qcov_min: $qcov >= $opt_qm
    	qcov_max: $qcov <= $opt_qx
    	scov_min: $scov <= $opt_sx
    	scov_max: $scov >= $opt_sm " if $opt_v;

    }

    #######################################################################################################
#}
    my $current_start = 0;
    my $current_end   = 0;
    my $Subject       = ();
    if ( lc $opt_k eq 'subject' ) {
        $current_start = $New_S_start;
        $current_end   = $New_S_end;
        $Subject       = $Subject_New;
		$Query_New=$Query;
    }
    elsif ( lc $opt_k eq 'query' ) {
        $current_start = $Q_start;
        $current_end   = $Q_end;
        $Subject       = $Query;
		$Query_New     = $Subject_New;
    }

    # process Sequence name if asked to
    if ( $opt_pd && $opt_pc ) {
        $Subject =~ s/^\s+//g;
        my @name = split( /$opt_pd/, $Subject );
        $Subject = $name[ $opt_pc - 1 ];
        $Subject =~ s/\s+//g;
    }

    #find the strand orientation of hits (plus or minus)
    my $strand;
    if   ( $current_start > $current_end ) { $strand = 'minus' }
    else                                   { $strand = 'plus' }

    #make the larger value end and lower value start.
    my $temp_send   = $current_end;
    my $temp_sstart = $current_start;

    $temp_send =~ s/\D//;
    $temp_sstart =~ s/\D//;
    $current_start = min( $temp_send, $temp_sstart );
    $current_end = max( $temp_send, $temp_sstart );
    my $hit_length = $current_end - $current_start;

    print "******\nComparing hit $current_start..$current_end region on $Subject\n"
      if $opt_v;

    #check if the hash element for the subject exists.
    if ( exists $blast{$Subject}{1}{'old_sstart'} ) {
        my $hit_found = "No";
        print "Other hits exists on $Subject.comparing sstart and send of new hit to hit from\n"
          if $opt_v;
        for ( my $hit = 1; $hit <= $uniq_hit{$Subject}{'hit'}; $hit++ ) {

            last if $hit_found eq "Yes";

            print "Region $hit:$blast{$Subject}{$hit}{'old_sstart'}:$blast{$Subject}{$hit}{'old_send'} with \n new hit:$current_start:$current_end\n"
              if $opt_v;

            #print "In side for loop current hit value:$hit\n";
            #print "In side for loop Max hit value:$uniq_hit{$Subject}{'hit'}\n";

            #if new hit is included in the existing hit-> Next;
            if ( ( ( $blast{$Subject}{$hit}{'old_sstart'} <= $current_start ) && ( $blast{$Subject}{$hit}{'old_send'} >= $current_end ) ) ) {
                print "New hit is included in the old hit\n" if $opt_v;
                $hit_found = "Yes";
                next;
            }

            #if new hit is partially/fully overlapping to existing hit; change the start and end values
            #1:check if start of new hit is in between the start and end of existing hit.
            #2: check if end of new hit is in between the start and end of existing hit.

            elsif (
                does_it_overlap( $blast{$Subject}{$hit}{'old_sstart'}, $blast{$Subject}{$hit}{'old_send'}, $current_start, $current_end ))
            {
                print "New hit overlaps with old hit. Reassigning start and End values\n"
                  if $opt_v;
                $blast{$Subject}{$hit}{'old_sstart'} = min( $blast{$Subject}{$hit}{'old_sstart'}, $current_start, $blast{$Subject}{$hit}{'old_send'}, $current_end );
                $blast{$Subject}{$hit}{'old_send'} = max( $blast{$Subject}{$hit}{'old_sstart'}, $current_start, $blast{$Subject}{$hit}{'old_send'}, $current_end );
                print "New values after reassignment:\nold:$blast{$Subject}{$hit}{'old_sstart'}:$blast{$Subject}{$hit}{'old_send'}\n "
                  if $opt_v;

                if ( $hit_length > $blast{$Subject}{$hit}{'largest_hit_length'} ) {
                    $blast{$Subject}{$hit}{'largest_hit_length'} = $hit_length;
                    $blast{$Subject}{$hit}{'strand'}             = $strand;
                }

                $hit_found = "Yes";
                next;
            }

            # check if new group includes existing group.
            elsif (( $blast{$Subject}{$hit}{'old_sstart'} >= $current_start )
                && ( $blast{$Subject}{$hit}{'old_send'} <= $current_end ) )
            {
                print "New hit includes old hit and is larger... Replacing old hit with new hit\n"
                  if $opt_v;
                $blast{$Subject}{$hit}{'old_sstart'} = min( $blast{$Subject}{$hit}{'old_sstart'}, $current_start, $blast{$Subject}{$hit}{'old_send'}, $current_end );
                $blast{$Subject}{$hit}{'old_send'} = max( $blast{$Subject}{$hit}{'old_sstart'}, $current_start, $blast{$Subject}{$hit}{'old_send'}, $current_end );

                if ( $hit_length > $blast{$Subject}{$hit}{'largest_hit_length'} ) {
                    $blast{$Subject}{$hit}{'largest_hit_length'} = $hit_length;
                    $blast{$Subject}{$hit}{'strand'}             = $strand;
                }

                $hit_found = "Yes";
                next;
            }
            else {
                next;
            }
        }

        # else create new category for new hit
        if ( $hit_found eq "No" ) {
            $uniq_hit{$Subject}{'hit'}++;
            print "No over lap found.....Creating new hit region:" . ( $uniq_hit{$Subject}{'hit'} ) . "  for $Subject\n"
              if $opt_v;
            $blast{$Subject}{ $uniq_hit{$Subject}{'hit'} }{'old_sstart'}         = $current_start;
            $blast{$Subject}{ $uniq_hit{$Subject}{'hit'} }{'old_send'}           = $current_end;
            $blast{$Subject}{ $uniq_hit{$Subject}{'hit'} }{'strand'}             = $strand;
            $blast{$Subject}{ $uniq_hit{$Subject}{'hit'} }{'largest_hit_length'} = $hit_length;
            next;
        }
    }

    # if $blast{$Subject} does not exists,create new hash element for the subject for the first time.
    else {
        print "Hit region does not exists. Creating hit region:1 on $Subject.\n"
          if $opt_v;
        $uniq_hit{$Subject}{'hit'} = 1;
        $blast{$Subject}{1}{'old_sstart'} = min( $current_end, $current_start );
        $blast{$Subject}{1}{'old_send'}           = max( $current_end, $current_start );
        $blast{$Subject}{1}{'strand'}             = $strand;
        $blast{$Subject}{1}{'largest_hit_length'} = $hit_length;
        next;
    }

}

print "...Done\n";

############################################################################################
# Reanalyze newly created regions to check if they overlap to other regions.
print "\n\n**************Reanalyzing regions for possible overlap as newly created regions might overlap after merging several small hits.*************\n\n"
  if $opt_v;
foreach my $subject ( keys %blast ) {

    # store the hit numbers in array @hits
    my @hits = keys %{ $blast{$subject} };

    for ( my $hit = 0; $hit < @hits; $hit++ ) {
        for ( my $nexthit = $hit + 1; $nexthit < @hits; $nexthit++ ) {

            if (
                does_it_overlap(
                    $blast{$subject}{ $hits[$hit] }{'old_sstart'},
                    $blast{$subject}{ $hits[$hit] }{'old_send'},
                    $blast{$subject}{ $hits[$nexthit] }{'old_sstart'},
                    $blast{$subject}{ $hits[$nexthit] }{'old_send'}
                )
              )
            {
                print
                  "Overlap found betwen the $subject hit:$hits[$hit] (start:$blast{$subject}{$hits[$hit]}{'old_sstart'}, End:$blast{$subject}{$hits[$hit]}{'old_send'}) and Hit:$hits[$nexthit](start:$blast{$subject}{$hits[$nexthit]}{'old_sstart'}, End:$blast{$subject}{$hits[$nexthit]}{'old_send'}) \nRessigning new coordinates\n\n"
                  if $opt_v;

                # reassign the start and end of the first hit.
                $blast{$subject}{ $hits[$hit] }{'old_send'} = max(
                    $blast{$subject}{ $hits[$hit] }{'old_send'},
                    $blast{$subject}{ $hits[$nexthit] }{'old_send'},
                    $blast{$subject}{ $hits[$hit] }{'old_sstart'},
                    $blast{$subject}{ $hits[$nexthit] }{'old_sstart'}
                );
                $blast{$subject}{ $hits[$hit] }{'old_sstart'} = min(
                    $blast{$subject}{ $hits[$hit] }{'old_send'},
                    $blast{$subject}{ $hits[$nexthit] }{'old_send'},
                    $blast{$subject}{ $hits[$hit] }{'old_sstart'},
                    $blast{$subject}{ $hits[$nexthit] }{'old_sstart'}
                );

                # remove $nexthit as it is assimilated in the first hit. Also delete related hash elements

                delete $blast{$subject}{ $hits[$nexthit] };    #{'old_send'};
                delete $blast{$subject}{ $hits[$nexthit] };    #{'old_sstart'};
                    # if spliced before deleting hash, gives error while deleting hash element, as there is no value $hits[$nexthit] to look for
                my $removed_arry = splice @hits, $nexthit, 1;

# the splice removes the element from array and other elements shift up so in order to check next element, loop need to check same position again. so decrease the value of $nexthit by one;
                $nexthit--;

            }
        }
    }
}
################################################################################################################################

# Print table and Extract hit from sequence file based on blast coordinates
my $hRef_genomic_seq;
if ($opt_s) { $hRef_genomic_seq = ReadFasta($opt_s) }
print "Printing uniq hits\n"                                   if !$opt_o;
print "Extracting hit regions from from genomic sequence\n"    if $opt_s;
print "Saving hit table to file $opt_b.$opt_k.summary.table\n" if $opt_o;
print {$fh_tabout} "Subject\tHit\#\tAln_Start\tAln_End\tAln_length\tSeq_length\n";
my $count = 0;
my %vectrim;
my %done;
foreach my $subject ( keys %blast ) {
    foreach my $hit ( keys %{ $blast{$subject} } ) {

        # print table
        my $length          = $blast{$subject}{$hit}{'old_send'} - $blast{$subject}{$hit}{'old_sstart'} + 1;
        my $final_qcoverage = $length * 100 / $sequence{$subject}{'len'};
        if (
            !$opt_nofilter && $opt_and
            && (   $length < $opt_final_length
                || $final_qcoverage < $opt_final_qcoverage_min
                || $final_qcoverage > $opt_final_qcoverage_max )
          )
        {

            print "AND: At the final filtering stage:
			Skipping this hit as one of these became True
			Length: $length < $opt_final_length
			fqcm: $final_qcoverage < $opt_final_qcoverage_min
			fqcx: $final_qcoverage > $opt_final_qcoverage_max
		\n" if $opt_v;

            next;
        }

        elsif ( !$opt_nofilter && $opt_or
            && !( $length >= $opt_final_length || $final_qcoverage >= $opt_final_qcoverage_min || $final_qcoverage <= $opt_final_qcoverage_max ) )
        {

            print "OR: At the final filtering stage:
			Skipping this hit as none of these are true
			Length: $length < $opt_final_length
			fqcm: $final_qcoverage < $opt_final_qcoverage_min
			fqcx: $final_qcoverage > $opt_final_qcoverage_max
		\n" if $opt_v;

            next;
        }

        print {$fh_tabout} "$subject\t$hit\t$blast{$subject}{$hit}{'old_sstart'}\t$blast{$subject}{$hit}{'old_send'}\t$length\t$sequence{$subject}{'len'}\n";

        $count++;

        # extract sequence

        if ($opt_s) {
            my ($head,$seq);

            # collecting information for vector trim
            push( @{ $vectrim{$subject}{'coords'} }, $blast{$subject}{$hit}{'old_sstart'}, $blast{$subject}{$hit}{'old_send'} );
            my($substr_start_site,$substr_length_extract);
            # calculations for extract sequences
            if ( $blast{$subject}{$hit}{'strand'} eq 'plus' ) {

                $substr_start_site     = ( $blast{$subject}{$hit}{'old_sstart'} - 1 - $opt_h );
                $substr_length_extract = ( $length + $opt_t + $opt_h );

                $substr_start_site = 0
                  if ( $blast{$subject}{$hit}{'old_sstart'} - 1 - $opt_h ) < 0;
                $substr_length_extract = ( length( $$hRef_genomic_seq{$subject} ) - $substr_start_site - 1 )
                  if ( $substr_start_site + $length + $opt_t + $opt_h ) > length( $$hRef_genomic_seq{$subject} );

                
            }
            elsif ( $blast{$subject}{$hit}{'strand'} eq 'minus' ) {

                $substr_start_site     = ( $blast{$subject}{$hit}{'old_sstart'} - 1 - $opt_t );
                $substr_length_extract = ( $length + $opt_t + $opt_h );

                $substr_start_site = 0
                  if ( $blast{$subject}{$hit}{'old_sstart'} - 1 - $opt_t ) < 0;
                $substr_length_extract = ( length( $$hRef_genomic_seq{$subject} ) - $substr_start_site - 1 )
                  if ( $substr_start_site + $length + $opt_t + $opt_h ) > length( $$hRef_genomic_seq{$subject} );

                
            }
            if ($opt_fs) {
                next if $done{$subject};
                $seq = $$hRef_genomic_seq{$subject};
                $head= $subject;
                
                $done{$subject}=1;
            }else{
                $seq = substr( $$hRef_genomic_seq{$subject}, $substr_start_site, $substr_length_extract );
                $seq = revcomp($seq) if $blast{$subject}{$hit}{'strand'} eq 'minus';
                my $seq_length = length($$hRef_genomic_seq{$subject});
                $head="$subject\_Hit:$hit Region:".($substr_start_site+1)." .. ". ($substr_start_site+$substr_length_extract) ." Extracted len:$substr_length_extract strand:$blast{$subject}{$hit}{'strand'} Seq_length:$seq_length";
            }
            
            
            print {$fh_seqout}
              ">$head\n$seq\n" if $seq;
        }

    }

}

# trim vector sequence
if ($vectrim) {
    foreach my $subject ( keys %vectrim ) {

        #push(@{$vectrim{$subject}{'coords'}});
        my $seq                   = $$hRef_genomic_seq{$subject};
        my $num_hits              = @{ $vectrim{$subject}{'coords'} };
        my $seqLength             = length($seq);
        my $substr_start_site     = 0;
        my $substr_length_extract = $seqLength;

        #$seq = revcomp($seq) if $blast{$subject}{$hit}{'strand'} eq 'minus';
        my ( $cleanseq, $clean_seq_start, $clean_seq_end ) = remove_vector( \@{ $vectrim{$subject}{'coords'} }, $seq, $from_end );
        print {$fh_vectrimout} ">$subject Start: $clean_seq_start End: $clean_seq_end length:", $clean_seq_end - $clean_seq_start, "nt\n$cleanseq\n"
          if $cleanseq;
        $vectrim{$subject}{'vectrimmed_seq'} = $cleanseq if $cleanseq;
    }

    ## print rest of the sequences in the vecttrimmed file.
    foreach my $remain_seq ( keys %{$hRef_genomic_seq} ) {
        print {$fh_vectrimout} ">$remain_seq\n$$hRef_genomic_seq{$remain_seq}\n" if !$vectrim{$remain_seq}{'vectrimmed_seq'};

    }

}

print "Total number of uniq hits:$count\n";
print STDOUT "Finished printing uniq hits table in file $opt_b.$opt_k.summary.table\n"
  if ($opt_o);
print STDOUT "Saved sequences in file $outfile" if ( $opt_s && $opt_o );

close($fh_seqout) if ( $opt_s && $opt_o );
close($fh_tabout) if $opt_o;
close(BLAST);
print "\n\n";
exit;

########################################################################################

sub remove_vector {

    my $Ref_coord_array           = shift;
    my $sequence                  = shift;
    my $allowed_distance_from_end = shift;

    my $seq_length = length($sequence);
    $allowed_distance_from_end = $allowed_distance_from_end ? $allowed_distance_from_end : 100;

    my @sorted_coord = sort { $a <=> $b } @$Ref_coord_array;

    # if only one hit region is detected at one end.
    if ( @sorted_coord == 2 ) {
        my ( $start, $length, $terminal ) = get_nonvect_coord( $sorted_coord[0], $sorted_coord[-1], $seq_length, $allowed_distance_from_end );
        return ( substr( $sequence, $start, $length ), $start + 1, $start + 1 + $length, );
    }

    # elsif two or more hit regions were detected at both end or one end. only use first and last hit as vector
    if ( @sorted_coord >= 4 ) {

        my ( $clean_start_1, $clean_length_1, $terminal_1 ) = get_nonvect_coord( $sorted_coord[0], $sorted_coord[1], $seq_length, $allowed_distance_from_end );

        my ( $clean_start_2, $clean_length_2, $terminal_2 ) = get_nonvect_coord( $sorted_coord[-2], $sorted_coord[-1], $seq_length, $allowed_distance_from_end );

        # if both hits are at two different ends extract sequence between first hit end and last hit start. returns clean_seq, start_coord, end_coord
        if ( $terminal_1 ne $terminal_2 ) {
            return ( substr( $sequence, $clean_start_1, $clean_length_2 - $clean_start_1 - 1 ), $sorted_coord[1] + 1, $sorted_coord[-2] - 1 );
        }

        # if both hits are at 5'end extract sequence from last hit end till the end of the sequence. returns clean_seq, start_coord, end_coord
        elsif ( $terminal_1 eq $terminal_2 && $terminal_1 eq 'five' ) {
            return ( substr( $sequence, $clean_start_2, $seq_length - $clean_start_2 - 1 ), $sorted_coord[-1] + 1, $seq_length );
        }

        # if both hits are at 3'end extract sequence from 0 till the first hit start. returns clean_seq, start_coord, end_coord
        elsif ( $terminal_1 eq $terminal_2 && $terminal_1 eq 'three' ) {
            return ( substr( $sequence, 0, $sorted_coord[0] - 1 ), 1, $sorted_coord[0] - 1 );
        }

        # else return whole sequence
        else { return ( $sequence, 1, $seq_length ) }
    }
}

# returns the start, length and terminal info of sequence region which is not part of vector
sub get_nonvect_coord {
    my $start                     = shift;
    my $end                       = shift;
    my $seq_length                = shift;
    my $allowed_distance_from_end = shift;
    $allowed_distance_from_end = $allowed_distance_from_end ? $allowed_distance_from_end : 100;

    # make sure the start is smaller value and end is larger value, in cases they are passed in reverse.
    my $hit_start = min( $start, $end );
    my $hit_end              = max( $start, $end );
    my $five_prime_leftover  = $hit_start;
    my $three_prime_leftover = $seq_length - $hit_end;

    #if hit is at 5' end and is within allowd distance from the start
    if ( $five_prime_leftover < $three_prime_leftover && $allowed_distance_from_end >= $hit_start ) { return ( $hit_end, $three_prime_leftover, 'five' ) }

    #if hit is at 3' end is within allowd distance from the end
    elsif ( $five_prime_leftover > $three_prime_leftover && $allowed_distance_from_end >= $three_prime_leftover ) { return ( 0, $hit_start, 'three' ) }

    # else return full sequence
    else { return ( 0, $seq_length, 'middle' ) }
}

sub min {
    @_ = sort { $a <=> $b } @_;
    return $_[0];

}

sub max {

    @_ = sort { $a <=> $b } @_;
    return $_[-1];

}

sub revcomp {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq;
}

#####subroutine to test overlap
sub does_it_overlap {
    my $start1 = shift;
    my $end1   = shift;
    my $start2 = shift;
    my $end2   = shift;

    if (
        ( ( $start1 <= $start2 ) && ( $end1 >= $end2 ) )    # sequence2 included in sequence 1
        || ( ( $start2 <= $start1 )        && ( $end2 >= $end1 ) )               # sequence1 included in sequence 2
        || ( ( $start1 <= $end2 + $opt_a ) && ( $start1 >= $start2 ) )           # sequence 1 overlaps and is on the right of sequence2.
        || ( ( $end1 <= $end2 )            && ( $end1 + $opt_a >= $start2 ) )    # sequence 1 overlaps and is on the left of sequence2.
        || ( ( $start2 <= $end1 + $opt_a ) && ( $start2 >= $start1 ) )           # sequence 1 overlaps and is on the right of sequence2.
        || ( ( $end2 <= $end1 )            && ( $end2 + $opt_a >= $start1 ) )    # sequence 1 overlaps and is on the right of sequence2.

      )
    {
        return 1;
    }
    else { return 0 }

}

# sub routines for reading fasta file

sub ReadFasta {                                                                  # to read fasta format files into hash. returns hash reference.

    my $seqfile = shift(@_);
    my $demo_header;

    my ( $header, @sequence );
    chomp $seqfile;
    open FASTA, "$seqfile";
    print "Reading Sequences from input file.....Plz wait...\n";
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
        my @scaffold_name = split( /$opt_d/, $header );
        $header = $scaffold_name[ $opt_c - 1 ];
        my $sequence = join( "", @sequence );
        $sequence =~ s/\s//g;
        $sequence =~ s/\n//g;

        if ( $header =~ /^\s*$/ ) { next; }
        $demo_header = $header;

        # Make sure no duplicate header exists else they will be overwritten and lost during hashing due to same key.
        if ( !exists $seq_hash{$header} ) {
            $seq_hash{$header} = $sequence;    #feed headers and sequences in hash.
                                               #$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
        }
        else {

            # find a uniq header name by adding a number at the end. If header still exists, increase the number by one
            while ( exists $seq_hash{$header} ) {
                $header = $header . $last_N;
                $last_N++;
            }

            $seq_hash{$header} = $sequence;

            #$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;

        }
    }

    my @seq_count = keys(%seq_hash);
    my $seq_count = @seq_count;

    print "Done....\nNumber of sequences read form input file = $seq_count\n\nExample of one the headers in sequence\n$demo_header\n";
    @seq_count = ();
    $/         = "\n";    # Record seperator set back to default newline.

    return ( \%seq_hash );

}

#-------------------------------------End ReadFasta---------------------------------------+
