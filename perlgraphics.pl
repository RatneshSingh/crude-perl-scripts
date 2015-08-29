#! /usr/bin/perl -w
use strict;
use GD::Simple;
use Getopt::Long;
our ( $help, @seq_lenths, %namesofseqs, $fil_name, $show_coords, $col_numHit );

# Declare default values for options
my $maxbarlen           = 4000;
my $maxbarwidth         = 150;
my $marginx             = 100;
my $marginy             = 100;
my $distancebetweenbars = 100;
my $domainraiseup       = 0;
my $domainraisedown     = 0;
my $maxseqlength        = 5000000;
my $rect_bg             = 'white';
my $rect_fg             = 'black';
my $feature_bg          = 'red';
my $feature_fg          = undef;
my $font                = 'Times:italic';
my $fontsize            = 50;
my $angle               = 0;
my $inputfile;
my $format       = 'BlastTable';
my $position     = 20;
my $col_name     = 1;
my $col_start    = 7;
my $col_end      = 8;
my $col_slen     = 15;
my $col_perid    = 3;
my $col_bitscore = 12;
my $col_alnlen   = 4;
my $col_qcov     = 13;
my $col_scov     = 14;
my $scale_min    = 1000;
my $scale_max    = 10000;
my $gradient     = 'none';
my $fil_aln_len  = 0;
my $per_idy_min  = 0;
my $per_idy_max  = 100;
my $fil_qcov     = 0;
my $fil_scov     = 0;
my $fil_numhit   = 1;
my $draw         = 'query';
my $size         = 'medium';
my $scalebar_pos = 50;
my $legend_size  = 50;
my $ps           = 5;              # pensize
my $max_numhit   = 0;
$distancebetweenbars += $position;
my $distancefromtop = $marginy + $fontsize + $position + $scalebar_pos + $legend_size;

# get optins from comand line
my $result = GetOptions(
    "table|input=s"               => \$inputfile,
    "format=s"                    => \$format,
    "maxbarlen|mbl=i"             => \$maxbarlen,
    "maxbarwidth|mbw=i"           => \$maxbarwidth,
    "marginx|mx=i"                => \$marginx,
    "marginy|my=i"                => \$marginy,
    "distance_between_bars|dbb=i" => \$distancebetweenbars,
    "domainraiseup|dru=i"         => \$domainraiseup,
    "domainraisedown|drd=i"       => \$domainraisedown,
    "maxseqlength|msl=i"          => \$maxseqlength,
    "distancefromtop|dft=i"       => \$distancefromtop,
    "rect_bg|rbg=s"               => \$rect_bg,
    "rect_fg|rfg=s"               => \$rect_fg,
    "feature_bg|fbg=s"            => \$feature_bg,
    "feature_fg|ffg=s"            => \$feature_fg,
    "font|f=s"                    => \$font,
    "fontsize|fs=i"               => \$fontsize,
    "angle|a=i"                   => \$angle,
    "position|p=i"                => \$position,
    "col_name|cn=i"               => \$col_name,
    "col_start|cs=i"              => \$col_start,
    "col_end|ce=i"                => \$col_end,
    "col_alnlen|cal=i"            => \$col_alnlen,
    "col_bitscore|cbs=i"          => \$col_bitscore,
    "col_slen|csl=i"              => \$col_slen,
    "col_perid|cp=i"              => \$col_perid,
    "col_qcov|cqc=i"              => \$col_qcov,
    "col_scov|csc=i"              => \$col_scov,
    "scale_min|sm=i"              => \$scale_min,
    "scale_max|sx=i"              => \$scale_max,
    "gradient|gd=s"               => \$gradient,
    "aln_len|al=i"                => \$fil_aln_len,
    "query_coverage|qcov|qc=i"    => \$fil_qcov,
    "subject_coverage|scov|sc=i"  => \$fil_scov,
    "fil_numhit|fnh=i"            => \$fil_numhit,
    "per_idy_max|pix=i"           => \$per_idy_max,
    "per_idy_min|pim=i"           => \$per_idy_min,
    "draw|d=s"                    => \$draw,
    "show_coords"                 => \$show_coords,
    "fil_name=s"                  => \$fil_name,
    "help|h"                      => \$help,
    "size=s"                      => \$size
);

# define usage for help
my $usage = "

usage: perl scriptname -table input_table_file -options


	-table|input|i	Input table containing coordinates

Input Table options [optional]
	-format	table format. e.g. BlastTable,Blast,FilteredTable [BlastTable]
	-col_name|cn	column of table containing name of the sequence[1]
	-col_start|cs	column of table containing start of the feature[7]
	-col_end|ce	column of table containing end of the feature[8]
	-col_perid|cp	Column containing percent identity of the alignment[3]
	-col_slen|csl	Column containing sequence length[15]
	-col_alnlen|cal	Column containing alignment length[4]
	-col_bitscore|cbs	Column containing bit score[12]
	-col_scov|csc	Column containing query coverage [15]
	-col_qcov|cqc	Column containing subject coverage[16]

Filter option [optional]
	-per_idy_min|pim	Exclude data with percent identity of alignment less than this [0]
	-per_idy_max|pix	Exclude data with percent identity of alignment more than this [100]
	-aln_len|al	Exclude data with alignment length smaller than this[0]
	-query_coverage|qcov|qc	Filter hits where query coverage (% of length) is less than this[0]
	-subject_coverage|scov|sc	Filter hits where subject coverage(% of length) is less than this[0]
	-fil_name	Only draw graphics for the sequence names matching this name
	-fil_numhit|fnh	Only draw the sequences which have this many hits [1].
					Counts uniq hits when used with filtered summdup NR table. Count all the hits
					when used with normal blast table. use large values for normal table as hit
					numbers will be when multiple queries hit same region of the subject.

Sequence [optional]
	-maxbarlen|mbl	length (in pixels) of the bar representing sequence[8000]
	-maxbarwidth|mbw	width (in pixels) of the bar representing sequence[150]
	-marginx|mx	Margin on x axis [100]
	-marginy|my	Margin on y axis [100]
	-distance_between_bars|dbb	Distance between each sequence in vertical position[100]
	-maxseqlength|msl	Max length of the sequence when this info is absent from blast table [100000]
	-distancefromtop|dft	Distance of the first sequence from the top [200]
	-rect_bg|rbg	Background color of the rectangle representing sequence [white]
	-rect_fg|rfg	Outline color of the rectangle representing sequence [black]
	-gradient|gd	bitscore|perid. Use color gradient based on bitscore/percent identity of the alignment[none]
	-draw|d	query|subject	Draw either query or subject.[query]

Features [optional]
	-domainraiseup|dru	Print the hit/domain raised above the bar [0]
	-domainraisedown|drd	Print the hit/domain down the bar [0]
	-feature_bg|fbg	Background color of the rectangle representing feature [red]
	-feature_fg|ffg	Outline color of the rectangle representing sequence [none]
	-show_coords	Print the coordinate of the hit on the figure [No]
	-size			Size of the figure drawn. verysmall, small, medium, large [medium]

Heading for sequence [optional]
	-font|f	Heading font in font:style format [Times:italic]
	-fontsize|fs	Size of the font [20]
	-position|p	Position of the heading wrt to sequence. positive to raise above the sequence or negative values to drop below the sequence start position. [10]
	-angle|a	Angle of the heading text [0]
	-help|h	Terminate the program and show help file.

Scale [optional]
	-scale_min|sm	Minor marks on the scale [10000].
	-scale_max|sx	Major marks on the scale [1000].
";

# die if input file is not provided or user invoked help tag
die "\n$usage\n"                                     if $help;
die "\nError:No input file detected....\n\n$usage\n" if !$inputfile;

# change the values if user provides blasttable and asked to draw subject instead of query
$gradient =~ s/\s+//g;
$draw     =~ s/\s+//g;
if ( lc $draw eq 'subject' && lc $format eq 'blasttable' ) {

    $col_name  = 2;
    $col_start = 9;
    $col_end   = 10;
    $col_slen  = 16;
}
elsif ( lc $format eq 'filteredtable' ) {
    $col_name   = 1;
    $col_numHit = 2;
    $col_start  = 3;
    $col_end    = 4;
    $col_alnlen = 5;
    $col_slen   = 6;

}

# define parameters based on size of the graph requested. Larger figures are harder to open on pc.
if ( $size eq 'verysmall' ) {

    $maxbarlen           = 1000;
    $maxbarwidth         = 10;
    $marginx             = 10;
    $marginy             = 10;
    $distancebetweenbars = 12;
    $domainraiseup       = 0;
    $domainraisedown     = 0;
    $distancefromtop     = 80;
    $font                = 'Times:italic';
    $fontsize            = 10;
    $position            = 1;
    $scalebar_pos        = 20;
    $legend_size         = 20;
    $ps                  = 1;
    $distancebetweenbars += $position;
    $distancefromtop = $marginy + $fontsize + $position + $scalebar_pos + $legend_size;

}
elsif ( $size eq 'small' ) {

    $maxbarlen           = 2000;
    $maxbarwidth         = 20;
    $marginx             = 20;
    $marginy             = 20;
    $distancebetweenbars = 30;
    $domainraiseup       = 0;
    $domainraisedown     = 0;
    $distancefromtop     = 100;
    $font                = 'Times:italic';
    $fontsize            = 25;
    $position            = 5;
    $scalebar_pos        = 40;
    $legend_size         = 40;
    $ps                  = 3;
    $distancebetweenbars += $position;
    $distancefromtop = $marginy + $fontsize + $position + $scalebar_pos + $legend_size;

}
elsif ( $size eq 'large' ) {

    $maxbarlen           = 8000;
    $maxbarwidth         = 200;
    $marginx             = 200;
    $marginy             = 200;
    $distancebetweenbars = 200;
    $domainraiseup       = 0;
    $domainraisedown     = 0;
    $font                = 'Times:italic';
    $fontsize            = 75;
    $position            = 50;
    $scalebar_pos        = 200;
    $legend_size         = 200;
    $ps                  = 8;
    $distancebetweenbars += $position;
    $distancefromtop = $marginy + $fontsize + $position + $scalebar_pos + $legend_size;
}
###########################################################################################################################################
# collect the names of the sequence and length info to draw the bar.
open FILE, "$inputfile";
while (<FILE>) {
    next if $_ =~ /^\s*$/;
    next
      if ( lc $format eq 'filteredtable'
        && ( $_ =~ /#/gi || $_ =~ /Aln_Start/i ) );
    my @line = split(/\s+/);

    #############################################
    # which column contains the information name , start and end of the feature
    my $nameseq       = $line[ $col_name - 1 ];
    my $start         = $line[ $col_start - 1 ];
    my $end           = $line[ $col_end - 1 ];
    my $seq_len       = $line[ $col_slen - 1 ];
    my $aln_per_id    = $line[ $col_perid - 1 ];
    my $aln_bitscore  = $line[ $col_bitscore - 1 ];
    my $alignment_len = $line[ $col_alnlen - 1 ];
    my $qcov          = $line[ $col_qcov - 1 ];
    my $scov          = $line[ $col_scov - 1 ];
    my $numhit        = $line[ $col_numHit - 1 ] if ( lc $format eq 'filteredtable' );
    ###################################################
# assign some values if table does not contain all the columns required for this script other wise hits wont pass filter and give error
    if ( lc $format eq 'filteredtable' ) {
        $aln_per_id   = 100;
        $aln_bitscore = 1000;
        $qcov         = 100;
        $scov         = 100;
    }

    # filter results based on user provided criteria
    next if ( $aln_per_id > $per_idy_max || $aln_per_id < $per_idy_min );
    next if ( $alignment_len < $fil_aln_len );
    next if ( $qcov < $fil_qcov          || $scov < $fil_scov );
    next if ( $fil_name && $nameseq !~ /$fil_name/gi );
    #####################################################

    ######################################################
    # calculate num_hit buy counting the name of subject or query if filtered table is not used
    # the count the occurence of $nameseq and use it as numHit
    if ( lc $format ne 'filteredtable' ) {
        $namesofseqs{$nameseq}{'numHit'}++;
        $numhit = 0;
        $max_numhit = $namesofseqs{$nameseq}{'numHit'} if $max_numhit < $namesofseqs{$nameseq}{'numHit'};

        #$numhit=$namesofseqs{$nameseq}{'numHit'};

    }



    chomp($nameseq);
    $nameseq =~ s/\s+//g;
    if ( $nameseq eq '' ) { next; }
    else {
        $namesofseqs{$nameseq}{'name'}   = $nameseq;
        $namesofseqs{$nameseq}{'slen'}   = $seq_len;
        $namesofseqs{$nameseq}{'numHit'} = $numhit
          if ( !$namesofseqs{$nameseq}{'numHit'} || $namesofseqs{$nameseq}{'numHit'} < $numhit );

        #$namesofseqs{$nameseq}{'length'} = $seq_len;
        push( @seq_lenths, $seq_len );
    }
}

# Adjust length and number of bar to draw if filtering by number of hits on th esequence
my $barcount = 0;
if ( lc $format eq 'filteredtable' ) {
    @seq_lenths = ();
    foreach my $key ( keys %namesofseqs ) {
        $barcount++ if ( $namesofseqs{$key}{'numHit'} >= $fil_numhit );
        push( @seq_lenths, $namesofseqs{$key}{'slen'} ) if ( $namesofseqs{$key}{'numHit'} >= $fil_numhit );

    }
}

# die if no hit passed the filter
die "No hits passed the filtering criteria. Use less stringent filter"
  if @seq_lenths <= 0;
close FILE;

# calculate the width and length of graphic window
my $pixelperbase = $maxbarlen / max(@seq_lenths);
my $num_of_bars  = scalar( keys %namesofseqs );
$num_of_bars = $barcount if lc $format eq 'filteredtable';
my $graphics_width =
  ( $num_of_bars * $maxbarwidth ) +
  ( ( $num_of_bars - 1 ) * $distancebetweenbars ) +
  ( $distancefromtop * 2 ) +
  ( $marginy * 2 );
my $graphic_length = $maxbarlen + ( $marginx * 2 );

# create a new image for calculated length and width
my $figure = GD::Simple->new( $graphic_length, $graphics_width );
##############################################################################################
# Draw legend on the top of the figure
$figure->penSize( $ps, $ps );
my $space = int( $legend_size / 5 );

if ( $gradient eq "perid" ) {

    #"red","pink","green","lime","blue","deepskyblue","aqua","yellow","black";
    #print "Gradient is :$gradient\n";
    my $one_block = $maxbarlen / 10;
    my %color     = (
        '0' => 'red',
        '1' => 'pink',
        '2' => 'green',
        '3' => 'lime',
        '4' => 'blue',
        '5' => 'skyblue',
        '6' => 'aqua',
        '7' => 'yellow',
        '8' => 'black'
    );
    my %name = (
        '0' => '100-90%',
        '1' => '90-80%',
        '2' => '80-70%',
        '3' => '70-60%',
        '4' => '60-50%',
        '5' => '50-40%',
        '6' => '40-30%',
        '7' => '30-20%',
        '8' => '20-10%'
    );
    for ( my $n = 0; $n <= 8; $n++ ) {
        $figure->bgcolor( $color{$n} );
        $figure->rectangle(
            $marginx + ( $n * $one_block ),
            $distancefromtop - $legend_size - $scalebar_pos - $fontsize,
            $marginx + ( $n * $one_block ) + $legend_size,
            $distancefromtop - $scalebar_pos - 1.5 * $fontsize
        );
        $figure->moveTo(
            $marginx + $legend_size + ( $n * $one_block ) + $space,
            $distancefromtop - $scalebar_pos - 1.5 * $fontsize
        );
        $figure->font('Times:italic');
        $figure->fontsize($fontsize);
        $figure->string( $name{$n} );
    }

}
elsif ( $gradient eq "bitscore" ) {

    #print "Gradient is :$gradient\n";
    my $one_block = $maxbarlen / 5;
    my %color     = (
        '0' => 'red',
        '1' => 'fuchsia',
        '2' => 'lime',
        '3' => 'blue',
        '4' => 'black'
    );
    my %name = (
        '0' => '>200',
        '1' => '80-200%',
        '2' => '50-80',
        '3' => '40-50',
        '4' => '<40'
    );
    for ( my $n = 0; $n <= 4; $n++ ) {
        $figure->bgcolor( $color{$n} );
        $figure->rectangle(
            $marginx + ( $n * $one_block ),
            $distancefromtop - $legend_size - $scalebar_pos - 1.5 * $fontsize,
            $marginx + ( $n * $one_block ) + $legend_size,
            $distancefromtop - $scalebar_pos - 1.5 * $fontsize
        );
        $figure->moveTo(
            $marginx + $legend_size + ( $n * $one_block ) + $space,
            $distancefromtop - $scalebar_pos - 1.5 * $fontsize
        );
        $figure->font('Times:italic');
        $figure->fontsize($fontsize);
        $figure->string( $name{$n} );
    }
}

else {
    print "Gradient is :$gradient\n";
    my $n         = 0;
    my $one_block = $maxbarlen / 5;
    $figure->bgcolor('red');
    $figure->rectangle(
        $marginx + ( $n * $one_block ),
        $distancefromtop - $legend_size - $scalebar_pos - 1.5 * $fontsize,
        $marginx + ( $n * $one_block ) + $legend_size,
        $distancefromtop - $scalebar_pos - 1.5 * $fontsize
    );
    $figure->moveTo( $marginx + $legend_size + ( $n * $one_block ) + $space,
        $distancefromtop - $scalebar_pos - 1.5 * $fontsize );
    $figure->font('Times:italic');
    $figure->fontsize($fontsize);
    $figure->string("Blast Hit");
}

###########################################################################################################
# Draw a scale bar using values from maxbarlen
$figure->penSize( $ps, $ps );
$figure->moveTo( $marginx, $distancefromtop - $scalebar_pos );
$figure->lineTo( $maxbarlen + $marginx, $distancefromtop - $scalebar_pos );

$maxseqlength = max(@seq_lenths) if scalar @seq_lenths >= 1;
$scale_max    = 100              if ( $maxseqlength > 100 );
$scale_max    = 1000             if ( $maxseqlength > 1000 );
$scale_max    = 10000            if ( $maxseqlength > 10000 );
$scale_max    = 100000           if ( $maxseqlength > 200000 );
$scale_max    = 1000000          if ( $maxseqlength > 2000000 );
$scale_max    = 10000000         if ( $maxseqlength > 20000000 );

#$scale_max=round_nearest($maxseqlength)/10 if scalar@seq_lenths >=1;
#$scale_min=round_nearest($maxseqlength)/100 if scalar@seq_lenths >=1;
$scale_min = $scale_max / 10;
for ( my $i = 0; $i <= $maxseqlength; $i = $i += $scale_min ) {

    $figure->moveTo( $marginx + int( $i * $pixelperbase ), $distancefromtop - $scalebar_pos );
    $figure->lineTo( $marginx + int( $i * $pixelperbase ), $distancefromtop - int( $fontsize / 2 ) )
      if ( ( $i % $scale_max ) == 0 );
    $figure->lineTo( $marginx + int( $i * $pixelperbase ), $distancefromtop - $fontsize )
      if ( ( $i % $scale_max ) != 0 );
    $figure->moveTo( $marginx + int( $i * $pixelperbase ), $distancefromtop - $scalebar_pos - $position );
    $figure->font('Times:italic');
    $figure->fontsize($fontsize);
    $figure->string($i) if ( ( $i % $scale_max ) == 0 );

}

###########################################################################################################
# Draw rectangles representing sequence
$figure->penSize( 1, 1 );

foreach my $value ( keys %namesofseqs ) {
    next
      if ( lc $format eq 'filteredtable'
        && $namesofseqs{$value}{'numHit'}
        && $namesofseqs{$value}{'numHit'} < $fil_numhit );

    print "\nSequence Name: $value\n";

    # Draw one bar for each sequence name
    $figure->bgcolor($rect_bg);
    $figure->fgcolor($rect_fg);

    $maxbarlen = $pixelperbase * $namesofseqs{$value}{'slen'};

    $figure->rectangle(
        $marginx,
        $distancefromtop + $marginy,
        $maxbarlen + $marginx,
        $marginy + $maxbarwidth + $distancefromtop
    );

    # save the relative position of each bar on y axis (from up to down)
    $namesofseqs{$value}{'y1'} = $marginy + $distancefromtop - $domainraiseup;
    $namesofseqs{$value}{'y2'} = $marginy + $maxbarwidth + $distancefromtop + $domainraisedown;

    # print the name of the sequence in the top corner
    $figure->moveTo( $marginx, $distancefromtop + $marginy - $position );
    $figure->font($font);
    $figure->fontsize($fontsize);
    $figure->angle($angle);
    $figure->string("$value        $namesofseqs{$value}{'slen'}bp");

    #print "Position for y1:",$namesofseqs{$value}{'y1'}."\n";
    #print "Position for y2:",$namesofseqs{$value}{'y2'}."\n";

    $distancefromtop = $distancefromtop + $distancebetweenbars + $maxbarwidth;

}
#############################################################################################################################
# Draw the rectangles for each hit
open FILE, "$inputfile";
while (<FILE>) {
    next if $_ =~ /^\s*$/;
    next
      if ( lc $format eq 'filteredtable'
        && ( $_ =~ /#/gi || $_ =~ /Aln_Start/i ) );
    my @line = split(/\s+/);

    s/\s+//g foreach @line;

    # which column contains the information name , start and end of the feature
    my $nameseq       = $line[ $col_name - 1 ];
    my $start         = $line[ $col_start - 1 ];
    my $end           = $line[ $col_end - 1 ];
    my $aln_per_id    = $line[ $col_perid - 1 ];
    my $aln_bitscore  = $line[ $col_bitscore - 1 ];
    my $alignment_len = $line[ $col_alnlen - 1 ];
    my $qcov          = $line[ $col_qcov - 1 ];
    my $scov          = $line[ $col_scov - 1 ];

    #my$numhit			=$line[$col_numHit-1];# if ( lc $format eq 'filteredtable' );
    ###################################################
# assign some values if table does not contain all the columsn required for this script other wise hits wont pass filter and give error
    if ( lc $format eq 'filteredtable' ) {
        $aln_per_id   = 100;
        $aln_bitscore = 1000;
        $qcov         = 100;
        $scov         = 100;
    }

    # filter results based on user provided criteria
    next if ( $fil_name && $nameseq !~ /$fil_name/gi );
    next if ( $aln_per_id > $per_idy_max || $aln_per_id < $per_idy_min );
    next if ( $alignment_len < $fil_aln_len );
    next if ( $qcov < $fil_qcov          || $scov < $fil_scov );
    next
      if ( lc $format eq 'filteredtable'
        && $namesofseqs{$nameseq}{'numHit'}
        && $namesofseqs{$nameseq}{'numHit'} < $fil_numhit );

    #print "AlnBitScore:$aln_bitscore\n";

#print "Drawing....Seq Aln_length:$alignment_len\tUserAlnLen:$fil_aln_len\t\tSeq Per_idy:$aln_per_id\tUser per idy Max:$per_idy_max\tMin:$per_idy_min\n";
    chomp($nameseq);
    $start =~ s/\D+//g;
    $end   =~ s/\D+//g;
    if ( $start !~ /\d+/ ) {
        print "start:$start is not numeric. skipping\n";
        next;
    }
    if ( $end !~ /\d+/ ) { print "end:$end is not numeric. skipping\n"; next }

    my $positionx1 = int( ( min( $start, $end ) * $pixelperbase ) ) + $marginx;
    my $positiony1 = $namesofseqs{$nameseq}{'y1'};                                #7;
    my $positionx2 = int( ( max( $start, $end ) * $pixelperbase ) ) + $marginx;
    my $positiony2 = $namesofseqs{$nameseq}{'y2'};                                #53;
                                                                                  #$seqname=$namesofseqs{$name};
                                                                                  #draw red rectangles on the gray bars
    $figure->bgcolor($feature_bg);
    $figure->bgcolor( get_color_perid( int($aln_per_id) ) )
      if $gradient eq 'perid';
    $figure->bgcolor( get_color_bitscore( int($aln_bitscore) ) )
      if $gradient eq 'bitscore';
    $figure->fgcolor($feature_fg);
    $figure->rectangle( $positionx1, $positiony1, $positionx2, $positiony2 );

    # if asked print the coordinates
    my $fg = my $bg = 'red';
    $fg = $bg = get_color_perid( int($aln_per_id) ) if $gradient eq 'perid';
    $fg = $bg = get_color_bitscore( int($aln_bitscore) )
      if $gradient eq 'bitscore';

    if ($show_coords) {
        print_label( $positionx1, $positiony1, min( $start, $end ), $bg, $fg, \$figure );
        print_label( $positionx2, $positiony1, max( $start, $end ), $bg, $fg, \$figure );
    }
}

# print the image out as jpeg and png files.
open OUT2, ">$inputfile.png";
print OUT2 $figure->png;
close OUT2;

#open OUT2, ">$inputfile.jpg";
#print OUT2 $figure->jpeg;
#close OUT2;

##########################################################

sub min {
    @_ = sort { $a <=> $b } @_;
    return $_[0];

}

sub max {

    @_ = sort { $a <=> $b } @_;
    return $_[-1];

}

sub get_color_perid {

    my $per_id = shift;
    $per_id =~ s/\s+//g;
    if    ( $per_id > 90 ) { return ('red') }
    elsif ( $per_id > 80 ) { return ('pink') }
    elsif ( $per_id > 70 ) { return ('green') }
    elsif ( $per_id > 60 ) { return ('lime') }
    elsif ( $per_id > 50 ) { return ('blue') }
    elsif ( $per_id > 40 ) { return ('skyblue') }
    elsif ( $per_id > 30 ) { return ('aqua') }
    elsif ( $per_id > 20 ) { return ('yellow') }
    else                   { return ('black') }
}

sub get_color_bitscore {

    my $bitscore = shift;
    $bitscore =~ s/\s//g;
    if    ( $bitscore >= 200 ) { return ('red') }
    elsif ( $bitscore > 80 )   { return ('fuchsia') }
    elsif ( $bitscore > 50 )   { return ('lime') }
    elsif ( $bitscore > 40 )   { return ('blue') }
    else                       { return ('black') }
}

sub print_label {

    my $pos_x         = shift;
    my $pos_y         = shift;
    my $text          = shift;
    my $bg_color      = shift;
    my $fg_color      = shift;
    my $ref_to_figure = shift;
    my $pos_y2        = $pos_y - 5;

    #print "\nRecieved Posx:$pos_x, posy:$pos_y, Text: $text\n";

    $bg_color = 'green' if !defined $bg_color;
    $fg_color = 'green' if !defined $fg_color;

    # Draw line to the feature start
    $$ref_to_figure->fgcolor($fg_color);
    $$ref_to_figure->bgcolor($bg_color);
    $$ref_to_figure->moveTo( $pos_x, $pos_y );
    $$ref_to_figure->lineTo( $pos_x, $pos_y2 );

    # Write text as label
    $$ref_to_figure->moveTo( $pos_x, $pos_y2 );
    $$ref_to_figure->angle($angle);          #Global angle
    $$ref_to_figure->font($font);            # Global $font
    $$ref_to_figure->fontsize($fontsize);    # Global $fontsize
    $$ref_to_figure->string("$text");
}

sub round_nearest {
    my $num = shift;
    return ( '1' . ( '0' x length( int($num) ) ) );

}
