#!/usr/bin/perl
use warnings;
use strict;
use GD;
use GD::Simple;
use Getopt::Long;


our($dagchainer_table,$outfile,$infile,%all_info,%lg_info,%contig_info,%cumulative_contig,%cumulative_lg,$verbose,$help,$lg_top,$cont_top,$svg,$order_by_homology,$referencex,$referencey,$refseqlenx,$refseqleny,$lgtable,$as_marker,$blast_table,$reg_table,$seq1,$loc1,$seq2,$loc2,$switch_axis,$fil_alen,$fil_eval,$fil_pid,$fill_name,$fil_qname,$fil_sname,$exact,$keep_same_order,$order_by_size,$order_by_alpha,$order_by_num,$print_list,$no_lines,$same_scale,$order_y,$order_x,$zoom_on_hits,$equal_block);
our($rot_x_lab,$rot_y_lab,$rot_xy_lab,$scale_x,$scale_y,$scale_x_dist,$scale_y_dist,$no_file_name);

### initialize some values
my$img_width||=1200;
my$img_height||=600;
my$margin=100;
my$font_size=10;
my$dotsize=3;
my$min_markers=1;
my$resolution=10000;
$rot_x_lab = -90;
$rot_y_lab = 0;
GetOptions(
    "out|o=s"=>\$outfile,
    "infile|input|in|i=s"=>\$infile,
    "lgtable|lgt=s"=>\$lgtable,
    "dagchainer_table|dt=s"=>\$dagchainer_table,
    ### Regular table related parameters
    "reg_table|rt=s"=>\$reg_table,
    "seq1|s1=i"=>\$seq1,
    "loc1|l1=i"=>\$loc1,
    "seq2|s2=i"=>\$seq2,
    "loc2|l2=i"=>\$loc2,

    ### blast file related parameters
    "blast_table|blasttable|blast|bt=s"=>\$blast_table,
    "fil_alen|fal=i"=>\$fil_alen,
    "fil_eval|fev=f"=>\$fil_eval,
    "fil_pid|fpi=f"=>\$fil_pid,
    "fil_name|fnm=s"=>\$fill_name,
    "fil_qname|fqn=s"=>\$fil_qname,
    "fil_sname|fsn=s"=>\$fil_sname,
    "exact"=>\$exact,
    "resolution|res=i"=>\$resolution,
    "keep_same_order|kso"=>\$keep_same_order,
    "keep_same_scale|ksc"=>\$same_scale,
    "rotate_x_label|rxl=i"=>\$rot_x_lab,
    "rotate_y_label|ryl=i"=>\$rot_y_lab,
    "rotate_xy_label|rxyl|ryxl=i"=>\$rot_xy_lab,
    "order_by_size|obs"=>\$order_by_size,
    "order_by_alpha|oba"=>\$order_by_alpha,
    "order_by_number|obn"=>\$order_by_num,
    "order_by_homology|olg|obh"=>\$order_by_homology,
    "user_defined_order_x|udox=s"=>\$order_x,
    "user_defined_order_y|udoy=s"=>\$order_y,
    "img_width|iw=i"=>\$img_width,
    "img_height|ih=i"=>\$img_height,
    "lg_top|lt=i" =>\$lg_top,
    "cont_top|ct=i"=>\$cont_top,
    "switch_axis|swa"=>\$switch_axis,
    "no_grid_lines|ngl"=>\$no_lines,
    "margin|m=i"=>\$margin,
    "font_size|fs=i"=>\$font_size,
    "dotsize|ds=i"=>\$dotsize,
    "scale_x|scx"=>\$scale_x,
    "scale_x_dist|dscx=i"=>\$scale_x_dist,
    "scale_y|scy"=>\$scale_y,
    "scale_y_dist|dscy=i"=>\$scale_y_dist,
    "no_file_name|nfn"=>\$no_file_name,
    "verbose|v"=>\$verbose,
    "svg"=>\$svg,
    "print_list|pl"=>\$print_list,
    "zoom_on_hits|zoh|zh|z"=>\$zoom_on_hits,
    "equal_block|eb"=>\$equal_block,
    "referencex|rfx=s"=>\$referencex,
    "referencey|rfy=s"=>\$referencey,
    "refseqlenx|rflx|rslx=s"=>\$refseqlenx,
    "refseqleny|rfly|rsly=s"=>\$refseqleny,
    "min_markers|mm=i"=>\$min_markers,
    "as_marker|am=s"=>\$as_marker,
    "help|h"=>\$help,
);

## change module to svg if requested
GD::Simple->class('GD::SVG') if $svg;

my$usage="
$0 (infile|lgtable|blastTable) options..
Options:
    -infile|i           infile from mstmap result out
    -lgtable|lgt        lg table with 4 columns.If -lgt is not provided with -i, script assumes that marker names in mst output file includes this info as \"MarkerName_Contig_number_location\"
                        eg. MarkerName\tContigName\tLocationOnContig\tLinkageGroup\tlocationOnLG
    -dagchainer_table|dt result produced by dagchainer.
    ********************************************************************
    -blast_table|bt     blast result in table format.
    -fil_alen|fal       only plot hits of this alignment length or higher
    -fil_eval|fev       only plot hits of this evalue or lower
    -fil_pid|fpi        only plot hits of this percent identity or higher

    -fill_name|fnm	only hits with this string as part of name in query or subject.
    -fil_qname|fqn	only hits with this string as part of name in query.
    -fil_sname|fqs	only hits with this string as part of name in subject.
    -exact		select exact matches instead of pattern match

    -resolution|res     resolution thrshold in nt to draw multiple dots if alignment is longer than resolution[10000]
    -keep_same_order|kso Keep the order of contigs on y and x axis same. Use it only when using blast result to self.
                        This opton will keep the broken diagonal to minimum when dot plot to self is plotted.
    -keep_same_scale|ksc  Keep scale of drawing on both axis same. Usefull when both axis are drawing sequences
                        and user want to keep the aspect ratios of lengths.
    ********************************************************************
    -reg_table|rt       Regular table with 4col of info,
    -seq1               column containing Sequence1 Name  [1]
    -loc1               column Sequence1   Match Location [2]
    -seq2               column containing Sequence2 Name  [3]
    -loc2               column Sequence2   Match Location [4]
    ********************************************************************
    -referencex|rfx     reference genome fasta file to get size for x axis.
    -referencey|rfy     reference genome fasta file to get size for y axis.
    -refseqlenx|rflx    file for x axis containing reference sequence length as seqName\tLength.
    -refseqleny|rfly    file for y axis containing reference sequence length as seqName\tLength.
    ********************************************************************
    -out|o              outfile [infile.png ]
    -print_list|pl      print the list of names showing homology
    -img_width|iw       img_width [ 1200 ]
    -img_height|ih      img_height [ 600 ]
    -lg_top|lt          lg_top [ 100 ]
    -cont_top|ct        cont_top [ 100 ]
    -margin|m           margin [ 100 ]
    -font_size|fs       font_size [ 10 ]
    -dotsize|ds         dotsize [ 3 ]
    -no_grid_lines|ngl  no_grid_lines
    -zoom_on_hits|zh    zoom plotting on the hit regions only.
    -equal_block|eb     assigne equal block for each query/subject irrespective of their size.
                        each entry will not on scale to each other.
    -verbose|v          verbose
    -svg                print image as svg. could not make it to print text.Only dotplot.

    ## order names on axis [default ordering is by size]
    -order_by_size|obs  order_by_size [default],
    -order_by_alpha|oba  order_by_alpha,
    -order_by_number|obn  order_by_num. order only by numbers in the names.
    -order_by_homology|olg  order_by_homology,
    -user_defined_order_x|udox user defined order on X axis from a file
    -user_defined_order_y|udoy user defined order on Y axis from a file
    -rotate_x_label|rxl rotate x axis labels [-90]
    -rotate_y_label|ryl rotate y axis labels [0]
    -rotate_xy_label|rxyl rotate xy axis labels
    -scale_x|scx  Add scale on x-axis (cannot be used with zoomed hits option)
    -scale_x_dist|dscx use provided scale for x axis
    -scale_y|scy  Add scale on y-axis (cannot be used with zoomed hits option)
    -scale_y_dist|dscy use provided scale for y axis
    -switch_axis|swa    switch_axis
    -order_by_lg|olg    order contigs based on similarity to lgs.
    -min_markers|mm     Minimum numbe of markers to present for a contig to be aligned on LG
    -as_marker|am       Use user provided char as marker for dotplot.eg. x
    -no_file_name|nfn   Dont print the outfilename at the bottom of  figure.
    -help|h             help

";

#### read linkage file
($infile || $lgtable || $blast_table || $reg_table || $dagchainer_table) or die "\n Please provide input file \n$usage\n";
$help && die "\n Please provide MSTMap linkage file \n$usage\n";



$outfile||=join "",$infile,
$order_by_homology?".lg_ordered":"",
$zoom_on_hits?".zoomedHits":"",
".png" if $infile;

$outfile||=join "",$lgtable,
$order_by_homology?".lg_ordered":"",
$zoom_on_hits?".zoomedHits":"",
".png" if $lgtable;

$outfile||=join "",$blast_table,
$order_by_homology ? ".Homology_ordered" : "",
$fil_alen ? ".alen$fil_alen" : "",
$fil_pid ? ".pid$fil_pid" : "",
$fil_eval ? ".eval$fil_eval" : "",
$fill_name ? ".Name.$fill_name" : "",
$fil_qname ? ".QName.$fil_qname" : "",
$fil_sname ? ".SName.$fil_sname" : "",
$keep_same_order?".kso":"",
$zoom_on_hits?".zoomedHits":"",
".png" if $blast_table;

$outfile||=join "",$reg_table,
$order_by_homology?".Homology_ordered":"",
$zoom_on_hits?".zoomedHits":"",
".png" if $reg_table;


$outfile||=join "",$dagchainer_table,
$order_by_homology ? ".Homology_ordered" : "",
$fil_alen ? ".alen$fil_alen" : "",
$fil_pid ? ".pid$fil_pid" : "",
$fil_eval ? ".eval$fil_eval" : "",
$fill_name ? ".Name.$fill_name" : "",
$fil_qname ? ".QName.$fil_qname" : "",
$fil_sname ? ".SName.$fil_sname" : "",
$keep_same_order?".kso":"",
$zoom_on_hits?".zoomedHits":"",
".png" if $dagchainer_table;









if ($reg_table) {

  $seq1||=1;
  $loc1||=2;
  $seq2||=3;
  $loc2||=4;

  $seq1-=$seq1;
  $loc1-=$loc1;
  $seq2-=$seq2;
  $loc2-=$loc2;
}



my%misc_count;
my%cont_order_info;
my%lg_cont_num;
my%linkage_map_info;
my$lg1;


### read linkage file and collect required info in hash %all_info
### %all_info needs contig_name, location on contig, linkage_Group, centimorghan.
read_mstmap($infile,\%all_info) if $infile;
read_table($lgtable,\%all_info) if $lgtable;
read_blast_table($blast_table,\%all_info,$resolution) if $blast_table;
read_dagchainer_table($dagchainer_table,\%all_info,$resolution) if $dagchainer_table;
read_reg_table($reg_table,\%all_info,$seq1,$loc1,$seq2,$loc2) if $reg_table;
###############################################################################################
###   Process linkage map data and plot graph.
###############################################################################################
print "\nProcessing collected information.....\n";

## turn off scale if $zoom_on_hits
if ($zoom_on_hits) {$scale_x=undef;$scale_y=undef}





### switch values if axis switch is requested.
if ($switch_axis) {
  print "\nSwitching axis for all the values.........";
  my%temp_all_info;

  foreach(keys %all_info){
    $temp_all_info{$_}{contig}=$all_info{$_}{lg};
    $temp_all_info{$_}{loc}=$all_info{$_}{cm};
    $temp_all_info{$_}{lg}=$all_info{$_}{contig};
    $temp_all_info{$_}{cm}=$all_info{$_}{loc};
    $temp_all_info{$_}{color}=$all_info{$_}{color} if $all_info{$_}{color};
  }
  undef %all_info;
  %all_info=%temp_all_info;
  undef %temp_all_info;

  print "done\n"
}


foreach my $mar_name(keys %all_info){
  my $contig_name = $all_info{$mar_name}{contig};
  my $location = $all_info{$mar_name}{loc};
  my $lg = $all_info{$mar_name}{lg};
  my $cm = $all_info{$mar_name}{cm};



  ### collect max , min info for linkage groups and LGs.
  $lg_info{$lg}{maxsize}||=$cm;  ## initialize
  $lg_info{$lg}{maxsize}=$cm if $lg_info{$lg}{maxsize}<$cm;  ## assign new if cm is larger

  $lg_info{$lg}{minsize}||=$cm;  ## initialize
  $lg_info{$lg}{minsize}=$cm if $lg_info{$lg}{minsize}>$cm;  ## assign new if cm is smaller


  $contig_info{$contig_name}{maxsize}||=$location;
  $contig_info{$contig_name}{maxsize}=$location if $contig_info{$contig_name}{maxsize} < $location;

  $contig_info{$contig_name}{minsize}||=$location;
  $contig_info{$contig_name}{minsize}=$location if $contig_info{$contig_name}{minsize} > $location;

  ### collect number of markers per contig per LG
  $lg_cont_num{$lg}{$contig_name}||=1;$lg_cont_num{$lg}{$contig_name}++;


  ### calculate average cm for each contig so that it can be ordered in dotplot in increasing or decreasing order.
  $misc_count{$contig_name}{$lg}||=0;$misc_count{$contig_name}{$lg}++;
  $cont_order_info{$lg}{$contig_name}{avgcm}||=$cm;
  $cont_order_info{$lg}{$contig_name}{avgcm} = ($cm + $misc_count{$contig_name}{$lg}*$cont_order_info{$lg}{$contig_name}{avgcm})/(1+$misc_count{$contig_name}{$lg});


  print "\n$lg\t$cm\t\t$contig_name\t$location" if $verbose;

}

### if the axis has been switched, reassign refseqlengths for x and y axis too.
if ($switch_axis) {
  print "\nSwitching axis information from reference information files the values.........";
  if ($refseqlenx && !$refseqleny){
    $refseqleny=$refseqlenx;
    undef $refseqlenx;
  }elsif (!$refseqlenx && $refseqleny){
      $refseqlenx=$refseqleny;
    undef $refseqleny;

  }elsif($refseqlenx && $refseqleny){
    my$temp_referencex=$refseqlenx;
    $refseqlenx=$refseqleny;
    $refseqleny=$temp_referencex;
    undef $temp_referencex;
  }elsif ($referencex && !$referencey){
    $referencey=$referencex;
    undef $referencex;
  }elsif (!$referencex && $referencey){
      $referencex=$referencey;
    undef $referencey;

  }elsif($referencex && $referencey){
    my$temp_referencex=$referencex;
    $referencex=$referencey;
    $referencey=$temp_referencex;
    undef $temp_referencex;
  }

}


##### get genome sizes from a sequence file or seqlen file to replace {maxsize} values in previously computed values.
if (!$zoom_on_hits) {  ### real length is not required if user wants to zoom on hit regions only
  if ($referencex) {print "\nGetting Sequence size information from Reference Sequence(referencex):$referencex\n";seq_length_fasta($referencex,\%contig_info);}
  if ($referencey) {print "\nGetting Sequence size information from Reference Sequence(referencey):$referencey\n";seq_length_fasta($referencey,\%lg_info);}
  if($refseqlenx){print "\nGetting Sequence size information from refseqlenx:$refseqlenx\n";seq_length_table($refseqlenx,\%contig_info);}
  if($refseqleny){print "\nGetting Sequence size information from refseqleny:$refseqleny\n";seq_length_table($refseqleny,\%lg_info);}
}


##### Sort the names to decide the order on axis. Default is by size.
my@contig_list=(sort {$contig_info{$b}{maxsize} <=> $contig_info{$a}{maxsize}} keys %contig_info);
my@lg_list=(sort {$lg_info{$b}{maxsize} <=> $lg_info{$a}{maxsize}} keys %lg_info);
print "\nSamples of values to be drawn on \nxaxis:\t",join ":",@lg_list[0..5],"\nYaxis:\t", join ":",@contig_list[0..5] if $verbose;
#foreach (keys %contig_info){print "\n ****This contig has info:$_:\tlength:$contig_info{$_}{maxsize}:";}
#### collect size sorted LGs and Contigs names
if($order_by_num){
  @contig_list=(sort {only_num($a) <=> only_num($b)} keys %contig_info);
  @lg_list=(sort {only_num($a) <=> only_num($b)} keys %lg_info);
}elsif ($order_by_alpha) {
  @contig_list=(sort {$a cmp $b} keys %contig_info);
  @lg_list=(sort {$a cmp $b} keys %lg_info);
}

if ($order_by_homology){
######################################################################################################
### sort contigs according to LG's relationship. closest first.
  my%temp_cinfo=%contig_info;  ## create a copy of contig info
  my@temp_clist;

  ## associate contigs to LG with max matches(dots)
  my%contig_for_lg;
  foreach my$lgt(@lg_list){
   foreach my$temp_cont(keys %{$lg_cont_num{$lgt}}){
     $contig_for_lg{$temp_cont}||=$lgt;
     $contig_for_lg{$temp_cont}=$lgt if $lg_cont_num{$lgt}{$temp_cont} > $lg_cont_num{$contig_for_lg{$temp_cont}}{$temp_cont};
   }
  }

#### sort contigs based on association and cm position to LGs



  print "\nSorting contigs by order\n";
  open(LIST,">$outfile.Names.list") if $print_list;
  foreach my$lgt(@lg_list){
    my@print_list;
    #my @tconarray=sort{$lg_cont_num{$lgt}{$b}<=>$lg_cont_num{$lgt}{$a}} keys %{$lg_cont_num{$lgt}};  ### sort lgs based on number of dots in each
      my @tconarray=sort{$cont_order_info{$lgt}{$a}{avgcm}<=>$cont_order_info{$lgt}{$b}{avgcm}} keys %{$cont_order_info{$lgt}};  ### sort lgs based on number of position of contigs in linkage group
    foreach (@tconarray){
      #push(@temp_clist,$_) if $temp_cinfo{$_}{maxsize};
      push(@temp_clist,$_) if ($contig_for_lg{$_} eq $lgt  && $lg_cont_num{$lgt}{$_} >= $min_markers);  ### collect contig name is if it is closest to this LG.
      push(@print_list,$_) if ($contig_for_lg{$_} eq $lgt  && $lg_cont_num{$lgt}{$_} >= $min_markers);

      delete $temp_cinfo{$_} if ($contig_for_lg{$_} eq $lgt && $lg_cont_num{$lgt}{$_} >= $min_markers && $temp_cinfo{$_});  ## delete contig entry from temp_cinfo hash
    }

    print LIST join ";",$lgt,@print_list,"\n"  if $print_list;
  }
  my@resid_clist=(sort {$temp_cinfo{$b}{maxsize} <=> $temp_cinfo{$a}{maxsize}} keys %temp_cinfo);
  @contig_list=(@temp_clist,@resid_clist);

}

if($keep_same_order){print "\nAssigning Names on Xaxis to same as Y axis to keep in raangement same order (-kso flag is on)." if $verbose; @lg_list=@contig_list;}

if ($order_y) {
  open(UDY,"$order_y");


  while (<UDY>) {
      s/^\s*|\s*$//;
    next if m/^\s*$/;
    my@temp_udl=split /,/;
    @lg_list=@temp_udl;
    }
  }

if ($order_x) {
  open(UDX,"$order_x");


  while (<UDX>) {
      s/^\s*|\s*$//;
    next if m/^\s*$/;
    my@temp_udl=split /,/;
    @contig_list=@temp_udl;
    }
}

######################################################################################################
print "\n2..Samples of values to be drawn on \nxaxis:\t",join ":",@lg_list[0..5],"\nYaxis:\t", join ":",@contig_list[0..5] if $verbose;
my$sum_contig=0.0001;
my$sum_lg=0.0001;

### assign number of lg and contigs to plot if not provided.
$lg_top||=scalar@lg_list;
$cont_top||=scalar@contig_list;
##readjust the top values in case there are not enough contigs or lgs
$lg_top=$lg_top < scalar@lg_list?$lg_top:scalar@lg_list;
$cont_top=$cont_top < scalar@contig_list?$cont_top:scalar@contig_list;
##### assign starting coordinates to each lg and contig
print "\nFinding contig boundries\n";
foreach(@contig_list[0..$cont_top-1]){
  my $zoom_mod=$zoom_on_hits?$contig_info{$_}{minsize}:0;
  $cumulative_contig{$_}{start}=$sum_contig; ## start of contig is after end of previous all contigs on axis.
  $sum_contig+=($contig_info{$_}{maxsize} - $zoom_mod);
  #$sum_contig-=$zoom_mod if $zoom_on_hits;  ## zoom on start by substracting start values from total
  $cumulative_contig{$_}{end}=$sum_contig;
  #print "\nContig:$_\tSize:$contig_info{$_}{maxsize}";
  }
foreach(@lg_list[0..$lg_top-1]){
  my $zoom_mod=$zoom_on_hits?$lg_info{$_}{minsize}:0;
  $cumulative_lg{$_}{start}=$sum_lg;
  $sum_lg+=($lg_info{$_}{maxsize} - $zoom_mod );
  #$sum_lg-=$zoom_mod if $zoom_on_hits; ## zoom on start by substracting start values from total ;
  $cumulative_lg{$_}{end}=$sum_lg;
  }

#### calculate cm/pixel and nt/pixel for drawing
my$nt_per_pix=$sum_contig/$img_width;
my$cm_per_pix=$sum_lg/$img_height;


if($same_scale){
  $nt_per_pix=$cm_per_pix > $nt_per_pix?$cm_per_pix:$nt_per_pix;
  $cm_per_pix=$nt_per_pix;
}

#### mpdify per pixel stat if user wants to draw off scale for equal space per query.
my%dist_per_pix;
my@temp_cont_list=@contig_list[0..$cont_top-1];
my$len_per_block_contig=$img_width/scalar@temp_cont_list;
my$effective_size_contig=$sum_contig;
foreach(@contig_list[0..$cont_top-1]){
  my $zoom_mod=$zoom_on_hits?$contig_info{$_}{minsize}:0;
  $effective_size_contig=$contig_info{$_}{maxsize} - $zoom_mod;# if $zoom_on_hits;
  $dist_per_pix{$_}=$nt_per_pix;
  $dist_per_pix{$_}=$effective_size_contig/$len_per_block_contig if $equal_block;

}

my@temp_lg_list=@lg_list[0..$lg_top-1];
my$len_per_block_lg=$img_height/scalar@temp_lg_list;
my$effective_size_lg=$sum_lg;

foreach(@lg_list[0..$lg_top-1]){
  my $zoom_mod=$zoom_on_hits?$lg_info{$_}{minsize}:0;
  $effective_size_lg=$lg_info{$_}{maxsize} - $zoom_mod;# if $zoom_on_hits;
  $dist_per_pix{$_}=$cm_per_pix;
  $dist_per_pix{$_}=$effective_size_lg/$len_per_block_lg if $equal_block;
}


#print "\nNt per pix:$nt_per_pix\nCm per pixel:$cm_per_pix\nSum contig:$sum_contig\nSum LG:$sum_lg\n";
print "\nDrawing on canvas\n";
# Create and print the output.

if(defined $rot_xy_lab){$rot_x_lab=$rot_xy_lab;$rot_y_lab=$rot_xy_lab;}


    my $img = GD::Simple->new($img_width+2*$margin,$img_height+2*$margin);
    my$white = $img->colorAllocate(255,255,255);
	my$black = $img->colorAllocate(0,0,0);
    my$red = $img->colorAllocate(255,0,0);
    my$gray = $img->colorAllocate(187,187,187);

    my$grid_color=$gray;
    $grid_color=$white if $no_lines;


    $img->interlaced('true');

    ## draw outline of the graph area
    $img->rectangle(0+$margin,0+$margin,$img_width+ $margin,$img_height+ $margin,$black);

    ## draw contig and lg boundry lines to graph and add labels
    print "\nDrawing line seperating contigs\n";
    $img->penSize(1,1);

    foreach(@contig_list[0..$cont_top-1]){
      my $zoom_mod=$zoom_on_hits?$contig_info{$_}{minsize}:0;
      $img->line($cumulative_contig{$_}{end}/$dist_per_pix{$_} + $margin, 0 + $margin,
                $cumulative_contig{$_}{end}/$dist_per_pix{$_} + $margin, $img_height + $margin ,
                $grid_color);

      $img->moveTo(($cumulative_contig{$_}{start} + (($contig_info{$_}{maxsize} - $zoom_mod)/2) + $font_size/2)/$dist_per_pix{$_} + $margin,
                   $margin-5);
      $img->font('gdMediumBoldFont');
      $img->fontsize($font_size);
      $img->angle($rot_x_lab);
      $img->string("$_");
      
    }
    
    ## add scale bar on x axis
        
    if ($scale_x||$scale_x_dist) {
      my$sbar_len=10;  ## length of scale bars
      my$scale_dist=$scale_x_dist ? $scale_x_dist : int(10**(length(int($sum_contig/10))));
      my$scale_mark=$scale_dist == 1000 ? "Kb" : $scale_dist == 1000000 ? "Mb": "(x$scale_dist bp)";
      
      
      print "\nScale used on X axis for $sum_contig bp is $scale_dist\n";
      for(my$x=0;$x<$sum_contig;$x+=$scale_dist){
        my$xloc1=$x/$nt_per_pix+$margin;
        my$yloc1=$img_height+$margin;
        my$xloc2=$x/$nt_per_pix+$margin;
        my$yloc2=$img_height + $margin + $sbar_len;
        
        $img->line($xloc1,$yloc1,$xloc2,$yloc2);
        $img->moveTo($xloc2,$yloc2+$font_size/2+$sbar_len);
        $img->font('gdMediumBoldFont');
        $img->fontsize($font_size);
        $img->angle(90);
        my$scale_num = $x/$scale_dist;
        $img->string(join(" ",$scale_num,$scale_mark));
      }
    }
    
    
    foreach(@lg_list[0..$lg_top-1]){
      my $zoom_mod=$zoom_on_hits?$lg_info{$_}{minsize}:0;
      $img->line(0+ $margin,  ## x1
                 $cumulative_lg{$_}{end}/$dist_per_pix{$_} + $margin,  ## y1
                 $img_width + $margin, ## x2
                 $cumulative_lg{$_}{end}/$dist_per_pix{$_} + $margin, ## y2
                 $grid_color);


      #$img->moveTo( $margin - 5 - $font_size * 3,  ## x
      $img->moveTo( int($margin/10),  ## x
                   ($cumulative_lg{$_}{start} + ($lg_info{$_}{maxsize}-$zoom_mod)/2 )/$dist_per_pix{$_}  + $margin + $font_size/2  ## y
      );
      $img->font('gdMediumBoldFont');
      $img->fontsize($font_size);
      $img->angle($rot_y_lab);
      $img->string("$_");
    }
    
    ### add scale y
              ## add scale bar on x axis
    if ($scale_y||$scale_y_dist) {
      my$sbar_len=10;  ## length of scale bars
      my$scale_dist=$scale_y_dist ? $scale_y_dist : int(10**(length(int($sum_lg/10))));
      my$scale_mark=$scale_dist == 1000 ? "Kb" : $scale_dist == 1000000 ? "Mb": "(x$scale_dist bp)";
      print "\nScale used on Y axis : $sum_lg bp is $scale_dist\n";
      for(my$x=0;$x<$sum_lg;$x+=$scale_dist){
        my$xloc1=$margin + $img_width;
        my$yloc1=$margin + $x/$cm_per_pix;
        my$xloc2=$margin +$img_width + $sbar_len; 
        my$yloc2=$margin + $x/$cm_per_pix;
        
        $img->line($xloc1,$yloc1,$xloc2,$yloc2);
        $img->moveTo($xloc2,$yloc2+$font_size/2+$sbar_len);
        $img->font('gdMediumBoldFont');
        $img->fontsize($font_size);
        $img->angle(0);
        my$scale_num = $x/$scale_dist;
        
        $img->string(join(" ",$scale_num,$scale_mark));
      }
    }
    
    

### start placing dots for each marker.
print "\nAdding dots to the canvas\n";
    $img->penSize(1,1);
foreach my$marker(keys %all_info){
  print "\nSkipping Marker :$marker on contig:$all_info{$marker}{contig}" if  !$cumulative_contig{$all_info{$marker}{contig}}{start}  && $verbose;
    next if !$cumulative_contig{$all_info{$marker}{contig}}{start};
    next if !$cumulative_lg{$all_info{$marker}{lg}}{start};

    ## zoom on hits modify the start values by substracting min values.
    my$mod_zoom_contig = $zoom_on_hits? $contig_info{$all_info{$marker}{contig}}{minsize} :0;
    my$mod_zoom_lg     = $zoom_on_hits? $lg_info{$all_info{$marker}{lg}}{minsize}         :0;

    my$xvalue=($all_info{$marker}{loc}+ $cumulative_contig{$all_info{$marker}{contig}}{start} - $mod_zoom_contig)/$dist_per_pix{$all_info{$marker}{contig}}  + $margin;
    my$yvalue=($all_info{$marker}{cm} + $cumulative_lg{$all_info{$marker}{lg}}{start} - $mod_zoom_lg)/$dist_per_pix{$all_info{$marker}{lg}}  + $margin;
    my$color||='red';
    $color=$all_info{$marker}{color} if $all_info{$marker}{color};
    #next if ($xvalue > $img_width || $yvalue > $img_height);

    $img->moveTo($xvalue,$yvalue);
    $img->bgcolor($color);
    $img->fgcolor($color);
    $img->ellipse($dotsize,$dotsize) if !$as_marker;
    $img->string("$as_marker") if $as_marker;
    #print "\nLoc:$all_info{$marker}{loc}\tContStart:$cumulative_contig{$all_info{$marker}{contig}}{start}\nprinting dot at $xvalue,$yvalue";
}
    # Print this thing out
    print "\nPrinting image to file.\n";
    $img->moveTo($margin,$img_height+$margin+20);
    $img->angle(0);
    $img->font('gdMediumBoldFont');
    ##$img->string(gdMediumBoldFont,$margin,$img_height+$margin+20,"$outfile",$red);
    ### print file name at the bottom
    $img->bgcolor('red');
    $img->fgcolor('red');
    $img->string("$outfile") if !defined $no_file_name;
    add_blast_legend(\$img, $margin, $img_height+ $margin + $margin/2) if $blast_table;
    add_dagchainer_legend(\$img, $margin, $img_height+ $margin + $margin/2) if $dagchainer_table;
    print_png($img, $outfile) if !$svg;
    print_svg($img, "$outfile.svg")  if $svg;

    print "Plot Complete!\n";

#########################################################
sub print_png{
    my ($img1, $out) = @_;
    $out||="dotPlot.png";
	open(OUT, ">$out") || die "Cannot write $out: $!\n";
	print OUT $img1->png;
    print "\nSaved PNG image as $out\n";
	close OUT;
}
sub print_svg{
    my ($img1, $out) = @_;
    $out||="dotPlot.svg";
	open(OUT, ">$out") || die "Cannot write $out: $!\n";
	print OUT $img1->svg;
    print "\nSaved SVG image as $out\n";
	close OUT;
}

sub seq_length_fasta{
  my$file=shift;
  my$hshRef=shift;
  $/="\n>";
  my$totcount=my$selcount=0;
  my$totseqlen=my$selseqlen=0;
  print "\nCannot find the file $file\n" if ! -f $file;
  open(FASTA,"$file") or return undef;
  while(<FASTA>){
    my($header,@sequence)=split(/\n/,$_);
    $header=~s/^>|^\s+|\s*$//g;
    my$sequence=join("",@sequence);
    $sequence=~s/\s+|>//g;
    my$len=0;
    $len=length($sequence);
    if (exists $$hshRef{$header}{maxsize}){
      $$hshRef{$header}{maxsize}=$len;
      $selcount++;
      $selseqlen+=$$hshRef{$header}{maxsize};
    }
    $totcount++;
    $totseqlen+=$len;
  }

  print "\n *******Mapping Summary*******\n$totcount Total Sequences read with total length $totseqlen";
  print "\n $selcount Sequences mapped with total length $selseqlen";
  print "\n",$selcount*100/$totcount,"% of total number was mapped";
  print "\n",$selseqlen*100/$totseqlen,"% of total length was mapped\n";
  close FASTA;
}


sub seq_length_table{
  my$file=shift;
  my$hshRef=shift;
  my$totcount=my$selcount=0;
  my$totseqlen=my$selseqlen=0;
  print "\nCannot find the file $file\n" if ! -f $file;
  open(TABLE,"$file") or return undef;
  while(<TABLE>){
    s/^\s+//g;
    next if m/^\s*$/;
    my($header,$len,$rest)=split(/\s+/,$_);
    $header=~s/>|^\s+|\s*$//g;
    $len=~s/\D+//g;
    ### only add info in the hash if it already exists. Otherwise it will add spurious names and mess up the plot
    if (exists $$hshRef{$header}){
      $$hshRef{$header}{maxsize}=$len;
      $selcount++;
      $selseqlen+=$$hshRef{$header}{maxsize};
      print "\n Updated the length for $header to $len (len):$$hshRef{$header}{maxsize}(from hash)" if $verbose;
    }else{
      print "\n Did not update the length for $header as it does not exists in hash" if $verbose;;
    }
    $totcount++;
    $totseqlen+=$len;
    print "\nUsing sequence length: $header\t$len" if $verbose;
  }
  print "\n";
  print "\n *******Mapping Summary*******\n$totcount Total Sequences read with total length $totseqlen";
  print "\n $selcount Sequences mapped with total length $selseqlen";
  print "\n",$selcount*100/$totcount,"% of total number was mapped";
  print "\n",$selseqlen*100/$totseqlen,"% of total length was mapped\n";
  close TABLE;
}


sub read_mstmap{
  my $file=shift;
  my $ref_hash=shift;
  open( LG,"$file") or die "\n Unable to open linkage file\n$usage\n";
  while (<LG>) {
    next if m/^;/;next if m/^\s*$/;s/^\s+//g;s/\s*$//g;
    $lg1 = $1 if m/^group\s+([\S]+)/; next if m/^group\s+([\S]+)/;

    my($contig1,$cm1)=split /\s+/;
    my@info=split /_/,$contig1;
    s/\s+//g foreach @info;
    next if $#info < 3;

    my $mark_name=$info[0];
    my $contig_name1=join "_",@info[1..2];
    my $location1=$info[3];

    $$ref_hash{$mark_name}{contig}=$contig_name1;
    $$ref_hash{$mark_name}{loc}= $location1;
    $$ref_hash{$mark_name}{lg}=$lg1;
    $$ref_hash{$mark_name}{cm}=$cm1;
  }
}

sub read_table{
  my $file=shift;
  my $ref_hash=shift;
  open( TABLE,"$file") or die "\n Unable to open linkage file\n$usage\n";
  while (<TABLE>) {
    next if m/^;/;next if m/^\s*$/;s/^\s+//g;s/\s*$//g;
    my($mark_name,$contig_name,$location,$lg,$cm,$ext)=split /\s+/;
    $$ref_hash{$mark_name}{contig}=$contig_name;
    $$ref_hash{$mark_name}{loc}= $location;
    $$ref_hash{$mark_name}{lg}=$lg;
    $$ref_hash{$mark_name}{cm}=$cm;
  }
}


sub read_blast_table{
  my $file=shift;
  my $ref_hash=shift;
  my $resolution=shift;
  $resolution||=10000;
  my$count=0;
  open( BLAST,"$file") or die "\n Unable to open linkage file\n$usage\n";
  print "\n Processing Blast file......";
  while (<BLAST>) {
    next if m/^#/;next if m/^\s*$/;s/^\s+//g;s/\s*$//g;
    my@bl_line=split /\s+/;

    next if ( $fil_pid && $bl_line[2] < $fil_pid );
    next if ( $fil_alen && $bl_line[3] < $fil_alen );
    next if ( $fil_eval && $bl_line[10] > $fil_eval );
    next if ( $fill_name  && ($bl_line[0] !~ m/$fill_name/i && $bl_line[1] !~ m/$fill_name/i));
    next if ( $fill_name && $exact && ($bl_line[0] ne $fill_name &&  $bl_line[1] ne $fill_name));
    next if ( $fil_qname && $bl_line[0] !~ m/$fil_qname/i);
    next if ( $fil_sname && $bl_line[1] !~ m/$fil_sname/i);

    next if ( $fil_qname && $exact && $bl_line[0] ne $fil_qname);
    next if ( $fil_sname && $exact && $bl_line[1] ne $fil_sname);



    $count++;
    my $mark_name= "blast_$count" ;
    my $contig_name=$bl_line[0]  ;
    my $location= $bl_line[6] ;
    my $lg=  $bl_line[1];
    my $cm= $bl_line[8] ;
    my $color=color_by_pid($bl_line[2]);

    if (abs($bl_line[6] - $bl_line[7])/$resolution >= 2 ) {
      line_to_dots($bl_line[6],$bl_line[7],$bl_line[8],$bl_line[9],$resolution,$ref_hash,$mark_name,$contig_name,$lg,$color);
    }else{
      $$ref_hash{$mark_name}{contig}=$contig_name;
      $$ref_hash{$mark_name}{loc}= $location;
      $$ref_hash{$mark_name}{lg}=$lg;
      $$ref_hash{$mark_name}{cm}=$cm;
      $$ref_hash{$mark_name}{color}=$color;
    }
  }

  print "\n Finished processing Blast file\n";
}


sub read_dagchainer_table{
  my $file=shift;
  my $ref_hash=shift;
  my $resolution=shift;
  $resolution||=10000;
  my$count=0;
  open( BLAST,"$file") or die "\n Unable to open linkage file\n$usage\n";
  print "\n Processing Blast file......";
  while (<BLAST>) {
    next if m/^#/;next if m/^\s*$/;s/^\s+//g;s/\s*$//g;
    my@bl_line=split /\s+/;

    #next if ( $fil_pid && $bl_line[2] < $fil_pid );
    #next if ( $fil_alen && $bl_line[3] < $fil_alen );
    next if ( $fil_eval && $bl_line[8] > $fil_eval );
    #next if ( $fill_name  && ($bl_line[0] !~ m/$fill_name/i && $bl_line[1] !~ m/$fill_name/i));
    #next if ( $fill_name && $exact && ($bl_line[0] ne $fill_name &&  $bl_line[1] ne $fill_name));
    next if ( $fil_qname && $bl_line[0] !~ m/$fil_qname/i);
    next if ( $fil_sname && $bl_line[4] !~ m/$fil_sname/i);

    next if ( $fil_qname && $exact && $bl_line[0] ne $fil_qname);
    next if ( $fil_sname && $exact && $bl_line[4] ne $fil_sname);



    $count++;
    my $mark_name= "blast_$count" ;
    my $contig_name=$bl_line[0]  ;
    my $location= $bl_line[2] ;
    my $lg=  $bl_line[4];
    my $cm= $bl_line[6] ;
    my $color=color_by_pval($bl_line[8]);

    if (abs($bl_line[2] - $bl_line[3])/$resolution >= 2 ) {
      line_to_dots($bl_line[2],$bl_line[3],$bl_line[6],$bl_line[7],$resolution,$ref_hash,$mark_name,$contig_name,$lg,$color);
    }else{
      $$ref_hash{$mark_name}{contig}=$contig_name;
      $$ref_hash{$mark_name}{loc}= $location;
      $$ref_hash{$mark_name}{lg}=$lg;
      $$ref_hash{$mark_name}{cm}=$cm;
      $$ref_hash{$mark_name}{color}=$color;
    }
  }

  print "\n Finished processing dagchainer file\n";



}






sub read_reg_table{
  my $file=shift;
  my $ref_hash=shift;
  my $tseq1=shift;
  my $tloc1=shift;
  my $tseq2=shift;
  my $tloc2=shift;

  my$count=0;
  open( BLAST,"$file") or die "\n Unable to open linkage file\n$usage\n";
  print "\n Processing Blast file......";
  while (<BLAST>) {
    next if m/^#|^;/;next if m/^\s*$/;s/^\s+//g;s/\s*$//g;
    my@bl_line=split /\s+/;
    $count++;
    my $mark_name= "blast_$count" ;
    my $contig_name=$bl_line[$tseq1]  ;
    my $location= $bl_line[$tloc1] ;
    my $lg=  $bl_line[$tseq2];
    my $cm= $bl_line[$tloc2] ;

    ## populate the hash %all_info.
    $$ref_hash{$mark_name}{contig}=$contig_name;
    $$ref_hash{$mark_name}{loc}= $location;
    $$ref_hash{$mark_name}{lg}=$lg;
    $$ref_hash{$mark_name}{cm}=$cm;

  }

  print "\n Finished processing Blast file\n";
}


sub line_to_dots{
  my $x1=shift;
  my $x2=shift;
  my $y1=shift;
  my $y2=shift;
  my $resolution=shift;
  my $ref_hash=shift;
  my $markname=shift;
  my $contigname=shift;
  my $lgname=shift;
  my $color=shift;

  my$numfragsx=($x2-$x1)/$resolution;
  my$orient_x=$numfragsx > 0?1:-1;
  my$numfragsy=($y2-$y1)/$resolution;;
  my$orient_y=$numfragsy > 0?1:-1;

  for(my$i=1;$i <= abs($numfragsx);$i++){
    my$tmarkname="$markname.$i";

    $$ref_hash{join "_",$markname,$i}{contig}=$contigname;
    $$ref_hash{join "_",$markname,$i}{loc}=$x1 + ($i * $resolution * $orient_x);
    $$ref_hash{join "_",$markname,$i}{lg}=$lgname;
    $$ref_hash{join "_",$markname,$i}{cm}=$y1 +  ($i * (abs($y2-$y1)/abs($numfragsx)) * $orient_y);
    $$ref_hash{join "_",$markname,$i}{color}=$color;


    #$$ref_hash{join "_",$markname,$i}{contig}=$contigname;
    #$$ref_hash{join "_",$markname,$i}{loc}=$orient_x>0?$x1:$x2 + ($i * $resolution * $orient_x);
    #$$ref_hash{join "_",$markname,$i}{lg}=$lgname;
    #$$ref_hash{join "_",$markname,$i}{cm}=$orient_y>0?$y1:$y2 +  ($i * (abs($y2-$y1)/abs($numfragsx)) * $orient_y);
    #$$ref_hash{join "_",$markname,$i}{color}=$color;


  }

  print "\ndividing blast hit between $contigname ($x1-$x2) and $lgname ($y1-$y2) into $numfragsx (x)/$numfragsy (y) dots" if $verbose;


}

sub color_by_pid{
  my $percent_id=shift;
  #print "\n Finding color for $percent_id";
  if ($percent_id >= 95 ) {
    return "red";
  }  elsif ($percent_id >= 90 && $percent_id < 95 ) {
    return "magenta";
  }  elsif ($percent_id >= 85 && $percent_id < 90 ) {
    return "orange";
  }  elsif ($percent_id >= 80 && $percent_id < 85 ) {
    return "lightgreen";
  }  elsif ($percent_id >= 75 && $percent_id < 80 ) {
    return "green";
  }  elsif ($percent_id >= 70 && $percent_id < 75 ) {
    return "skyblue";
  }  elsif ($percent_id > 65 && $percent_id < 70 ) {
    return "steelblue";
  }  elsif ($percent_id > 60 && $percent_id < 65 ) {
    return "blue";
  }  elsif ($percent_id > 55 && $percent_id < 60 ) {
    return "yellow";
  }  elsif ($percent_id > 50 && $percent_id < 55 ) {
    return "khaki";
  }  elsif ($percent_id > 45 && $percent_id < 50 ) {
    return "black";
  }
  else{ print "\nUnable to assign color for value:$percent_id\n";return undef}

}


sub color_by_pval{

  my $percent_id=shift;

#print "\n Finding color for $percent_id";
  if ($percent_id <= 1e-150 ) {
    return "red";
  }  elsif ($percent_id >= 1e-150 && $percent_id < 1e-100 ) {
    return "magenta";
  }  elsif ($percent_id >= 1e-100 && $percent_id < 1e-50 ) {
    return "orange";
  }  elsif ($percent_id >= 1e-50 && $percent_id < 1e-20 ) {
    return "lightgreen";
  }  elsif ($percent_id >= 1e-20 && $percent_id < 1e-10 ) {
    return "green";
  }  elsif ($percent_id >= 1e-10 && $percent_id < 1e-5 ) {
    return "skyblue";
  }  elsif ($percent_id >= 1e-5 && $percent_id < 0.0001 ) {
    return "steelblue";
  }  elsif ($percent_id >= 0.0001 && $percent_id < 0.01 ) {
    return "blue";
  }  elsif ($percent_id >= 0.01 && $percent_id < 0.1 ) {
    return "yellow";
  }  elsif ($percent_id >= 0.1 && $percent_id < 1 ) {
    return "khaki";
  }  elsif ($percent_id >= 1 && $percent_id < 10 ) {
    return "black";
  }
  else{ print "\nUnable to assign color for value:$percent_id\n";return undef}


}


sub add_blast_legend{
  my $img_ref=shift;
  my $loc_x=shift;
  my $loc_y=shift;

  my@text=('100-95%','95-90%','90-85%','85-80%','80-75%','75-70%','70-65%','65-60%','60-55%','55-50%','50-45%');
  my@color=qw(red magenta orange lightgreen green skyblue steelblue blue yellow khaki black);
  my$count_num=0;
  my $cur_locx=$loc_x;
  for(my$i=0;$i<=$#text;$i++){
    my$count_num_y=0;
    #print "\nAdding Text:$text[$i] with color:$color[$i]";
    $$img_ref->moveTo($cur_locx,$loc_y);
    $$img_ref->angle(0);
    $$img_ref->font('gdMediumBoldFont');
    ##$img->string(gdMediumBoldFont,$margin,$img_height+$margin+20,"$outfile",$red);
    ### print file name at the bottom
    $$img_ref->bgcolor($color[$i]);
    $$img_ref->fgcolor($color[$i]);
    $$img_ref->string("$text[$i]");
    $cur_locx+=70;
    if ($cur_locx+70 > $img_width+$margin) {
      $cur_locx=$loc_x;
      $loc_y+=25
    }
  }
}


sub add_dagchainer_legend{
  my $img_ref=shift;
  my $loc_x=shift;
  my $loc_y=shift;

  my@text=('<1e-150','<1e-100','<1e-50','<1e-20','<1e-10','<0.0001','<0.001','0.01','<0.1','<1','<10');
  my@color=qw(red magenta orange lightgreen green skyblue steelblue blue yellow khaki black);
  my$count_num=0;
  my $cur_locx=$loc_x;
  for(my$i=0;$i<=$#text;$i++){
    my$count_num_y=0;
    $$img_ref->moveTo($cur_locx,$loc_y);
    $$img_ref->angle(0);
    $$img_ref->font('gdMediumBoldFont');

    $$img_ref->bgcolor($color[$i]);
    $$img_ref->fgcolor($color[$i]);
    $$img_ref->string("$text[$i]");
    $cur_locx+=70;
    if ($cur_locx+70 > $img_width+$margin) {
      $cur_locx=$loc_x;
      $loc_y+=25
    }
  }

}







sub only_num{
  my$string=shift;
  $string=~s/[^\d\.]+//g;
  return $string;

}
