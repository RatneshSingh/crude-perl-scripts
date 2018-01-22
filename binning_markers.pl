#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
our($loc_file,$bin_size,$out_file,$help,$max_marker_percontig,$verbose);



my$usage="
perl $0 -l joinmap_file  -b bin_size  -o outfile
-b  bin_size to binmarkers[100000]
-m max marker per contig.This will deactivate size
    binning and will pick max marker accross the length.
";

GetOptions(
  "loc_file|l|i=s"=>\$loc_file,
  "bin_size|binsize|b=i"=>\$bin_size,
  "max_marker_percontig|m=i" =>\$max_marker_percontig,
  "out|o=s"=>\$out_file,
  "verbose|v"=>\$verbose,
  "help|h"=>\$help
) or die "\nError in command line arguments\n";

die "$usage" if  $help||!$loc_file;
$bin_size||=100000;

$out_file||="$loc_file.$bin_size.binned.loc";

my%markers;

##Read Joinmap
open(LOC,"$loc_file") or die "\nUnable to find loc file:$loc_file\n";
my($name,$popt,%bin,@array_loc,@array_ind);
my$start_ind=0;
while (<LOC>) {
  next if m/^\s*$/;
  s/^\s+//g;
  if (m/^name\s*=([\S]+)/){
    $name=$1;
    print "\nFound Name:$name on line $.";
  }
  if (m/^popt=\s*([\S]+)/){
    $popt=$1;
    print "\nFound Pop type:$popt on line $.";
  }
  next if m/^popt\s*=/i;
  next if m/^nloc\s*=/i;
  next if m/^nind\s*=/i;
  next if m/^name\s*=/i;
  #### 40_contig_0_1034818 <nnxnp>

  my$line=$_;
  $line=~s/\s+$//g;
  $line=~s/^\s+//g;
  $line=~s/\n+//g;
  if ($start_ind < 1) {
    $start_ind = 1 if m/^individual\s*names\s*:/;
    next if m/^individual\s*names\s*:/;
    #print "\n$.\t ****processing $line";
    my@array=split /\s+/,$line;
    my($markername,$contigname,$contignumber,$contiglocation) = split /_/,$array[0];

    my$mtype=$array[1];
    my$mtype1=$1 if $mtype=~m/<([\w]{2})x([\w]{2})>/i;
    my$mtype2=$2 if $mtype=~m/<([\w]{2})x([\w]{2})>/i;
    my@tarray=@array[2..$#array];
   # print "\n$line\nM1:$mtype1\tM2:$mtype2\tchisq for $markername: Calc:\t";
    my$chisq=get_chi($mtype1,$mtype2,\@tarray);

   # print "$chisq";# List:",$markers{$markername}{qual};

    $markers{$markername}{qual}=$chisq;

    my$contig=join "_",$contigname,$contignumber;
    $markers{$markername}{contig}=$contig;
    $markers{$markername}{loc}=$contiglocation;
    $markers{$markername}{line}=$line;
    my$bin=int($contiglocation/$bin_size)+1;
    $bin{$contig}{$bin}||=[];
    push(@{$bin{$contig}{$bin}},$markername);
  }

  if ($start_ind == 1) {
    push(@array_ind,$line);
  }






}

my%unused_markers=%markers;



### bin filtering.
my%frac_bin;
my@max_array_loc=();
foreach my$cont(keys %bin){
  my@all_markers;

  my$max_location=0;
  foreach my$bin(sort {$a <=> $b} keys %{$bin{$cont}}){
    my$best_mark=${$bin{$cont}{$bin}}[0];
    push(@all_markers,$best_mark) if (scalar@{$bin{$cont}{$bin}} == 1);
    $max_location=$markers{$best_mark}{loc};
    if (scalar@{$bin{$cont}{$bin}} > 1) {
      foreach my$mar(@{$bin{$cont}{$bin}}){
        $max_location=$markers{$mar}{loc} if $markers{$mar}{loc} > $max_location;
        $best_mark=$mar if (exists $markers{$mar}{qual} && $markers{$mar}{qual} < $markers{$best_mark}{qual});

        push(@all_markers,$mar);  ### collecting marker names for "max_mark_per_contig" option
      }
    }

    push(@array_loc,$markers{$best_mark}{line});
  }

  ## if there are more than $max_marker_percontig markers in each contig, reduce the markers to user assigned number.
  ## create new bin structure based on user assigned number and allocate markers in bins
  if ($max_marker_percontig) {

    my $frag_size=$max_location/($max_marker_percontig-1);
    my$num_mark=scalar@all_markers;
    if($num_mark > $max_marker_percontig ) {
      print "\n\n****Found $num_mark on contig:$cont. Binning to reduce the number of marker to $max_marker_percontig" if $verbose;
      foreach (@all_markers){
       push(@{$frac_bin{$cont}{int($markers{$_}{loc}/$frag_size)}},$_); ### collect and process later to pick required number of markers.
       print "\n\t\tplacing marker $_ \( at $markers{$_}{loc} on $cont \) in bin # ",int($markers{$_}{loc}/$frag_size) if $verbose;
     }
    }else{
      print "\n**** Existing ", scalar@all_markers," :Markers on $cont are less than $max_marker_percontig. Using all the markers" if $verbose;
      push(@max_array_loc, map {$markers{$_}{line} } @all_markers);
    }
  }
}





#### limit number of markers to $max_mark_per_contig if asked.

if ($max_marker_percontig) {
  foreach my$cont(keys %frac_bin){
     my$num_bins=scalar(keys %{$frac_bin{$cont}});

      if($num_bins > $max_marker_percontig )  { print "\n*****The number of bins:$num_bins seems more than requested $max_marker_percontig marker per contig\n" if $verbose; }

      foreach my$bin(sort {$a <=> $b} keys %{$frac_bin{$cont}}){        ### pick one best marker from each bin
          #my$best_mark=${$frac_bin{$cont}{$bin}}[0];
          #if (scalar@{$frac_bin{$cont}{$bin}} > 1) {
          my $best_mark=get_best_marker(\@{$frac_bin{$cont}{$bin}},\%markers,\%unused_markers) ;
          #}
          push(@max_array_loc,$markers{$best_mark}{line});                 ### collected selected best marker per bin
          print "\nSelected $best_mark for bin:$bin on contig:$cont" if $verbose;
      }

      if($num_bins < $max_marker_percontig )  {   ## select more markers if not enough markers
        print "\n*****The number of bins:$num_bins seems less than requested $max_marker_percontig marker for contig:$cont. Picking ",$max_marker_percontig-$num_bins," Random markers." if $verbose;
         for (1..($max_marker_percontig-$num_bins)){
           my @tbin = sort {$a <=> $b} keys %{$frac_bin{$cont}};
           my $ext_mark;
           my $try=0;
              while(!$ext_mark){
                $try++;
                my$rand_elem=$tbin[rand @tbin];
                $ext_mark = get_best_marker(\@{$frac_bin{$cont}{$rand_elem  }},\%markers,\%unused_markers);
                print "\nTry $try ........No Markers left to select from for $cont in bin $rand_elem" if ($verbose && scalar @{$frac_bin{$cont}{ $rand_elem }} ==0 );
                last if $try > $max_marker_percontig; ### kill the loop after few tries.


                } ## pick markers randomly for the rest.
              push(@max_array_loc,$markers{$ext_mark}{line});
              print "\n\tRandomly Selected $_ marker: $ext_mark for contig:$cont\n" if $verbose;
         }
      }
  }
}






if ($max_marker_percontig){

  write_joinmap($name,$popt,\@max_array_loc,\@array_ind,$out_file);

}else{

  write_joinmap($name,$popt,\@array_loc,\@array_ind,$out_file);
}
###### #################
sub get_best_marker{
  my $ref_array_markers=shift;
  my $ref_hash_marker_qual=shift;
  my $ref_hash_unused_markers=shift;

  ### remove elements which are already used (not in %unused_markers)
  print "\n\n....Selecting Best marker from:",join ":",@{$ref_array_markers} if $verbose;
  print "\nNo markers left to search for contig..." if ($verbose && scalar@{$ref_array_markers} == 0);
  return undef if scalar@{$ref_array_markers} == 0;


  for (reverse(0 .. scalar@{$ref_array_markers}-1)){
    if (!$$ref_hash_unused_markers{$$ref_array_markers[$_]}){
    print "\n\t\t........$_ marker $$ref_array_markers[$_] is already been used. Removing from list." if $verbose;
    splice(@$ref_array_markers,$_,1);
    }
  }
  ### return undef if no markers left in the list.
  print "\nNo markers left to search for contig..." if ($verbose && scalar@{$ref_array_markers} == 0);
  return undef if scalar@{$ref_array_markers} == 0;

  ### find best one in unused set.
  my $best_mark=$$ref_array_markers[0];
  foreach my$mar(@{$ref_array_markers}){
    $best_mark=$mar if ($$ref_hash_unused_markers{$mar}{qual} && $$ref_hash_marker_qual{$mar}{qual} && $$ref_hash_marker_qual{$mar}{qual} < $$ref_hash_marker_qual{$best_mark}{qual} ) ;
  }

  delete $$ref_hash_unused_markers{$best_mark} if $$ref_hash_unused_markers{$best_mark}; ## remove marker from unused list hash
  return($best_mark);
}

sub write_joinmap{
  my$name=shift;
  my$pop_type=shift;
  my$refarray_loc=shift;
  my$refarray_ind=shift;
  my$outfile=shift;

  $name||=$loc_file;
  $pop_type||="CP";


  open(LOCOUT,">$outfile") or die "Unable to open";
  print LOCOUT "name=$name\n";
  print LOCOUT "popt=$pop_type\n";
  print LOCOUT "nloc=",scalar@$refarray_loc,"\n";
  print LOCOUT "nind=",scalar@$refarray_ind,"\n";
  print LOCOUT "\n";

  print LOCOUT join "\n",@$refarray_loc;
  print LOCOUT "\n";
  print LOCOUT "\nindividual names:\n";
  print LOCOUT join "\n",@$refarray_ind;


}

sub get_chi{
  my$mtype1= shift;
  my$mtype2= shift;
  my$ref_array=shift;

  my%count;
  foreach(@$ref_array){
    $count{$mtype1}||=1;
    $count{$mtype1}++ if m/$mtype1/i;
    $count{$mtype2}||=1;
    $count{$mtype2}++ if m/$mtype2/i;
  }
  my$exp=($count{$mtype1} + $count{$mtype2})/2;


  my$chisq=
            (abs($count{$mtype1} -  $exp))**2 / $count{$mtype1}
             +
            (abs($count{$mtype2} -  $exp))**2 /$count{$mtype2}
            ;
  return $chisq;


}