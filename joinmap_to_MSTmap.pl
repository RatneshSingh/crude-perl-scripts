#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
our(
$jmp,$outfile,
$population_type,
$population_name ,
$distance_function ,
$cut_off_p_value ,
$no_map_dist ,
$no_map_size ,
$missing_threshold ,
$estimation_before_clustering ,
$detect_bad_data ,
$objective_function ,
$number_of_loci ,
$number_of_individual ,
$info
);


my $usage=join "","
perl $0 options
OPTIONS:
infile |i                             input join map file name.
outfile|o                             outfile name [",$jmp?$jmp:'Inputfilename',".mstmap.output]
population_type|pt                    population_type[DH]
population_name|pn                    population_name[input_file_name]
distance_function|df                  distance_function[kosambi]
cut_off_p_value|cop                   cut_off_p_value[0.0000001]
no_map_dist|npd                       no_map_dist[15]
no_map_size|nps                       no_map_size[2]
missing_threshold|mt                  missing_threshold[0.25]
estimation_before_clustering|ebc      estimation_before_clustering[yes]
detect_bad_data|dbd                   detect_bad_data[yes]
objective_function|of                 objective_function[COUNT]
info                                  use informative naming system for files[off].
                                      If turned On, file names will have parameters values added to.
";


GetOptions(
"infile|i=s"=>\$jmp,
"population_type|pt=s"=>\$population_type,
"population_name|pn=s"=>\$population_name ,
"distance_function|df=s"=>\$distance_function ,
"cut_off_p_value|cop=f"=>\$cut_off_p_value ,
"no_map_dist|npd=f"=>\$no_map_dist ,
"no_map_size|nps=i"=>\$no_map_size ,
"missing_threshold|mt=f"=>\$missing_threshold ,
"estimation_before_clustering|ebc=s"=>\$estimation_before_clustering ,
"detect_bad_data|dbd=s"=>\$detect_bad_data ,
"objective_function|of=s"=>\$objective_function ,
"info"=>\$info
);



$jmp || die "\nPlease provide an input file in joinmap format\n$usage\n";

$population_type||="DH";          #< para1 >
$population_name ||="$jmp";       #< para2 >
$distance_function ||="kosambi";  #< para3 >
$cut_off_p_value ||="0.0000001";   #< para4 >
$no_map_dist ||="15";             #< para5 >
$no_map_size ||="2";              #< para6 >
$missing_threshold ||="0.25";     #< para7 >
$estimation_before_clustering ||="yes";   #< para8 >
$detect_bad_data ||="yes";                #< para9 >
$objective_function ||="COUNT";            #< para10 >

my@ind_name;
my@loci;


open(JMP,$jmp) or die "Unable to open file $jmp";
our($name,$popt,$nloc,$nind,$individual_names,$printname,$spr);
while (<JMP>) {
  next if m/^\s*$/;
  $name=$1 if m/^name\s*=\s*([^\s]+)/;
  $popt=$1 if m/^popt\s*=\s*([^\s]+)/;
  $nloc=$1 if m/^nloc\s*=\s*([^\s]+)/;
  $nind=$1 if m/^nind\s*=\s*([^\s]+)/;
  $individual_names=$1 if m/^(individual\s*names:)/;
  if ($nind && $spr && !$individual_names) {
    s/\<[\w]+\>//g;
    s/\s+/\t/g;
    s/lm/B/gi;
    s/ll/A/gi;
    s/np/B/gi;
    s/nn/A/gi;
    s/--/U/gi;
    push(@loci,$_);
  }
  elsif($individual_names && $printname){my$loci=$_; $loci=~s/\s+//g;push(@ind_name, $loci)}

  elsif($individual_names){$printname=1}
  elsif($nind){$spr=1}
}



my$rand=time;
my$mstfile="$jmp.$rand.mstmap";
$mstfile="$jmp.mstmap.DF$distance_function.COP$cut_off_p_value.NMD$no_map_dist.NMS$no_map_size.MT$missing_threshold.EBC$estimation_before_clustering" if $info;

open(OUT,">$mstfile.input");


print OUT "population_type $population_type
population_name $population_name
distance_function $distance_function
cut_off_p_value $cut_off_p_value
no_map_dist $no_map_dist
no_map_size $no_map_size
missing_threshold $missing_threshold
estimation_before_clustering $estimation_before_clustering
detect_bad_data $detect_bad_data
objective_function $objective_function
number_of_loci ", scalar@loci,"
number_of_individual ", scalar@ind_name,"\n\n";
print OUT join "\t","locus_name",@ind_name;
print OUT "\n";
print OUT join "\n",@loci;

close OUT;
sleep(1);

open(LOG,">$mstfile.log");
print LOG "
\#MSTMap was run with following options.
\# inputfile: $mstfile.input
\# output file:$mstfile.output
\# log file:$mstfile.log
\#population_type $population_type
\#population_name $population_name
\#distance_function $distance_function
\#cut_off_p_value $cut_off_p_value
\#no_map_dist $no_map_dist
\#no_map_size $no_map_size
\#missing_threshold $missing_threshold
\#estimation_before_clustering $estimation_before_clustering
\#detect_bad_data $detect_bad_data
\#objective_function $objective_function
\#number_of_loci ", scalar@loci,"
\#number_of_individual ", scalar@ind_name,"\n
\#Output results from MSTMap run:\n
\n
";
close LOG;
sleep(1);
my$cmd=join " ","MSTMap","$mstfile.input","$mstfile.output",">>","$mstfile.log","2>&1";
print "\nRunning MSTmap with converted input file with following command\n";
system($cmd) || die "\nUnable to run MSTmap command\n";