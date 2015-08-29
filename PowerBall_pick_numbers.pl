#!/usr/bin/perl -w
use strict;
use Getopt::Long;


our($num_ticket,$win_num,$win_pb,$max_num,$max_PowerBall,$tp_price,$check_drawing,$generate_number,$num_digits,$help,$use_weight,$weight_file);

my$result=GetOptions(

  "num_ticket|nt:i"=>\$num_ticket,
  "winning_number|wn:s"=>\$win_num,
  "winning_powerball|wp:i"=>\$win_pb,
  "max_num_in_pool|mn:i"=>\$max_num,
  "max_num_in_ball|mb:i"=>\$max_PowerBall,
  "ticket_price|tp:f"=>\$tp_price,
  "check_drawing|cd"=>\$check_drawing,
  "generate_number|gn" =>\$generate_number,
  "num_digits|nd:i"=>\$num_digits,
  "use_weight|uw"=>\$use_weight,
  "weight_file|freq_file|ff|wf=s"=>\$weight_file,
  "help|h"=>\$help
);


my$usage="
###############################################################################################
 perl $0 options....
 OPTIONS:
  num_ticket|nt        Number of random tickets to generate
  winning_number|wn    winning number. (Seperate numbers with :)
  winning_powerball|wp winning powerball,
  max_num_in_pool|mn   Max number in pool to choose from[59]
  max_num_in_ball|mb   Max ball to choose from[34]
  ticket_price|tp      Ticket price to purchase[2]
  check_drawing|cd     Only check the drawing and winnings. dont print ticktets[default].
  generate_number|gn   Print ticket numbers generated randomly.
  num_digits|nd        length of ticket number[5]
                       .e.g Choose 5 number from pool of 1..75.
  weighted|uw          use freq list to pick random numbers [inbuilt freq list].
  wf                   Use provided freq list to to determine weights of each number picked [inbuilt freq list].
                       eg. List of frequencies of numbers picked in last couple years.
                       three columns:Number Number_Freq PowerBall_Freq
                       1\t4\t2
                       2\t18\t3
                       ..
                       ....
                       35\t13\t
  help|h               Print help and exit.
###############################################################################################
";

die "$usage" if defined $help;







$num_ticket=$num_ticket?$num_ticket:10;
$tp_price=$tp_price?$tp_price:2;
$max_num=$max_num?$max_num:59;
$max_PowerBall=$max_PowerBall?$max_PowerBall:34;
$num_digits=$num_digits?$num_digits:5;

$num_ticket=$num_ticket?$num_ticket:10;
my%batch;
my $big_prizes=0;
my%weight;
my%pbweights;
##### create weighted hash for number in pool.
if ($use_weight && $weight_file){
  open(WEIGHT,"$weight_file");
  while (<WEIGHT>) {
    s/^\s+//g;
    my($num,$prob,$pbweight)=split(/\s+/,$_);
    $weight{$num}=$prob;
    $pbweights{$num}=$pbweight if $pbweight !~/^\s*$/;
  }
}

elsif($use_weight && !$weight_file){
  my@PNFreqList=(1,27,2,29,3,29,4,23,5,33,6,25,7,33,8,33,9,30,10,31,11,31,12,26,13,31,14,34,15,24,16,24,17,28,18,22,19,31,20,21,21,17,22,26,23,34,24,27,25,28,26,32,27,18,28,34,29,27,30,24,31,27,32,25,33,26,34,27,35,22,36,32,37,22,38,23,39,31,40,26,41,24,42,23,43,23,44,28,45,29,46,26,47,21,48,28,49,34,50,21,51,25,52,26,53,25,54,36,55,33,56,30,57,26,58,28,59,31);
  my@PBFreqList=(1,11,2,7,3,8,4,6,5,10,6,8,7,9,8,4,9,7,10,10,11,8,12,7,13,11,14,11,15,7,16,12,17,11,18,11,19,9,20,13,21,5,22,9,23,11,24,9,25,8,26,9,27,10,28,10,29,18,30,5,31,4,32,10,33,13,34,8,35,13);


  %weight=@PNFreqList;
  %pbweights=@PBFreqList;
}
else{
  foreach(1 .. $max_num){$weight{$_}=1}
  foreach (1 .. $max_PowerBall){$pbweights{$_}=1}
}




#### pick numbers for tickets
for (1..$num_ticket){
  my @ticket;
  my%temp_hash=%weight;

  for (1..$num_digits){
    my$number=weighted_rand(%temp_hash);
    push(@ticket,$number);
    delete $temp_hash{$number};
  }

   $batch{$_}{'n'}=join ":",sort {$a<=>$b}@ticket;
   $batch{$_}{'p'}=weighted_rand(%pbweights);
}


print "\n\nTickets:" if $generate_number;
foreach my$ticks(keys %batch){
  #print "\nAll values:Tick$ticks:$batch{$ticks}{'n'}\tPB:$batch{$ticks}{'p'}\tWT:$win_num\tWPB:$win_pb\t";
  $big_prizes+=count_prize($batch{$ticks}{'n'},$batch{$ticks}{'p'},$win_num,$win_pb) if $win_num;
  print "\n$batch{$ticks}{'n'}\t$batch{$ticks}{'p'}" if $generate_number;
}
 print "\n";

print "\nTotal money won:$big_prizes\n";
print "\nTotal money Spent:",$num_ticket*$tp_price,"\n";
print "$usage" if !$win_num;










sub count_prize{
  my $ticket_number=shift;
  my $power_number=shift;
  my $winning_number=shift;
  my $winning_PB=shift;

  my$match=0;
  my$pmatch=0;
  my $prize=0;

  foreach my$tnum(split /\:/,$ticket_number){
    foreach my$wnum(split /\:/,$winning_number){
      $tnum=~s/\D+//g;
      $wnum=~s/\D+//g;
      $match++ if $tnum==$wnum;
    }
  }

     $prize=4 if ($match==1 && $power_number==$winning_PB);
     $prize=4 if ($match==0 && $power_number==$winning_PB);
     $prize=7 if ($match==2 && $power_number==$winning_PB);
     $prize=7 if ($match==3 && $power_number!=$winning_PB);
     $prize=100 if ($match==3 && $power_number==$winning_PB);
     $prize=100 if ($match==4 && $power_number!=$winning_PB);
     $prize=10000 if ($match==4 && $power_number==$winning_PB);
     $prize=1000000 if ($match==5 && $power_number!=$winning_PB);
     $prize="9999999999999999" if ($match==5 && $power_number==$winning_PB);

     print "\n*****************\nYou have the winning number\n*****************\n$ticket_number\t$power_number\n*****************\n" if ($match==5 && $power_number==$winning_PB);
     #print "\nLucky Ticket:$ticket_number $power_number\t$winning_number  $winning_PB\tPrize:$prize" if $prize >0;

return ($prize);



}

###########################
#from Perl Cookbook 2.10
#take in a hash and via a weighted random, pick one of the
#keys based on the value and return it
###########################
sub weighted_rand {
    my %dist = @_;
    my $rand;

    my @bucket;
    foreach my $key (keys %dist)
    {
        push @bucket, ($key) x $dist{$key};
    }

    choose_weighted(\@bucket);
}


sub choose_weighted
{
        my $bucket = shift;
        return $bucket->[rand(@$bucket)];
}