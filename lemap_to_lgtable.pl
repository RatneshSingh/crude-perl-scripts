
my$mtype='male';
open NAMES,$ARGV[0];
my$count=0;
while(<NAMES>){
next if m/^\s*$/;
$count++;
s/^\s+//g;
s/\s*$//g;
$mar_name{$count}=$_;
}

open LEMAP, $ARGV[1];
my$lgroup=0;
while(<LEMAP>){
#$lgroup=$1 if m/\#\*\*\*\s*LG\s*=\s*(\d+)\s+/;
if(m/^\#\*\*\*\s*LG\s*=\s*(\d+)\s+/){
$lgroup=$1;
print "\n\n;lowerbound 0.0";
print "\ngroup LG$lgroup";
print "\n;BEGINOFGROUP";
print "\nLocus   Position";

}elsif(m/^#COUNT\s*=\s*\d+/){
print "\n\n;ENDOFGROUP";
}

next if m/^#/;
#print "\n Found Linkage group $lgroup";
#marker_number  male_position   female_position ( error_estimate )[ duplicate* OR phases for each family]
my($mar_num,$male_pos,$female_pos,$error,$phase)=split(/\s+/,$_);
###;lowerbound 0.0
#group LG10
#;BEGINOFGROUP
#Locus   Position
#

print "\n$mar_name{$mar_num}\t",$mtype eq 'male'?$male_pos:$female_pos;




}



