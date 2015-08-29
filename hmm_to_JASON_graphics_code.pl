#!/usr/bin/perl -w
use strict;
use Getopt::Long;

our($hmm_file,$out,$feval,$ceval,$ieval,$num_domain,$name_domain,$complex,$help);

GetOptions( 'f|hmm:s' => \$hmm_file,
           'o|out:s'=>\$out,
		   'fe|feval:f'=>\$feval,
           'ce|ceval:f'=>\$ceval,
           'ie|ieval:f'=>\$ieval,
           'num|num_domain:n'=>\$num_domain,
           'nam|name_domains:s'=>\$name_domain,
           'complex'=>\$complex,
           'help|h'=>\$help
           );


usage() if $help ;
exit if $help;

$feval=$feval?$feval:10;
$ceval=$ceval?$ceval:10;
$ieval=$ieval?$ieval:10;
$num_domain=$num_domain?$num_domain:4;
$name_domain=$name_domain?$name_domain:"";

my(%all_seqs);
my%color;
open HMM,"$hmm_file" or die "\n\n*******Error:cannon open hmm file:$hmm_file**********\n\n";
while(<HMM>){
  next if $_=~/^\#/;
  my@lines=split /\s+/,$_;
  s/\s+//g foreach @lines;
  next if $lines[6] > $feval;
  next if $lines[11] > $ceval;
  next if $lines[12] > $ieval;
  next if $lines[0] !~ /$name_domain/i;
  #print "feval:$feval\tceval:$ceval\tieval$ieval\n";
  #$domains{$lines[3]}{'target'}=$lines[0];
  #$domains{$lines[3]}{'tacc'}=$lines[1];
  #$domains{$lines[3]}{'tlen'}=$lines[2];
  #$domains{$lines[3]}{'qname'}=$lines[3];
  #$domains{$lines[3]}{'qacc'}=$lines[4];
  #$domains{$lines[3]}{'qlen'}=$lines[5];
  #$domains{$lines[3]}{'fs_eval'}=$lines[6];
  #$domains{$lines[3]}{'fs_score'}=$lines[7];
  #$domains{$lines[3]}{'fs_bias'}=$lines[8];
  #$domains{$lines[3]}{'td_num'}=$lines[9];
  #$domains{$lines[3]}{'td_numtotal'}=$lines[10];
  #$domains{$lines[3]}{'td_ceval'}=$lines[11];
  #$domains{$lines[3]}{'td_ieval'}=$lines[12];
  #$domains{$lines[3]}{'td_score'}=$lines[13];
  #$domains{$lines[3]}{'td_bias'}=$lines[14];
  #$domains{$lines[3]}{'hmm_from'}=$lines[15];
  #$domains{$lines[3]}{'hmm_to'}=$lines[16];
  #$domains{$lines[3]}{'ali_from'}=$lines[17];
  #$domains{$lines[3]}{'ali_to'}=$lines[18];
  #$domains{$lines[3]}{'env_from'}=$lines[19];
  #$domains{$lines[3]}{'env_to'}=$lines[20];
  #$domains{$lines[3]}{'dom_reliability'}=$lines[21];
  #$domains{$lines[3]}{'description'}=join " ", @lines[22..$#lines];

   my$json_dom=create_domain(\@lines,$complex);
   push(@{$all_seqs{$lines[3]}{'doms'}},$json_dom);
   $all_seqs{$lines[3]}{'len'}=$lines[5];

}


#### Print the domain information
open OUT,">$out" if $out;
my$fh = *STDOUT;
$fh = *OUT if $out;
foreach my$molecule(keys %all_seqs){

  my$i=0;
  print $fh "\{
    \"length\":\"$all_seqs{$molecule}{'len'}\",
    \"regions\":\[";
  foreach my$doms(@{$all_seqs{$molecule}{'doms'}}){

    print $fh $doms;
    print $fh "," if $i != @{$all_seqs{$molecule}{'doms'}}-1;

    $i++;

  }

  print $fh "\]
  \}






  ";

}







################################
#### subroutine for random color generation
sub random_color {
    my ($r, $g, $b) = (int(rand(256)), int(rand(256)), int(rand(256)));
    my $color1 = sprintf("#%02x%02x%02x", $r, $g, $b);
    my $color2 = sprintf("#%02x%02x%02x", ($r + 128) % 256, ($g + 128) % 256, ($b + 128) % 256);
    return ($color1, $color2);
}

sub fix_color{
  my$domain=shift;
  return $color{$domain} if exists $color{$domain};
  $color{$domain}=random_color();
  #return ["#9999ff","#399","#1fc01f","#c00f0f","#f00","#0f0","#00f","#0ff","#f0f","#ff0"];
  return $color{$domain};
}

sub start_style{
  my $hmmstart=shift;
  return 'curved' if $hmmstart <=1;
  return 'jagged' if $hmmstart>1;
}

sub end_style{
  my $tlen=shift;
  my $hmmend=shift;

  return 'curved' if $hmmend>=$tlen;
  return 'jagged' if $hmmend<$tlen;


}

sub get_family{
  my$acc=shift;
  ###print "\told family:$acc";
  $acc=~s/(\.\d+)//g;
  ####print "\tnew family:$acc\n";
  return $acc;

  }

sub create_domain{
  my$ref_array=shift;
  my$complex=shift;
  my@lines=@{$ref_array};

  my$description=join " ", @lines[22..$#lines];

  #### create json script to draw this domain
  my$json_dom="
      \{
      \"type\" : \"pfama\",
      \"text\" : \"$lines[0]\",
      \"colour\" : \"".fix_color($lines[0])."\",
      \"display\": \"true\",
      \"startStyle\" : \"".start_style($lines[15])."\",
      \"endStyle\" : \"".end_style($lines[2],$lines[16])."\",
      \"start\" : \"$lines[19]\",
      \"end\" : \"$lines[20]\",
      \"aliEnd\" : \"$lines[18]\",
      \"aliStart\" : \"$lines[17]\"

      \}
  " ;

  ### complex models
  $json_dom="

  \{
      \"modelStart\" : \"$lines[19]\",
      \"modelEnd\" : \"$lines[20]\",
      \"colour\" : \"".fix_color($lines[0])."\",
      \"endStyle\" : \"".end_style($lines[2],$lines[16])."\",
      \"startStyle\" : \"".start_style($lines[15])."\",
      \"display\" : true,
      \"end\" : \"$lines[20]\",
      \"aliEnd\" : \"$lines[18]\",
      \"href\":\"/family/".get_family($lines[1])."\",
      \"text\" : \"$lines[0]\",
      \"modelLength\" : \"$lines[2]\",
      \"metadata\" : \{
        \"scoreName\" : \"e-value\",
        \"score\" : \"$lines[11]\",
        \"description\" : \"$description\",
        \"accession\" : \"".get_family($lines[1])."\",
        \"end\" : \"$lines[20]\",
        \"database\" : \"pfam\",
        \"aliEnd\" : \"$lines[18]\",
        \"identifier\" : \"$lines[0]\",
        \"type\" : \"Domain\",
        \"aliStart\" : \"$lines[17]\",
        \"start\" : \"$lines[19]\"
      },
      \"type\" : \"pfama\",
      \"aliStart\" : \"$lines[17]\",
      \"start\" : \"$lines[19]\"
    }
  " if $complex ;

  return $json_dom;

 }


sub usage{

  print "
  This script reads the hmm perdomain output and generates a script to draw graphics using
  http://pfam.janelia.org/generate_graphic
  server.

  usage: script_name options
  options:
  '-f|-hmm:s'         hmm per domain output file from hmmpfam,
  '-o|-out'           save json script in outfile
  '-fe|-feval:f'      filter domains by  full sequence evalue,
  '-ce|-ceval:f'      filter domains by  conditional evalue,
  '-ie|-ieval:f'      filter domains by  independent evalue
  '-num|-num_domain:n'  number of domains to print. Only top doains are printed,
  '-nam|-name_domains' filter domains with names.
  '-complex'         complex format of graphics,
  '-help|h'          help

  ";






}