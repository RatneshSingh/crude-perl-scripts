# perl
#use 'pl2bat file.pl' command to convert perl file to a .bat file
use strict;
# use LWP::Simple;
use LWP::UserAgent;

my $ua = new LWP::UserAgent;
$ua->timeout(120);
#system (cls);

my $url='http://www1.pasteur.fr/cgi-bin/pmtg/fasta.pl?form=html&cmd=FST&attri_1=';


my ($list2,$output);
if ($ARGV[0]){$list2=$ARGV[0];}
else{print"write the name of List file\n";
	$list2=<STDIN>;
	}
open LIST,"$list2" or die "Cannot find list file\n";

if ($ARGV[1]){$output=$ARGV[1];}
else{print "Print the name of output file \n";
my $output=<STDIN>;
	}

open OUT,">$output" or die " Cannot open output file\n";

open OUT2,">SeqsNotFound.txt";


foreach(<LIST>){

	chomp;
	## generataion of final url
	my$newfinalurl=$_;
	$newfinalurl=~s/ //g;
	#print "url:$newfinalurl\n";
                  
   print "$newfinalurl\n";
  my $seq=pullseq($newfinalurl);
  #print "$seq\n";
      		
 # deciding when to end loop;
 if ($seq ne 'noseq'){print OUT "$seq\n";} 
 else{print "Sequence not found:$newfinalurl\n trying next\n";
      print OUT2 "$newfinalurl\n";
     # last;
      }

      }


print "Sequences were saved in file named $ARGV[1]\n";

sleep(10);

exit;
#****************************************************************
sub pullseq{
	my $url= shift(@_);
	my $request = new HTTP::Request('GET', $url);
	my $response = $ua->request($request);
	my $content = $response->content();
#        print "content:$content\n" ;
     	if($content=~/>(>\d+[\s\S\w\W]+)<br><\/tt>/){
          my$seq=$1;
          $seq=~s/<br>/\n/g;
          $seq=~s/\t+/ /g;
          return($seq);
       	  #print "$seq";
		}
        else{return ('noseq');}
}