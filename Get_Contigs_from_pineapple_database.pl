# perl
#use 'pl2bat file.pl' command to convert perl file to a .bat file
use strict;
# use LWP::Simple;
use LWP::UserAgent;

my $ua = new LWP::UserAgent;
$ua->timeout(120);
#system (cls);


my $url='http://genet.imb.uq.edu.au/cgi-bin/Pineapple/single.pl?id=Contig_';


open OUT,">Pineapple_database_contigs.fasta" or die " Cannot open output file\n";

for(my$i=1;$i<=224;$i++){
        my$seq=();
        my$newfinalurl=$url.$i;
          $newfinalurl=~s/\s+//g;
			#print "url:$newfinalurl\n";
                  
           print "$newfinalurl\n";
           $seq=pullseq($newfinalurl);
      		#print "$seq\n";
		print OUT ">Pineapple_Contig_$i\n$seq\n";
		print "******\nSaved sequence for Pineapple_Contig_$i\n$seq\n\n******"
      		

}

print "Sequences were saved in file named Pineapple_database_contigs.fasta\n";

sleep(10);

exit;
#****************************************************************
sub pullseq{
	my $url= shift(@_);
	my $request = new HTTP::Request('GET', $url);
	my $response = $ua->request($request);
	my $content = $response->content();
#        print "content:$content\n" ;
     	if($content=~/Contig sequence<\/td><td><tt>([\w<br>\s]+)<\/tt><\/td><\/tr>/){
          my$seq=$1;
          $seq=~s/<br>//g;
	  $seq=~s/<td>//g;
	  $seq=~s/<tt>//g;
          $seq=~s/\s+//g;
	  $seq=~s/\///g;
          return($seq);
       	  #print "$seq";
		}
        else{return ('noseq');}
}