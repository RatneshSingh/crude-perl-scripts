# perl
#use 'pl2bat file.pl' command to convert perl file to a .bat file
use strict;
# use LWP::Simple;
use LWP::UserAgent;

my $ua = new LWP::UserAgent;
$ua->timeout(120);
system (cls);

open LIST,"$ARGV[0]";

while(<LIST>){
	s/\s+//g;
	next if $_ eq '';
	my $url1='http://www.dna.affrc.go.jp/sigscan/disp.cgi?';
	
	$url1.=$_;
	my $line=();
	
	$line=fetch_descrition($url1);
	
	print "$_\t$line\n";
}



exit;




#****************************************************************
sub fetch_descrition{
	my $url= shift(@_);
	my $request = new HTTP::Request('GET', $url);
	my $response = $ua->request($request);
	my $content = $response->content();
#	print $content;

	my@page= split(/\n/,$content);
	my $line_req= my$line1=();
	
	foreach (@page){
                                               #print "$line\n\n\n";
		if($_=~/^DE\s+([\s\w\W\d]+)/){
		                                       #print "found line $line\n";
			$line1=$1;
		                                       #print "This is line: $line1\n";
		}
	$line_req.=$line1;
	$line1=();
	
	}
	$line_req =~ s/\s+/ /g;
	
	                                           #print "This is joined line: $line_req\n";
	return($line_req);
}