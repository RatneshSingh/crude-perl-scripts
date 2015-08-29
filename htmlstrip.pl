#!/usr/bin/perl -w
 
use strict;

package HTMLStrip;
use base "HTML::Parser";
 
 my $p = new HTML::Parser;
    $p->parse_file('http://www1.pasteur.fr/cgi-bin/pmtg/fasta.pl?form=html&cmd=FST&attri_1=100.10103.1');

 
 
  sub text {
     my ($self, $text) = @_;
      print $text;
  }
  
  my$ps = new HTMLStrip;
  # parse line-by-line, rather than the whole file at once
 while (<>) {
     $ps->parse($_);
  }
     # flush and parse remaining unparsed HTML
      $ps->eof;