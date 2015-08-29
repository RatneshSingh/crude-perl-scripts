#!/usr/bin/perl -w

## PatMatch
## unjustify_fasta.pl v. 1.1 (updated to fix 1 record skip bug)
## fix also gets rid of first empty line
## Copyright (C) 2004 The Arabidopsis Informatin Resource
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

use strict;



if (@ARGV == 0) {
    unjustify(\*STDIN, \*STDOUT);
} else {
    for my $filename (@ARGV) {
	open(IN, $filename) || die;
	open(OUT, ">$filename.prepared") || die;
	unjustify(\*IN, \*OUT);
	close(IN);
	close(OUT);
    }
}



sub unjustify {
    my ($INFILE, $OUTFILE) = @_;
    my @buffer;
    while (my $line = <$INFILE>) {
	chomp ($line);
	if ($line =~ /^>/) {
	    if(scalar(@buffer) > 0) {
	   	 print $OUTFILE join("", @buffer);
	    	 print $OUTFILE "\n";
	    }
	    @buffer = ("$line\n");
	}
	else {
	    $line =~ s/[^a-zA-Z]//;
	    push @buffer, $line;
	}
    }
    if(scalar @buffer > 0) {
	    print $OUTFILE join("", @buffer);
	    print $OUTFILE "\n";
    }
}
