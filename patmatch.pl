#!/usr/bin/perl -w
use strict;

use FindBin;

## PatMatch
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


## Used for running PatMatch as a stand-alone program.  This is different from
## runPatmatch.pl which is used by the CGI version of PatMatch.  This program
## is a simple pipeline that runs patmatchPatternChecker.pl to check the
## syntax of a pattern, then runs scan_pipeline.pl to find matches.
##
## Arguments:
##
##     [residue_type] [pattern] [sequence_file] [number of mismatches] [mismatch types]
##
## residue_type can be:
##
##    -n: nucleotide
##    -c: complementary (crick) strand
##    -p: protein
##
## mismatch types can be:
##
##    ids: insertions, deletions, substitutions
##    i:   insertions only
##    id:  insertions and deletions
##    etc...
##
## Author: Thomas Yan
## Date: 07/01/2005


my $residue_type = $ARGV[0] || die "Must give residue type";
my $pattern = $ARGV[1] || die "Must give pattern";;
my $sequence_file = $ARGV[2] || die "Must give sequence file name";
my $mismatches = ($ARGV[3] || 0);
my $mismatch_types = ($ARGV[4]  || 'ids');

my $class = '';
if (($residue_type eq "-n") || ($residue_type eq "-c")) { # Nucleotide pattern
    $class = "dna";
} else { # Peptide pattern
    $class = "pep";
}

my $SYNTAX_CHECKER_BIN = "perl " . $FindBin::Bin . 
    '/patmatchPatternChecker.pl';

my $patStatus = `$SYNTAX_CHECKER_BIN $class '$pattern'`;
chomp($patStatus);
if ($patStatus eq "OK") { # syntax OK, run PatMatch
    my $SCAN_PIPELINE = "perl " . $FindBin::Bin . '/scan_pipeline.pl';
    open(OUTPUT, "$SCAN_PIPELINE '$residue_type' '$pattern' '$sequence_file' '$mismatches' '$mismatch_types' |");
    while (my $line = <OUTPUT>) {
	print $line;
    }
    close(OUTPUT);
}
else { # pattern syntax error
    print $patStatus;
    exit();
}
