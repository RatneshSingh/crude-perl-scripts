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

## Wrapper to simulate scan_for_matches.

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


my $residue_type = $ARGV[0] || die "Must give residue type";
my $pattern = $ARGV[1] || die "Must give pattern";;
my $sequence_file = $ARGV[2] || die "Must give sequence file name";
my $mismatches = ($ARGV[3] || 0);
my $mismatch_types = ($ARGV[4]  || 'ids');

# We set the maximum buffer size of any particlar sequence to be 10
# megabytes.  You may need to change this if your sequences are
# longer, to avoid buffer splits.
my $MAX_BUFFER_SIZE = 10000000;

## The utility executables are assumed to be in subdirectories relative to this program.
my $NRGREP_PATTERN_CONVERTER_BIN = "perl " . $FindBin::Bin . '/patmatch_to_nrgrep.pl';
my $pattern_for_nrgrep = `$NRGREP_PATTERN_CONVERTER_BIN $residue_type '$pattern'`;
my $NRGREP_BIN = $FindBin::Bin . "/nrgrep_coords/nrgrep_coords";
my $NRGREP_PIPELINE = "$NRGREP_BIN -i -b $MAX_BUFFER_SIZE -k $mismatches$mismatch_types '$pattern_for_nrgrep' '$sequence_file'";

my $SEQUENCE_INDEX_BIN = "perl " . $FindBin::Bin . "/generate_sequence_index.pl";


# Generate the index (should this be already generated into a file?)
debug("Starting to generate index");
open(INDEX_OUTPUT, "$SEQUENCE_INDEX_BIN < $sequence_file |");
my @file_offsets; # array of file byte offsets
my %loci; # hash table with file byte offset as key and locus name as value
debug("Generating index");
while (my $line = <INDEX_OUTPUT>) {
    my @fields = split(/\s/, $line);
    push(@file_offsets, $fields[0]);
    $loci{$fields[0]} = $fields[1];
}
close(INDEX_OUTPUT);


## Print out the results.
debug ("About to call: $NRGREP_PIPELINE");
open(NRGREP_OUTPUT, "$NRGREP_PIPELINE |");
debug ("done");
print_results($residue_type);
close(NRGREP_OUTPUT);


if ($residue_type eq '-c') # now we need to search the forward strand as well
{    
    my $forward_nrgrep_pattern = `$NRGREP_PATTERN_CONVERTER_BIN -n '$pattern'`;
    my $forward_pipeline = "$NRGREP_BIN -b $MAX_BUFFER_SIZE -k $mismatches$mismatch_types '$forward_nrgrep_pattern' '$sequence_file'";
    debug ("About to call $forward_pipeline");
    open(NRGREP_OUTPUT, "$forward_pipeline |");
    debug ("done");
    my $residue = '-n';
    print_results($residue);
    close(NRGREP_OUTPUT);
}

exit 0;



sub print_results
{
    debug("print_results");
    my $residue = shift;
    while (my $line = <NRGREP_OUTPUT>)
    {
	if ($line =~ /^\[/)
	{
	    $line =~ s/[\[\]:,]//g;
	    my @fields = split(/\s+/, $line);
	    my $offset = get_file_offset($fields[0]);

	    my $seq_beg = $fields[0] - $offset + 1;  ## 1-based indicing
	    my $seq_end = $fields[1] - $offset;
	    my $locus_name = $loci{$offset};

	    if (looksLikeFastaHeader($locus_name)) {
		next;
	    }

	    if ($residue eq '-c')
	    {
		my $rev_complement = get_reverse_complement($fields[2]);
		print ">$locus_name:[$seq_end,$seq_beg]\n$rev_complement \n"; # trailing space is intentional
	    }
	    else
	    {
		print ">$locus_name:[$seq_beg,$seq_end]\n$fields[2] \n";  ## trailing space is intentional to keep compatible with scan_for_matches output.
	    }
	}
    }
}


## Returns true if the locus looks like a FASTA header hit.
sub looksLikeFastaHeader {
    my ($locus_name) = @_;
    if ($locus_name =~ m/^>/) {
	return 1;
    }
    return 0;
}



## Test routine to get the offset of the beginning of the sequence.
## inputs: $file_offset
## and global variable @file_offsets.
sub get_file_offset_linear 
{
    debug("get_file_offset_linear");
    my $file_offset = shift;
    my $last_offset = $file_offsets[0];
    for my $offset (@file_offsets) {
	if ($offset > $file_offset) {
	    return $last_offset;
	}
	$last_offset = $offset;
    }
    return $last_offset;
}



# Get the file offset of the beginning of the sequence that contains the given
# file offset
sub get_file_offset 
{
    debug("get_file_offset");
    my $file_offset = shift;

    # print "Searching for $file_offset\n";
    my $array_length = @file_offsets;

    # Perform binary search
    my $low = 0;
    my $high = $array_length - 1;
    while ($high > $low)
    {
	## precondition: $file_offset is in @file_offsets[$low .. $high]

	my $middle = int(($low + $high) / 2);
	if ($file_offsets[$middle] == $file_offset)
	{
	    return $file_offset;
	}
	## handle degenerate case of two elements
	elsif ($high - $low == 1) { 
	    if ($file_offset >= $file_offsets[$high]) {
		return $file_offsets[$high];
	    } else {
		return $file_offsets[$low]; 
	    }
	}

	elsif ($file_offsets[$middle] < $file_offset)
	{
	    $low = $middle;
	}
	elsif ($file_offsets[$middle] > $file_offset)
	{
	    $high = $middle - 1;
	}
	## postcondition: the interval is smaller.  Do we preserve the
	## invariant?  Yes.
    }
    return $file_offsets[$low];
}
    
# Returns the reverse complement of the given nucleotide sequence
sub get_reverse_complement
{
    debug("get_reverse_complement");
    my $sequence = shift; # a nucleotide sequence
    $sequence = reverse($sequence);
    $sequence = complement($sequence);
    return $sequence;
}


# Get the complement of a nucleotide sequence
sub complement
{
    debug("complement");
    my $sequence = shift;
    $sequence =~ tr/atcgATCG/tagcTAGC/;
    return $sequence;
}

sub debug 
{
    my $msg = shift;
#    print STDERR $msg, "\n";
}
