#!/usr/bin/perl
# Author: thomasw
# Date: 2007-10-29 - Added outfile existence check
# Date: 2006-08-15
# Description: Takes in FASTA and CSV with names of ids to extract and
# extracts the FASTAs.

use warnings;
use strict;

use Bio::Seq;
use Bio::SeqIO;
use Text::CSV;

unless (@ARGV) {
  #die &usage;
  &usage;
}

my $in_file = shift @ARGV;
my $csv_file = shift @ARGV;
my $out_file = shift @ARGV;

if (-e $out_file) {
  print "Outfile: $out_file already exists.\n";
  &usage;
}

my %ids_hash;

{ # Input ids to extract
  my $csv = Text::CSV->new();
  open (IDS, "<$csv_file");
  while (my $line = <IDS>) {
    if ($csv->parse ($line)) {
      my @fields = $csv->fields();
      $ids_hash{$fields[0]} = 1;
    }
    else {
      my $err = $csv->error_input();
      die "Error: split() failed to gerate CSV input line: $err\n";
    }
  }
  close (IDS);
}

{
  my $in_sio = Bio::SeqIO->new('-file' => "$in_file", '-format' => 'fasta');
  my $out_sio = Bio::SeqIO->new('-file' => ">$out_file", '-format' => 'fasta');

  while (my $in_seq_obj = $in_sio->next_seq()) {
    my $in_seq_id = $in_seq_obj->primary_id();
    #print "DEBUG: inseqid = $in_seq_id\n";
    if (defined ($ids_hash{$in_seq_id})) {
      my $in_seq = $in_seq_obj->seq();
      my $out_seq = Bio::Seq->new('-id' => $in_seq_id, '-desc' => $in_seq_obj->desc(), '-seq' => $in_seq);
      $out_sio->write_seq($out_seq);
    }
  }
}

sub usage {
  print "Usage: $0 fastaInputFile CSVFileWithIdsToExtract fastaOutputFile\n";
  print "Note: Does not use regex, so \| values are not considered \"or\"s.\n";
  exit;
}
