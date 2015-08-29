#!/usr/bin/perl -w
use strict;

print "\n\n\n: Usage ==> perl script\tname\tlist_of_folder_name\n\n";

open LIST,$ARGV[0];

while (<LIST>){
my$folder=$_;
$folder=~s/\s//g;
`mkdir $folder`;
$folder=();
}
close LIST;
exit;