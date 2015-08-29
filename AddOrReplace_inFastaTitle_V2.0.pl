#!/usr/bin/perl -w
use strict;
use Getopt::Std;


#This script is made to modify Fasta Header as per user requirement.
# for Bugs contact
# Author: Ratnesh.Singh@ag.tamu.edu

our ( $inputfile, $outputfile, %seq_hash, @pattern, @notfoundList, @foundList, $manual_pattern, $manual_replacement );
our (%pattern_key_hash);
our ( $opt_s, $opt_t, $opt_d, $opt_c, $opt_o, $opt_m, $opt_n, $opt_p, $opt_r, $opt_i, $opt_v );

$opt_d = "_";
$opt_m = "exact";
$opt_c = "a";
$opt_t = "length";
$opt_n = "no";
$opt_i = "no";
$opt_v = "no";

getopt('stdconmpriv');

my $usage = "\n
This script will modify the fasta headers in a fasta file based on the following options selected.

\n\nusage : perl script_name [options]

options:
-s	Sequence file containing sequences in fasta format
-t	table_file_name|length|manual. [length]
	table_file_name: table of two columns, 1:pattern to look for, 2: New name to be used for modification.
			The pattern and new name should be seperated by tab, Pattern in first row.
			e.g. pattern1	New_name1
			pattern2	New_name2
	length: Use length of the sequence as New_name.
	manual: provide pattern with '-p' and replacement string with '-r' flags. If not, will be asked for pattern and
	replacement string manually. Only append/prepend/substitutes the pattern with the replacement string in the header.
-d	delimiter to use to append or prepend the New_name ['_']
-c	a|r|p|s. [default is append]
	[a]ppend:add to the end,
	[r]eplace: replace the whole header with New_Name,
	[p]repend: add to the start.
	[s]ubstitute: only substitute the patern with with New_Name.
-o	Output file to store resulting sequences.[Mod_inputFileName]
-n	Yes|no. Print list of patterns not found in sequence.
-m	Search mode to be used. Use following options:[exact]
	'match' to use pattern matching mode.
	'exact' to use exact match mode.
-p	Manual pattern to look for. provide only one string. works with manual mode only.
-r	Manual replacement string to replace/append/prepend when manual pattern is found. works with manual mode only.
-i	yes|no. ignore case while pattern match [no].
-v	yes|no. Print the list of all the replacements done.[no]
";

#opening inputfile containing sequences in fasta format.
if ($opt_s) { $inputfile = $opt_s; }
else        { die "\nPlease provide sequence filename.\n $usage"; }
chomp($inputfile);
open INPUT, "<$inputfile" or die "Cannot open $inputfile.....\n\n";

my $filename = getfilename($opt_s);

#opening outputfile
if   ($opt_o) { $outputfile = $opt_o; }
else          { $outputfile = 'Mod_' . $filename; }
chomp($outputfile);
open OUT, ">$outputfile" or die "cannot create $outputfile.....\n\n";

#opening file containing conversion list pattern and new name seperated by white spaces.
if ( lc $opt_t ne 'length' && lc $opt_t ne "manual" ) {
    open TABLE, "$opt_t";
    while (<TABLE>) {
        next if /^\s*$/;
        my ( $pattern, $mod_key ) = split( /\s+/, $_ );

        $pattern_key_hash{$pattern} = $mod_key;
        push( @pattern, $pattern );

    }
}
elsif ( lc $opt_t eq "manual" ) {
    if ($opt_p) { $manual_pattern = $opt_p; }
    else        { print "\nType the pattern to look for:"; $manual_pattern = <STDIN>; }
    chomp $manual_pattern;

    if ($opt_r) { $manual_replacement = $opt_r; }
    else        { print "Type the string to replace manual pattern with:"; $manual_replacement = <STDIN>; }
    chomp $manual_replacement;
}

#*************************************************************************************************************************
#to check if pattern is properly formed and working or not.

my $count = @pattern;    #just to check if array is formed or not

#~ print "Elements in pattern array \n@pattern\n";
print "Done........\nNumber of Pattern in pattern file =$count\n\n\n";

#~ print "$pattern[0],$pattern[1],$pattern[2],$pattern[3],$pattern[4]";

#*************************************************************************************************************************
# Read sequences in hash.
%seq_hash = ReadFasta($inputfile);

#*************************************************************************************************************************
#read all the pattern or key words into array @pattern. and search for the matching one in abovemade hash.
my $count_found = my $count_notfound = 0;
print "Looking for patterns...Plz wait....\n";
##########################################################################################################################
####if ( lc $opt_t ne "length" && lc $opt_t ne "manual" ) {
####    foreach my $ta ( keys %pattern_key_hash ) {
####        next if $ta =~ /^\s*$/;
####        chomp($ta);
####        my ( $header, $sequence, $headerNew ) = ();
####        if   ( $opt_m eq 'match' ) { ( $header, $sequence ) = searchPattern($ta); }
####        else                       { ( $header, $sequence ) = searchExact($ta); }
####        if ( $header && $sequence ) {
####            if    ( $opt_c eq "a" ) { $headerNew = $header . $opt_d . $pattern_key_hash{$ta} }
####            elsif ( $opt_c eq "r" ) { $headerNew = $pattern_key_hash{$ta} }
####            elsif ( $opt_c eq "p" ) { $headerNew = $pattern_key_hash{$ta} . $opt_d . $header }
####            elsif ( $opt_c eq "s" ) { $headerNew = $header; $headerNew =~ s/$ta/$pattern_key_hash{$ta}/; }
####            $seq_hash{$headerNew} = $sequence;
####            delete $seq_hash{$header};
####            $count_found++;
####            push( @foundList, $ta );
####            if ( lc $opt_v eq 'yes' ) { print "Replaced $header \<----- with-----\> $headerNew\n"; }
####        }
####        else {
####            $count_notfound++;
####            push( @notfoundList, $ta );
####        }
####    }
####}
########################################################################################################################
if ( lc $opt_t ne "length" && lc $opt_t ne "manual" ) {

    #foreach my$ta(keys %pattern_key_hash){
    foreach my $header ( keys %seq_hash ) {
        next if $header =~ /^\s*$/;
        chomp($header);
        our ( $ta, $mod_pattern, $sequence, $headerNew );
        if   ( $opt_m eq 'match' ) { ( $ta) = searchAllPattern($header); }
        else                       { ( $ta) = searchAllExact($header); }
        if ( $ta ) {
            if    ( $opt_c eq "a" ) { $headerNew = $header . $opt_d . $pattern_key_hash{$ta} }
            elsif ( $opt_c eq "r" ) { $headerNew = $pattern_key_hash{$ta} }
            elsif ( $opt_c eq "p" ) { $headerNew = $pattern_key_hash{$ta} . $opt_d . $header }
            elsif ( $opt_c eq "s" ) { $headerNew = $header; $headerNew =~ s/$ta/$pattern_key_hash{$ta}/; }
            $sequence=$seq_hash{$header};
			$seq_hash{$headerNew} = $sequence;
            delete $seq_hash{$header};

            $count_found++;
            push( @foundList, $ta );
            if ( lc $opt_v eq 'yes' ) { print "$ta\tReplaced $header \<----- with-----\> $headerNew\n"; }

        }
        else {
            $count_notfound++;
            push( @notfoundList, $header );
        }

    }

}

#########################################################################################################################
# if sequence length has to be added into the header.
elsif ( lc $opt_t eq "length" ) {
    foreach my $header ( keys %seq_hash ) {
        my $seq_length = "SeqLlength:" . length( $seq_hash{$header} );
        my ($headerNew);
        if    ( $opt_c eq "a" ) { $headerNew = $header . $opt_d . $seq_length }
        elsif ( $opt_c eq "r" ) { $headerNew = $seq_length }
        elsif ( $opt_c eq "p" ) { $headerNew = $seq_length . $opt_d . $header }
        else                    { $headerNew = $header . $opt_d . $seq_length }

        $seq_hash{$headerNew} = $seq_hash{$header};
        delete $seq_hash{$header};

        $count_found++;
        push( @foundList, $header );
        if ( lc $opt_v eq 'yes' ) { print "Replaced $header \<----- with-----\> $headerNew\n"; }
    }
}
#########################################################################################################################
# if Manual entry method has to be used.
elsif ( lc $opt_t eq "manual" ) {
    print "Finding and replacing $manual_pattern with $manual_replacement in each sequence headers\n";
    foreach my $header ( keys %seq_hash ) {
        my ($headerNew);

        if ( $opt_c eq "a" ) { $headerNew = $header . $opt_d . $manual_replacement }

        #elsif($opt_c eq "r"){$headerNew=$pattern_key_hash{$ta}} # does not make sense to replace whole header with the string. uncomment if required.
        elsif ( $opt_c eq "p" ) { $headerNew = $manual_replacement . $opt_d . $header }
        else                    { $headerNew = $header; $headerNew =~ s/$manual_pattern/$manual_replacement/ig; }

        $seq_hash{$headerNew} = $seq_hash{$header};
        delete $seq_hash{$header};

        $count_found++;
        push( @foundList, $header );
        if ( lc $opt_v eq 'yes' ) { print "Replaced $header \<----- with-----\> $headerNew\n"; }

    }

}
######################################################################################################################
# Print modified sequences.
foreach ( keys %seq_hash ) { print OUT">$_\n$seq_hash{$_}\n"; }
if ( lc $opt_n ne "no" ) { open OUT2, "NotFound_" . $inputfile . ".list"; print OUT2 join( "\n", @notfoundList ) }

print "\nNumber of Patterns match found for = $count_found\n";
print "Number of Patterns match not found for = $count_notfound\n";
print "Following patterns were nor found in the sequences:" . join( "\t", @notfoundList ) . "\n";

#*************************************************************************************************************************
#exit the program and close all files.
close INPUT;
close OUT;
#close PATTERN;
exit;

#######################################################
#-------------------------------------Start ReadFasta---------------------------------------+

sub ReadFasta {    # to read fasta format files into hash. returns hash.

    my $seqfile = shift(@_);

    my ( $header, @sequence );
    chomp $seqfile;
    open FASTA, "$seqfile";
    print "reading Sequences from input file.....Plz wait...\n";
    my %seq_hash = ();

    #$seq_hash{'RS_Concatenated'}="";

    $/ = "\n>";    # Change record seperator to read Fasta
    my $last_N = 1;
    while (<FASTA>) {
        chomp;
        ( $header, @sequence ) = split( "\n", $_ );

        $header =~ s/>//;       # Remove Leading > from Header
        $header =~ s/\s*$//;    # Remove trailing spaces from header
        $header =~ s/^\s*//;    # Remove Leading spaces from Header

        my $sequence = join( "", @sequence );
        $sequence =~ s/\s//g;
        $sequence =~ s/\n//g;

        if ( $header =~ /^\s*$/ ) { next; }

        # Make sure no duplicate header exists else they will be overwritten and lost during hashing due to same key.
        if ( !exists $seq_hash{$header} ) {
            $seq_hash{$header} = $sequence;    #feed headers and sequences in hash.
                                               #$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
        }
        else {

            # find a uniq header name by adding a number at the end. If header still exists, increase the number by one
            while ( exists $seq_hash{$header} ) { $header = $header . $last_N; $last_N++; }

            $seq_hash{$header} = $sequence;

            #$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;

        }
    }

    my @seq_count = keys(%seq_hash);
    my $seq_count = @seq_count;

    print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
    @seq_count = ();
    $/         = "\n";    # Record seperator set back to default newline.

    return (%seq_hash);

}

#-------------------------------------End ReadFasta---------------------------------------+

#SUBROUTINE TO PRSETHROUGH HASH LOKING FOR PATTERN AND RETURN THE SEQUENCE IF PATTERN MATCHES
sub searchAllPattern {
    my $header = shift;
    chomp($header);
    my $mod_header = cleanString($header);
    foreach my $key ( keys %pattern_key_hash ) {
        my $mod_key = cleanString($key);
        if ( lc $opt_i eq "yes" ) { $mod_key = lc $mod_key; $mod_header = lc $mod_header; }
        if ( $mod_header =~ m/$mod_key\D/ ) {return ($key);}
        next;
    }
}
########################################################################################
sub searchPattern {
    my $pattern = shift;
    chomp($pattern);
    my $mod_pattern = cleanString($pattern);
    foreach my $key ( keys %seq_hash ) {
        my $mod_key = cleanString($key);
        if ( lc $opt_i eq "yes" ) { $mod_key = lc $mod_key; $mod_pattern = lc $mod_pattern; }
        if ( $mod_key =~ m/$mod_pattern\D/ ) {    #use this line for matching pattern

            # 			print "For pattern:\t$pattern\t\tPicked seq-->\t$key\n";
            return ( $key, $seq_hash{$key} );
        }else{ next;}
    }
}

###########################################################################################
#SUBROUTINE TO PRSETHROUGH HASH LOKING FOR PATTERN AND RETURN THE SEQUENCE IF PATTERN MATCHES
sub searchExact {
    my $pattern = shift;
    chomp($pattern);
    if ( exists $seq_hash{$pattern} ) { return ( $pattern, $seq_hash{$pattern} ); }
    else                              { print "Cannot find $pattern\n "; }
}
###########################################################################################
#SUBROUTINE TO PRSETHROUGH HASH LOKING FOR PATTERN AND RETURN THE SEQUENCE IF PATTERN MATCHES
sub searchAllExact {
    my $header = shift;
    chomp($header);
    if ( exists $pattern_key_hash{$header} ) { return ( $header ); }
    else                              { print "Cannot find $header\n "; return 0; }
}

#########################################################################################
sub getfilename {
    my $rawfilename = shift;

    #remove folder names from rawfilename names.
    my @names = split( /[\/\\]/, $rawfilename );
    return $names[-1];

}


#################################################################################

  sub cleanString {
    my$mod_pattern=shift;
	$mod_pattern=~s/\s*$//g;					# Remove trailing spaces from header
	$mod_pattern=~s/^\s*//g;					# Remove Leading spaces from Header
#	$mod_pattern=~s/([^A-Za-z0-9\s\-\+\=\_])/\\$1/g; # not working. Escapes special characters
	$mod_pattern=~s/([^A-Za-z0-9\s\-\+\=\_])/_/g;
	return($mod_pattern);
}
