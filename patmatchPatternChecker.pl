#!/usr/bin/perl -w
use strict;
use Parse::RecDescent;

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

## Checks the syntax of a PatMatch pattern.  This should be used before
## running patmatch_to_nrgrep to convert a PatMatch pattern into an NR-Grep
## pattern.
## 
## Author: Thomas Yan

# message used by error checker
my $MESSAGE = "Please check your pattern and try again.\n";
my $USAGE = "Usage: patmatchPatternChecker.pl <class> <pattern>\n";

my $class = $ARGV[0] || die "$USAGE" . "Valid classes: " .
    "\"dna\" for nucleotides or \"pep\" for peptides.\n";
my $pattern = $ARGV[1] || die "$USAGE" . "Give me a pattern to check.\n";
checkPattern($pattern);

#############################################################################
##
##  SUBROUTINE NAME
##    checkPattern()
##
##  SYNOPSIS 
##    checkPattern()
##
##  DESCRIPTION
##    Performs a test of the validity of the pattern provided by the user.
##
##  ARGUMENTS
##    $pattern - (in) The pattern to check.
##
##  RETURN VALUE
##    none
##
#############################################################################

sub checkPattern
{
    my $pattern = shift;

    # Context free grammar for the Patmatch pattern syntax
    my $grammar = q{
       startrule: pattern /^\Z/
       pattern: left_anchor query right_anchor
       left_anchor: '<' | ''
       right_anchor: '>' | ''
       query: single query | single
       single: literal_list repeat | notany repeat | any repeat | any repeat 
       | group repeat | literal | notany | any | group
       literal: /[A-Za-z\.]/
       literal_list: literal literal_list | literal
       notany: '[^' literal_list ']'
       any: '[' literal_list ']'
       group: '(' query ')'
       repeat: '{' /\d+/ ',' /\d+/ '}'
       | '{' /\d+/ ',' '}'
       | '{' ',' /\d+/ '}'
       | '{' /\d+/ '}'
    };

    my $parser = new Parse::RecDescent($grammar);
    if (defined $parser->startrule($pattern))
    {
	checkCharacters($pattern);
	# checkMinimumLength($pattern);
    }
    else
    {
	errorReport("Invalid pattern syntax.  $MESSAGE");
    }

    # The pattern syntax is correct.
    print "OK\n";

} # checkPattern()

#############################################################################
##
##  SUBROUTINE NAME
##    checkCharacters()
##
##  SYNOPSIS 
##    checkCharacters()
##
##  DESCRIPTION
##    Check for invalid characters in the pattern
##
##  ARGUMENTS
##    $pattern - (in) The pattern to check.
##
##  RETURN VALUE
##    none
##
#############################################################################

sub checkCharacters
{
    my $pattern = shift;
    
    if ($class eq "dna")
    {
	checkNucleotides($pattern);
    }
    else
    {
	checkPeptides($pattern);
    }
} # checkCharacters()

#############################################################################
##
##  SUBROUTINE NAME
##    checkNucleotides()
##
##  SYNOPSIS 
##    checkNucleotides()
##
##  DESCRIPTION
##    Check for invalid characters in the pattern
##
##  ARGUMENTS
##    $pattern - (in) The pattern to check.
##
##  RETURN VALUE
##    none
##
#############################################################################

sub checkNucleotides
{
    my $pattern = shift;

    if ($pattern =~ m/[EFIJLOPQZefijlopqz]/)
    {
	errorReport("Invalid nucleotide character found in pattern.  " .
		    "$MESSAGE");
    }
} # checkNucleotides()

#############################################################################
##
##  SUBROUTINE NAME
##    checkPeptides()
##
##  SYNOPSIS 
##    checkPeptides()
##
##  DESCRIPTION
##    Check for invalid peptide characters in the pattern
##
##  ARGUMENTS
##    $pattern - (in) The pattern to check.
##
##  RETURN VALUE
##    none
##
#############################################################################

sub checkPeptides
{
    my $pattern = shift;
    
    if ($pattern =~ m/[uU]/)
    {
	errorReport("Invalid peptide character found in pattern.  " .
		    "$MESSAGE");
    }
} # checkPeptides()

#############################################################################
##
##  SUBROUTINE NAME
##    checkMinimumLength()
##
##  SYNOPSIS 
##    checkMinimumLength()
##
##  DESCRIPTION
##    Check to make sure that the pattern is not shorter than 3 tokens.
##    Anything between { }, ( ), or [ ] counts as one token.
##    Assumes that the pattern syntax is correct.
##
##  ARGUMENTS
##    $pattern - (in) The pattern to check.
##
##  RETURN VALUE
##    none
##
#############################################################################

sub checkMinimumLength
{
    my $pattern = shift;
    
    my @patArray = split(//, $pattern);
    my $tokens = 0;
    my $countingMode = 1; # indicates if characters are being counted or not
    my $prevCharOpenBracket = 0; # previous character was '{'

    my $patLength = @patArray;
    for (my $i = 0; $i < $patLength; $i++)
    {
	my $char = $patArray[$i];
	if ($prevCharOpenBracket)
	{
	    $prevCharOpenBracket = 0; # no longer true
	    $tokens += getIntValueAfterBracket($pattern, $i);
	}
	elsif ($char eq '{')
	{
	    $prevCharOpenBracket = 1;
	    $countingMode = 0;
	}
	elsif ($char eq '[')
	{
	    $tokens = incrementToken($countingMode, $tokens);
	    $countingMode = 0;
	}
	elsif (($char eq '(') || ($char eq ')') || ($char eq ']') 
	       || ($char eq '}'))
	{
	    $countingMode = 1;
	}
	else
	{
	    $tokens = incrementToken($countingMode, $tokens);
	}
    }
    
    if ($tokens < 3)
    {
	errorReport("Your pattern is shorter than the minimum number of 3 "
		    . "residues.  $MESSAGE");
    }
} # checkMinimumLength()

#############################################################################
##
##  SUBROUTINE NAME
##    incrementToken()
##
##  SYNOPSIS 
##    incrementToken()
##
##  DESCRIPTION
##    Check for invalid peptide characters in the pattern
##
##  ARGUMENTS
##    $countingMode - 1 if should increment, 0 if should not increment
##    $tokenValue - the current value of the token
##
##  RETURN VALUE
##    incremented value of the token
##
#############################################################################

sub incrementToken
{
    my $countingMode = shift;
    my $tokenValue = shift;

    if ($countingMode)
    {
	$tokenValue++;
    }
    return $tokenValue;
} # incrementToken()

#############################################################################
##
##  SUBROUTINE NAME
##    getIntValueAfterBracket()
##
##  SYNOPSIS 
##    getIntValueAfterBracket()
##
##  DESCRIPTION
##    Get the value of the integer after the '{' in the pattern whose
##    location is specified by the given index.
##
##  ARGUMENTS
##    $pattern - the PatMatch pattern
##    $index - index of the '{' character being analyzed
##
##  RETURN VALUE
##    value of the integer after the '{' in the pattern
##
#############################################################################

sub getIntValueAfterBracket
{
    my $pattern = shift;
    my $index = shift;

    my @patArray = split(//, $pattern);
    my $patLength = @patArray;
    my $charStr = ""; # Characters after the '{' being analyzed
    for (my $i = $index; $i < $patLength; $i++)
    {
	my $char = $patArray[$i];
	if ($char =~ m/\d/)
	{
	    $charStr .= $char;
	}
	else
	{
	    last;
	}
    }
    if ($charStr =~ /\d+/) # $charStr is an integer
    {
	return $charStr;
    }
    else
    {
	return 0;
    }
} # getIntValueAfterBracket

# Report error and exit
sub errorReport
{
    my $msg = shift;
    print "$msg\n";
    exit();
}
