#!/usr/bin/perl -w

#############################################################################
#                                                                           #
#   copyright            : (C) 2000-2002 by David J. Studholme              #
#   email                : ds2@sanger.ac.uk                                 #
#                                                                           #
#############################################################################

#############################################################################
#*                                                                          #
#*   This program is free software; you can redistribute it and/or modify   #
#*   it under the terms of the GNU General Public License as published by   #
#*   the Free Software Foundation; either version 2 of the License, or      #
#*   (at your option) any later version.                                    #
#*                                                                          #
#############################################################################
#This is the Perl script used to scan genome sequences against the ModE binding site matrix.#
### VERSION / UPDATE LOG ###

# 7/11/02 Corrected bugs caught by Richard Pau
# 21/11/02 Changed to log base 2 rather than natural log
# 07/03 Enabled reporting of distance upstream 
# enabled normalising of scores.
# 27/08/03 Switched format of matrix files.
# 27/08/03 Integrated Hit object into same file.

use strict ;
use Getopt::Std ;


$| = 1 ;

### Parse command line arguments
### -m specify matrix file, -g specify genome name, -c specify cut-off

my %option ;
getopts ( "hsnim:g:c:" , \%option ) ;

if ( $option{h} ) {
  show_help()  ;
  exit 
}

my ($seqread);
my ($seq_length, @seq_window);
my ($mean_score, $total_dev, $sd, $variance);
my $GC_content; 

### Set threshold score above which counts as a hit
my $cut_off = 80 ;
if ( $option{c} ) {
  $cut_off = $option{c} ;  
}

### Set input/output filenames
my $genome = "ecoli" ;
if ( $option{g} ) {
  $genome = $option{g}
}

### Remove any filename extensions
if ( $genome =~ m/([\w\d\.]+)\.fna/ ) {
  $genome = $1
}

### Name of file containing sequence to be scanned
my $sequence_file = "$genome.fna";  

### Name of file containing details of ORFs
my @ptt_file ;
unless ( $option{n} ) {
	my $ptt_file = "$genome.ptt" ;
	open (PTT, "<$ptt_file") or die "Failed to open ptt file: $genome.ptt\n" ;
	while (my $readline = <PTT> ) {
	push @ptt_file, $readline ;
	}
}

### Name of matrix file
my $matrix_file = $option{m} or die "You must specify a matrix file\n" ;


### Begin the scans!!
promscan_kl($matrix_file, $sequence_file, \@ptt_file)  ;
print "\n\n\n\n" ;
exit ;



sub promscan_kl {
  
  my ($matrix_file, $sequence_file, $ptt_file) = @_;
  
  my %all_scores_forward ;
  my %all_scores_reverse ;
  my @hits ;
  
  
  my ($sequence, $header) ;
  print STDERR "Reading the sequence file into memory ... " ;
  open (SEQFILE, "<$sequence_file") or die "Failed to open $sequence_file\n" ;
  
  while (my $readline = <SEQFILE> ) {
    chomp $readline ;
    $sequence .= $readline if $readline =~ m/^[A-Za-z]+$/ ;
    if ( $readline =~ m/>(.*)/ ) {
      $header = $1 ;
    }
  }
  close SEQFILE ;
  print STDERR "done.\n" ;
  print STDERR "\nSequence: $header\n\n" ;
  print "\nSequence: $header\n\n" ;
  

  print STDERR "Splitting the sequence into a list ... " ;
  my @sequence = split //, $sequence ;
  print STDERR "done.\n" ;
  

  ### Calculate the GC content of the sequence
  my $GC_content = calculate_GC($sequence_file) ;
  print "G+C = $GC_content\n";
  print STDERR "G+C = $GC_content\n";
  my %content ;
  $content{G} = $GC_content / 2 ;
  $content{C} = $GC_content / 2 ;
  $content{A} = (1 - $GC_content) / 2 ;
  $content{T} = (1 - $GC_content) / 2 ;
  
  ### Read matrix file
  my $pmatrix = read_matrix_file($matrix_file) ;
  my %matrix = %$pmatrix ;

  ### Do sanity check on matrix, and calculate width of the matrix
  die "Matrix corrupted\n" unless ( @{$matrix{A}} == @{$matrix{C}} 
				    and @{$matrix{C}} == @{$matrix{G}}
 				    and @{$matrix{A}} == @{$matrix{T}} ) ;
  my $seq_window_length = @{$matrix{A}} ;
  
  
  ### Calculate maximum possible score
  my $overall_highest = 0 ;
  
  for(my $i = 1; $i < $seq_window_length; $i++) {
    my $highest_score = 0 ;
    my $highest_base ;
    foreach my $base ("A", "C", "G", "T" ) {
      if ($matrix{$base}[$i] > $highest_score) {
	$highest_score = $matrix{$base}[$i];
	$highest_base = $base ;
      }
    }
    $overall_highest +=  ( $highest_score * log($highest_score/$content{$highest_base}) / log(2) ) if $highest_score ;
  }
  #print STDERR "(Highest possible score before normalisation = $overall_highest)\n" ;
  #print "(Highest possible score before normalisation = $overall_highest)\n" ;
  

  print STDERR "Sites scoring > $cut_off\n";
  print "Sites scoring > $cut_off\n";
  
  ### Initialise the first window of sequence to be scanned
  for (my $i = 0; $i < $seq_window_length; $i++)  {
    $seq_window[$i] = $sequence[$i] ;
  }
  
  ### read sequence one base at a time
  my $n = 0 ;
  
  while (my $seqread = uc(shift @sequence)) {
    my $seq_window ;
    if ($seqread =~ m/[ACGTN]/) { # sanity-check for legal bases 
      
      ### move the sequence window along one residue
      shift (@seq_window);               # removes first residue
      push (@seq_window,$seqread);       # adds next residue
      
      ### Assign a score to current window
      ### This score, known as the Kullback-Liebler distance, 
      ### reflects the theoretical binding energy of the 
      ### DNA-protein interaction.
      ### For more information see: Stormo, G.D. 2000.
      ### DNA binding sites: representation and discovery. 
      ### Bioinformatics. 16:16-23.
      
      ### First check the forward strand
      my $strand = "+";
     
      ### Calculate score. Remember, divide by log2 to get log to the base of 2
      my $current_score = 0;
      for (my $i = 0; $i < $seq_window_length; $i++) {
	my $base = $seq_window[$i] ;
	$current_score = $current_score + 
	  ( ($matrix{$base}[$i]) * log($matrix{$base}[$i] / $content{$base}) / log(2) )
	    if $matrix{$base}[$i] ;
	}
      
      
      ### Normalise the score
      $current_score /= $overall_highest ;
      $current_score *= 100 ;
      

      if ( $option{s} ) {
	### Record each score in a matrix %all_scores_forward
	$all_scores_forward{$n} = $current_score
      }
      
      ### If this window scores greater than the cut off, then record details ....
      if($current_score > $cut_off) {
	
	### Convert @seq_window into a string, $seq_window
	$seq_window = "" ;
	my $q = 0 ;
	while ($seq_window[$q]) {
	  $seq_window .= $seq_window[$q] ;
	  $q++
	}
	
	### Is hit in coding or non-coding region?
	my ($locate_match, $intergenic_size)  = locate_match($ptt_file, $n, $strand, $sequence, $seq_window);
	
	### Round score to nearest whole number
	$current_score = int ($current_score) ;
	

	
	### Print details of this window
	if ( $locate_match or $option{n} ) {

	  ### Create a new Hit object and set its values
	  my $hit = Hit->new() ;
	  $hit->position($n) ;
	  $hit->strand($strand) ;
	  $hit->score($current_score) ;
	  $hit->strand($strand) ;
	  $hit->sequence($seq_window) ;
	  $hit->locus($locate_match) ;

	  ### Record the existence of this Hit object in an array
	  push @hits, $hit ;

	  printf STDERR "%8s|%1s|%s|%5d|%s\n",
	    $n, $strand, $seq_window, $current_score, $locate_match ;
	}
      }
      


      ### now check the Reverse strand
      $strand = "-";

      ### Calculate score. Remember divide by log2 to get log to the base of 2
      $current_score = 0;
      for (my $i = 0; $i < $seq_window_length; $i++) {
	my $j = $seq_window_length - $i -1 ;
   	my $base = revcomp($seq_window[$j]) ;
	$current_score = $current_score + 
	  ( ($matrix{$base}[$i]) * log($matrix{$base}[$i] / $content{$base}) / log(2) )
	    if $matrix{$base}[$i] ;
      }
      
      ### Normalise the score
      $current_score /= $overall_highest ;
      $current_score *= 100 ;
      
      if ( $option{s} ) {
	### Record each score in a matrix %all_scores_forward
	$all_scores_reverse{$n} = $current_score
      }
            
      
      #### print current score if greater than threshold value
      if($current_score > $cut_off) {
	
	### Convert @seq_window into a string, $seq_window
	$seq_window = "";
	my $q = 0;
	while ($seq_window[$q]) {
	  $seq_window = $seq_window.$seq_window[$q] ;
	  $q++
	}
	
	### Convert $seq_window to Reverse Complement
	my $revcomp =  reverse ($seq_window);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	$seq_window =   $revcomp;
	
	

	### Is hit in coding or non-coding region?
	my ($locate_match, $intergenic_size) = locate_match($ptt_file, $n, $strand, $sequence, $seq_window);
	
	### Round the score to nearest whole number
	$current_score = int ($current_score) ;
	

	
	### Print details of this window
	### ... but only if it is intergenic
	if ( $locate_match or $option{n} ) {
	  
	  ### Create a new Hit object and set its values
	  my $hit = Hit->new() ;
	  $hit->position($n) ;
	  $hit->strand($strand) ;
	  $hit->score($current_score) ;
	  $hit->strand($strand) ;
	  $hit->sequence($seq_window) ;
	  $hit->locus($locate_match) ;
	  
	  ### Record the existence of this Hit object in an array
	  push @hits, $hit ;
	  
	  printf STDERR "%8s|%1s|%s|%5d|%s\n",
	  $n, $strand, $seq_window, $current_score, $locate_match ;
	}
      }
      
      
      ### The sequence pointer proceeds by one base
      $n ++;
      
    }; # end of if($seqread =~ /[ACGTN]/)
    
  }; # end of while (read $sequence_file, $seqread, 1) loop
  
  ### The final value of $n, the sequence-pointer, 
  ### gives the number of readable bases 
  $seq_length = $n ;
  print "N = $seq_length\n";

  ###  Create a hash of arrays, to index Hit objects by score
  ###  Cannot use a simple hash because multiple Hits may have
  ###   the same score.
  my %hitscores ;
  foreach my $hit ( @hits ) {
    my $score = $hit->score() ;
    unless ( $hitscores{$score} ) {
      $hitscores {$score} = [] 
    }
    my $ref = $hitscores{$score} ;
    my @list = @$ref ;
    push @list, $hit ;
    @$ref = @list ;
  }
  
  ### Sort the Hit objects into descending order of scores  
  my @unordered_scores = keys ( %hitscores )  ;
  my @ordered_scores  = sort { $b <=> $a } @unordered_scores ;
  my @ordered_hits ;
  foreach my $score ( @ordered_scores) {
    my $ref = $hitscores{$score} ;
    my @hits = @$ref ;
    foreach my $hit ( @hits ) {
      push @ordered_hits, $hit
    } 
  }
  
  ### Now print out the hits ...
  foreach my $hit ( @ordered_hits ) {
    my ( $pos, $str, $seq, $score, $locus ) = 
      ( $hit->position(), $hit->strand(), $hit->sequence(), 
	$hit->score(), $hit->locus() ) ;
    printf "%8s|%1s|%s|%5d|%s\n",
    $pos, $str, $seq, $score, $locus	  
  }
  

  ### If we have opted to calculate statistics ...
  if ( $option{s} ) {
    
    ### Calculate mean score for all possible windows in the 
    ### query sequence
    my $total_score = 0;
    for(my $i = 0; $i < $seq_length ; $i++) {
      $total_score = $total_score 
	+ $all_scores_forward{$i}
	  +  $all_scores_reverse{$i}
	}
    
    $mean_score = $total_score / ($seq_length * 2 ) ;
    print "Mean score = $mean_score\n";
    
    ### Calculate variance and sd of scores for all 
    ### possible windows in the query sequence
    $total_dev = 0;
    for( my $i = 0; $i < $seq_length; $i++) {
      $total_dev = $total_dev
	+ (($all_scores_forward{$i} - $mean_score)**2)
	  + (($all_scores_reverse{$i} - $mean_score)**2) ;
    }
    
    $variance = $total_dev / ($seq_length * 2) ;
    print "Variance = ", $variance, "\n\n\n";
    $sd = $variance**(0.5)
  }
  
  close SEQFILE;
  
} # end of promscan_kl()


sub locate_match  {
  my $ptt_file = shift or die "Failed to pass a PTT file" ;
  my $promoter_position = shift  or die "Failed to pass a position\n" ;
  my $promoter_strand = shift or die "Failed to pass a strand\n";
  my $sequence = shift or die "Failed to pass a sequence\n" ;
  my $seq_window = shift or die "Filed to pass a seq_window ($ptt_file,$promoter_position,$promoter_strand\n" ;

  unless (defined $option{n} ) {
    
    my ($start_pos1, $end_pos1, $strand1, $product1) = (0, 0, "?", "ORIGIN");
    my ($start_pos2, $end_pos2, $strand2, $product2) = (0, 0, "?", "ORIGIN");
    
    
    my $promoter_description = "";
    ### $promoter_description = "not found";
    foreach my $read_line (@$ptt_file) {  
    
      # Read each window of the protein table file
      my  @fields = read_ptt_line($read_line) ;
    
      unless ($fields[0] eq "failed")    {

	($start_pos1, $end_pos1, $strand1, $product1) = ($start_pos2, $end_pos2, $strand2, $product2) ;
	($start_pos2, $end_pos2, $strand2, $product2) = @fields ;
   
	my $size_of_intergenic_region =  0  ;
	my $upstream_sequence = "" ;
	
	### Check wether promoter falls within this window 
	if ($promoter_position > $end_pos1 && $promoter_position < $start_pos2)  {
	  
	  ###  promoter is in intergenic region...
	  my $size_of_intergenic_region_temp = $start_pos2 - $end_pos1 ;
	  
	  
	  if ($promoter_strand eq "+" && $strand2 eq "+") {
	    $size_of_intergenic_region = $size_of_intergenic_region_temp ;
	    $upstream_sequence = substr($sequence, ($end_pos1-50), ($size_of_intergenic_region+49)) ;
	    
	    my $spacer = $start_pos2 - $promoter_position - length($seq_window) ;
	    $promoter_description = "-N($spacer)- $product2|$size_of_intergenic_region" ; 
	    
	  }

	  if ($promoter_strand eq "-" && $strand1 eq "-") {
	    $size_of_intergenic_region = $size_of_intergenic_region_temp ;
	    my $seq_tmp = substr($sequence, ($end_pos1-2), ($size_of_intergenic_region+52))
	      or die "Failed to get seq_tmp! end_pos1=$end_pos1, size_of_intergenic_region=$size_of_intergenic_region\n" ;
	   	
	    my $spacer = $promoter_position - $end_pos1 ;
	    $promoter_description = "-N($spacer)- $product1|$size_of_intergenic_region" ; 
	    
	  }
	}
      
	elsif  ( ( $promoter_position >= $start_pos1 )
		 and ( $promoter_position <= $end_pos1)
		 and ( $option{i} ) ) {
	  $promoter_description = "internal to $product1"
	}
      }
      
    }
    my $size_of_intergenic_region = $start_pos2 - $end_pos1 ;
    return($promoter_description, $size_of_intergenic_region);
    
  }
}  # End of sub locate_match



sub read_ptt_line {
  
  # Expects a scalar value containing one line of a .ptt file
  # Extracts the 'location', strand, and 'product'
  #   from that line and returns them, in that order
  
  my $read_line = $_[0];
  my $i;
  my $product;
  my @fields;
  
  if ($read_line =~ /(\d+)\.\.(\d+)\s+([+-])/){
    # Extract 'location' and 'strand'    
    @fields = split(/\s+/,$read_line);
    $i = 5;
    $product = "";
    while ($fields[$i]) {
      # Extract 'product'
      $product = $product." ".$fields[$i];
      $i++;
    } 
    return ($1, $2, $3, $product);
  } ; # End of 'if( $read_line =~ ....' block
  
  ### print ("variable read_line:$read_line\n");
  return ("failed");
  
} ; # end of sub read_ptt_line;

sub calculate_GC {
  
  ### Expects the following parameter:
  ### $sequence_file - name of file containing 
  ###   sequence data. 
  ### Expects NCBI .fna format
  ### Returns G+C content in the input sequence
  
  my ($sequence_file) = @_;
  my ($n, $GC, $AT) = (0, 0, 0); # initialise counters
  my $seqread  ;
  
  open(SEQFILE, $sequence_file) ||
    die "Could not open Sequence file: $sequence_file\n" ;
  
  ### ignore the first line of file
  $seqread = <SEQFILE>;   
  
    print STDERR "calculating G+C content\n" ;
  
  while (read SEQFILE, $seqread, 1) {
    # read SEQFILE on residue at a time
    if ($seqread =~ /[GCgc]/) {$GC++};
    if ($seqread =~ /[ATat]/) {$AT++};
  }
  
  close SEQFILE ;
  $n = $GC + $AT ;
  my $GC_content = $GC / $n ;
  print "\nnumber of readable bases = $n\n" ;
  return $GC_content ;
  
}; # end of sub calculate_GC /






sub revcomp {
  my $string = shift @_ or print STDERR "Called revcomp on null string!" ;
  my $seq = "" ;
  my @array = split(//, $string) ;
  while (my $char = pop @array) {
    $char =~ tr/ACGT/TGCA/ ;
    $seq = $seq.$char 
  }
  #print STDERR "Revcomp of $string -> $seq\n" ;
  return $seq 
}




sub read_matrix_file {
  
  my $matrix_file = shift ;
  my %matrix ;

  open(MATRIXFILE, $matrix_file)
    and print STDERR "Opened matrix file: $matrix_file\n"
      or die "Could not open Configuration file: $matrix_file" ;


  print "\nScoring matrix:\n\n" ;

  while (my $readline = <MATRIXFILE>) {
    chomp ($readline) ;
    
    if ( $readline =~ m/([ACGT]):([\s\d]+)/ ) {
      print "$readline\n" ;
      print STDERR "$readline\n" ;  
      
      my $base = $1 ;
      my @values = split /\s+/, $2 ;
      
      $matrix{$base} = [] ;
      foreach my $value (@values) {
	push @{$matrix{$base}} , $value ;
      }
    }
  }
      
  print "\n\n" ;
  print STDERR "\n\n" ;  

  close(MATRIXFILE);
  print STDERR "Closed matrix file\n" ;
  
return \%matrix ;
}





































sub show_help {
    print STDERR <<EOF;


=====================================================================================================
$0: scan a DNA sequence against a motif specified as a scoring matrix
=====================================================================================================

Usage: $0 <options>

    Available options 

        -h            : Show this help.
	-l            : Use linear scoring formula ignoring G+C content,
                        rather than the Kullback-Leibler distance.
        -s            : Calculate statistics (this requires a lot of memory!)
        -i            : Show matches that are internal to ORFs as well as intergenic ones.
        -n            : Do not use a .ptt file (ie no ORF information)
	-g <file>     : Specify the input file (genome)
                        You can specify the .ppt file or the .fna file
                        or just the base name:
                        eg. $0 -g NC_000962.ptt
                            $0 -g NC_000962.fna
                            $0 -g NC_000962
    	-m <file>     : Specify the matrix file.
	-c <number>   : Specify the threshold (cut-off) score (default is 1000).
                        You may wish to experiment with various threshold values
                         to get the optimum balance between sensitivity and accuracy.

======================================================================================================

Output is printed to the output file as well as to STDOUT (ie the screen). 
For each hit, ie a close match to the consensus sequence that is being sought,
the output is arranged into six columns:

Column 1              : The position of the hit in the input sequence.
Column 2              : The strand on which the hit is found (+/-).
Column 3              : Sequence of the hit.
Column 4              : A score for the hit.  The higher the score, the better the match!
Column 5              : The position of the hit relative to adjacent/downstream Open Reading Frames
Column 6              : The size of the intergenic region in which the match is found.

For more information see the PromScan website http://www.promscan.uklinux.net

======================================================================================================
EOF
}










### Hit.pm class of objects representing a hit,
###  i.e. a good match to the sought sequence
###  defined by the frequency matrix

package Hit ;

sub new {
  my $self = {} ;
  bless $self, "Hit" ;
  return $self ;
}


sub position {
  my $self  =  shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{POSITION} = $specified 
  }
  return $$self{POSITION}
}


sub strand {
  my $self = shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{STRAND} = $specified 
  }
  return $$self{STRAND}
}


sub score {
  my $self = shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{SCORE} = $specified 
  }
  return $$self{SCORE}
}



sub sequence {
  my $self = shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{SEQ} = $specified 
  }
  return $$self{SEQ}
}

sub locus {
  my $self = shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{LOCUS} = $specified 
  }
  return $$self{LOCUS}
}

sub intergenic_size {
  my $self = shift @_ ;
  if ( my $specified = shift @_ ) {
    $$self{SIZE} = $specified 
  }
  return $$self{SIZE}
}

return 1 ; # end of package

