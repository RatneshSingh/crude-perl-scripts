
#!/usr/bin/perl -w
use strict;
my (@files,$seqName);
#**************************************************************
#This script is exclusively for parsing output text files from "PCR NOW"
#program(web interface). All the text files should be in one folder. this script will read all the text files and parse through them one by one 
#reading information about primers. it will print them in a text file names "PCR_NOW_Primers.txt". this file can be directly be viewed in excel
#program. 
#**************************************************************


print " \n\n\n*****Version 1.2 This script is exclusively for PCR NOW out put. All the Primer information will be written in 'PCR_NOW_Primer.txt' file.******\n\n"; 



#open the directory containing files(current directory".").
opendir(DIR, ".");
#read filesnames (readdir(DIR)  having .txt pattern in an array
my @files = grep(/output.*\.txt$/,readdir(DIR));
closedir(DIR);



#open out put file to keep primer information.
open OUT,">PCR_NOW_Primers.txt";
print OUT"SeqName\tPrimerName\tPrimerSeq\tstart\tlen\ttm\tgc\tProd_Size";
open OUT2,">PCR_NOW_Primer.fasta";

#processing files and printing information.
foreach my $file (@files) {
print "Processing file $file. Plz Wait......\n";

open FILE,"$file" or die "Can't open file $file";

$/="\n";
#reading each file one a a time and parsing through.
while (<FILE>) {
        #print "This is record $_ \n";
                                                              
  
  
#looking for sequence name.   
    if (/^PRIMER\s+PICKING\s+RESULTS\s+FOR\s+(\w+)/) {
                  $seqName=$1;
              #~ print "Sequence Name is $seqName\n";
                  print OUT"\n$seqName";                                                
    }                                                      
  #looking for line containing information about left primer. Only first line will be taken.                                                      
    elsif(/^LEFT\s+PRIMER\s+.*/) {
                  my$left_line=$_;
                    $left_line=~s/\s+/,/g;
              #~ print "after replacing spaces with tabs $left_line\n";
                    
                  my($a,$b,$start,$len,$Tm,$GC,$any,$threeprime,$rep,$L_primer)=split(/,/,$left_line);

                  print "L_primer start = $start\n";                                               
                  print "L_primer length = $len\n";                                                
                  print "L_primer Tm = $Tm\n";                                                
                  print "L_primer GC = $GC\n";                                                
                  print "L_primer any = $any\n";                                                
                  print "L_primer length = $threeprime\n";                                                
                  print "L_primer length = $rep\n";                                                
                  print "L_primer sequence = $L_primer \n";                                                
                  print OUT"\n\t$seqName",'_L',"\t$L_primer\t$start\t$len\t$Tm\t$GC\t";
                  print OUT2">$seqName",'_L',"\n$L_primer\n";
    }        
                
#looking for right primer information containing line.        
    elsif(/^RIGHT\s+PRIMER\s+.*/) {
                  my  $right_line=$_;
                    $right_line=~s/\s+/,/g;
#                    print "after replacing spaces with tabs $left_line\n";
                  my ($a,$b,$start,$len,$Tm,$GC,$any,$threeprime,$rep,$R_primer)=split(/,/,$right_line);

                  print "R_primer start = $start\n";                                           
                  print "_primer length = $len\n";                                                
                  print "R_primer Tm = $Tm\n";                                                
                  print "R_primer GC = $GC\n";                                                
                  print "R_primer any = $any\n";                                                
                  print "R_primer length = $threeprime\n";                                                
                  print "R_primer length = $rep\n";                                                
                  print "R_primer sequence = $R_primer \n";                                                
                  print OUT"\n\t$seqName",'_R',"\t$R_primer\t$start\t$len\t$Tm\t$GC\t";
                 print OUT2">$seqName",'_R',"\n$R_primer\n";
    }
#looking for the product size.
    if(/^PRODUCT\s+SIZE:\s+(\d+).*/){
                  my $prod_size=$1;
                #~ print "Product Size is=$prod_size\n";
                  print OUT"$prod_size\t";}            
          }
            }
#**************************************************************

print "Done......................\n";
exit;
