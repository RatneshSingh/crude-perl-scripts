# This script will convert your DNA sequence to PROTEIN Sequence in 3 frames. Stop codons will be replaced by _ and ambigous codons with 'X'.

# While executing this script it asks for the file name of the DNA sequence. If the sequence file is not available in the same directory of this script, enter the name of the file along with the path.  eg.In windows:  c:\dnafile.txt, In Linux: /home/user/sequence/dnafile.txt

print "\n\n\t\#################### DNA 2 PROTEIN #################### \n\n";
print "This script will convert your DNA sequence to PROTEIN Sequence\n\n";
#print "ENTER THE FILENAME OF THE DNA SEQUENCE:= ";

$DNAfilename = $ARGV[0];
chomp $DNAfilename;
unless ( open(DNAFILE, $DNAfilename) ) {
    print "Cannot open file \"$DNAfilename\"\n\n";
}

open(OUT,">$DNAfilename.translated.fasta");


print "reading Sequences from input file.....Plz wait...\n";
$/="\n>";

while(<DNAFILE>){#
    chomp;
    ($header,@sequence)=split("\n",$_);
    $header=~s/>//;
    $DNA= join("",@sequence);
    $DNA=~ s/\s//g;
    $DNA=~s/\n//g;
    #print "$header\n";
    #feed headers and sequences in hash.
    $seq_hash{$header}=$DNA;

    #print "$seq_hash{$header}\n\n";

my $protein_frame1='';
my $protein_frame2='';
my $protein_frame3='';
my $codon;
#for Frame 1
for(my $i=0;$i<(length($DNA)-2);$i+=3)
{
$codon=substr($DNA,$i,3);
$protein_frame1.=&codon2aa($codon);
}
print OUT">$header"."_Prot-Frame1"."\n$protein_frame1\n";


#for Frame 2
for(my $i=1;$i<(length($DNA)-2);$i+=3)
{
$codon=substr($DNA,$i,3);
$protein_frame2.=&codon2aa($codon);
}
print OUT">$header"."_Prot-Frame2"."\n$protein_frame2\n";

#for Frame 2
for(my $i=2;$i<(length($DNA)-2);$i+=3)
{
$codon=substr($DNA,$i,3);
$protein_frame3.=&codon2aa($codon);
}
print OUT">$header"."_Prot-Frame3"."\n$protein_frame3\n";


}#


#####################################################################
sub codon2aa{
my($codon)=@_;
$codon=uc $codon;
my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
if(exists $g{$codon})
{
return $g{$codon};
}
else
{
print STDERR "Bad codon \"$codon\"!!\n";
return 'X';
}
}