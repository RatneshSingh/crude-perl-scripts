#!/usr/bin/perl
use strict;
use warnings;

#print "*****USAGE*****\n\nScript\tseqfile\tspidey_output_file\tgene_outputfile\tsummary_output\n\n";

my (%seq_hash);

#my$seqfile=$ARGV[0];
my$spidey_output_file=$ARGV[0];
#my$gene_output="$seqfile.gene.out";
#my$summary_output="$seqfile.summary.out";
#my$gene_output=$ARGV[2];
#my$summary_output=$ARGV[3];

#open GENE,">$gene_output" or die "gene_output not defined\n";
#open FASTA,"$seqfile" or die "Seqfile not defined\n";
#open ,">$summary_output" or die "summary_output not defined\n";

#####################################################################################################
## Read database sequence in to a hash %seq_hash{}, remove ">" from header.
## remove any white spaces and newline from sequence.
#
#print "reading Sequences from input file.....Plz wait...\n";
#$/="\n>";
#
#while(<FASTA>){#
#    chomp;
#    my($header,@sequence)=split("\n",$_);
#    $header=~s/>//;
#    $header=~s/\s//g;
#
#    my$sequence= join("",@sequence);
#    $sequence=~ s/\s//g;
#    $sequence=~s/\n//g;
#    #print "$header\n";
#    #feed headers and sequences in hash.
#    $seq_hash{$header}=$sequence;
#
#    #print "$seq_hash{$header}\n\n";
#}#
#my @seq_count=keys (%seq_hash);
#my $seq_count=@seq_count;
#
#print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
#@seq_count=();
############################################################################################

#########################################################################################################
# Parsing the spidey output file
#print "reading Spidey out put file.........\n";
spidey_parse($spidey_output_file);

###########################################################################################
# Subroutine for parsing spidey output file for the gene name, co-ordinates of exon and
# introns.Output table with all the co-ordinates arranged in a nice format.

sub spidey_parse{
open SPIDYOUT,(my$spidey_output= shift(@_));
$/= "\nGenomic:";
#print "Spidey output file is:$spidey_output\n";
    print"##gff-version 3\n";
	while(<SPIDYOUT>){
      #print"\nprocessing New Record\n";
		my ($genomic,$mRNA,$strand,$missing_mRNA_ends,%gene,@mRNA_info,$gene_start,$gene_end)=();
		my ($exon_N,$splice_site,$mRNA_coverage,%exon_start,%exon_end,%mRNA_start,%mRNA_end,%exon_percent_identity,%exon_gap,$percent_identity,%exon_mismatch)=0;
		@mRNA_info= split(/\n/,$_);
        $mRNA_info[0]=~s/Genomic://;
        $genomic=$1 if $mRNA_info[0]=~/^\s+([^\s]+)\s+/;
        $mRNA=$1 if $mRNA_info[1]=~/^mRNA:\s+([^\s]+)\s+/;
        foreach my$mRNA_info(@mRNA_info){
            chomp($mRNA_info);
            if($mRNA_info=~/^([\d]+)\s+([\d]+)/){$gene_start=$1;$gene_end=$2;}
            $strand=$gene_start>$gene_end?"-":"+" if ($gene_start && $gene_end);

            if($mRNA_info=~/^Exon\s\d+\s+([\d-]+\s+)/){
              my$exoendstart=$1;
              $exoendstart=~s/\s+//g;
              my($exon_start,$exon_end)=split(/-/,$exoendstart)if defined $exoendstart;
              $exon_start{$exon_start}=$exon_start;
              $exon_end{$exon_end}=$exon_end;
              $strand=$exon_start>$exon_end?"-":"+";
            }
        next;
        }# foreach loop.


	my@a= sort {$a <=> $b} values(%exon_start) if $strand eq "+";
	my@b= sort {$a <=> $b} values(%exon_end) if $strand eq "+";
    @a= sort {$b <=> $a} values(%exon_start) if $strand eq "-";
	@b= sort {$b <=> $a} values(%exon_end) if $strand eq "-";
	my@c= sort(@a,@b);
	my$l=@c;
	#print "\n$a[$l-1]\n";
	#print "@a";



	#print "\ngene_start=$exon_start{1}\tgene_end=$exon_end{$exon_N}";
#
	#($gene{$mRNA})= extract_gene($c[0],$c[$l-1],$genomic,$strand);
	#print ">Gene_$mRNA\n$gene{$mRNA}\n" if defined $gene{$mRNA} ;
    foreach($genomic,$gene_start,$gene_end,$strand,$mRNA){s/s+//g;}
    print "$genomic\tspidey\tgene\t$gene_start\t$gene_end\t.\t$strand\t.\tID=$mRNA\n";
    for(my$i=0;$i<scalar@a;$i++){

      print"$genomic\tspidey\tCDS\t$a[$i]\t$b[$i]\t.\t$strand\t.\tID=$mRNA.CDS",$i+1,";Parent=$mRNA\n";

    }


	}# while loop.
    print"\n";


#close(GENE);
#close();
#close(FASTA);
#close(GENE);

# information collected: Genomic DNA=>$genomic,mRNA Name=>$mRNA,Strand=>$strand,Number of exons=>$exon_N,
# splice sites=>$splice_site,mRNA coverage=>$mRNA_coverage,Missing mRNA ends: $missing_mRNA_ends
# Exon start=>$exon_start{$i},Exon_end=>$exon_end{$i},mRNA start=>$mRNA_start{$i},
# mRNA end=>$mRNA_end{Si},$exon_percent_identity{$i},$exon_gap{$i}




}# End of subroutine spidey_parse
