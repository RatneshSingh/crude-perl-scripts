#!/usr/bin/perl
#### added sorting of coordinates before extracting mRNA sequence to avoid problems due to misarranged gff files.
####
####
####
####
####
####
####
####
####
####
####
####
####
####
####







use strict;
use Text::Wrap;
$Text::Wrap::columns = 59;
$Text::Wrap::separator2="\n                    ";
$Text::Wrap::break=",";

my %Gene;
my %gene_child;
my $gene_id;
my $mRNA_id;
my %trace_parent;
my %colour_list;
my$opt_delimiter=" ";
my$opt_column=1;
my$add_colour='true';


####################################
die "\n\n$0\tgff_file\tfasta_file\tcolour_file[optional]\n\n" if !$ARGV[0];


####################################
### collect colour information from table of two colums: gene_name   colour_info
if($ARGV[2]){

     open TABLE,"$ARGV[2]";
     while(<TABLE>){
          s/^\s+//g;
          my($gene,@colour)=split /\s+/;
          $gene=~s/\s+//g;
          my$colour=join(" ",@colour);
          $colour=~s/^\s+|\s+$//g;
          $colour_list{$gene}=$colour;

     }

}

open GFF, "$ARGV[0]";



my ($genome_seq,$seq_len) = ReadFasta($ARGV[1]);
my$chromosome;
##### collect parent child info first.
while ( my $line = <GFF> ) {
     my ( $seqid, $source, $type, $start, $end, $score, $strand, $phase, @attributes ) = split /\s+/, $line;
    next if $type!~/mRNA/i;
    chomp( $seqid, $source, $type, $start, $end, $score, $strand, $phase, @attributes );
    s/^\s*//g foreach ( $seqid, $source, $type, $start, $end, $score, $strand, $phase );
    s/\s*$//g foreach ( $seqid, $source, $type, $start, $end, $score, $strand, $phase );
    $type =~ s/\s+//g;
    my $attribute = join( "", @attributes );
    my $id     = $1 if $attribute =~ /ID=([^;]+)/i;
    my $parent = $1 if $attribute =~ /Parent=([^;]+)/i;
    next if $id eq "" && $parent eq "";
    if ( $type =~ /mRNA/i ) { $trace_parent{$id} = $parent }
    $chromosome=$seqid;

}
close GFF;

######collect gene names in GFF file
open GFF, "$ARGV[0]";
while ( my $line = <GFF> ) {

    my ( $seqid, $source, $type, $start, $end, $score, $strand, $phase, @attributes ) = split /\s+/, $line;
    chomp( $seqid, $source, $type, $start, $end, $score, $strand, $phase, @attributes );
    s/^\s*//g foreach ( $seqid, $source, $type, $start, $end, $score, $strand, $phase );
    s/\s*$//g foreach ( $seqid, $source, $type, $start, $end, $score, $strand, $phase );
    $type =~ s/\s+//g;
    my $attribute = join( "", @attributes );
    my $id     = $1 if $attribute =~ /ID=([^;]+)/i;
    my $parent = $1 if $attribute =~ /Parent=([^;]+)/i;
    next if $id eq "" && $parent eq "";

    if ( $type =~ /gene/i ) { push (@{$Gene{$id}{'gene'}{$id}},$start, $end );$Gene{$id}{'strand'}{'strand'}=$strand }

    if ( $type =~ /mRNA/i ) { push (@{$Gene{$parent}{'mRNA'}{$id}},$start, $end ); $trace_parent{$id} = $parent }
    if ( $type =~ /CDS/i )  {
      push (@{$Gene{ $trace_parent{$parent} }{'CDS'}{$parent}},[$start, $end]) if $strand eq "+";
      unshift (@{$Gene{ $trace_parent{$parent} }{'CDS'}{$parent}},[$start, $end]) if $strand eq "-";
      }
}

###sort array containing CDS coordinates to avoid problems due to unarranged gff files
foreach my $gene_name ( keys %Gene ) {
     foreach my$mRNA(keys %{$Gene{$gene_name}{'CDS'}}){
          for(my$coord=0;$coord <= scalar@{$Gene{$gene_name}{'CDS'}{$mRNA}}-1;$coord++){
               @{${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}= sort {$a<=>$b} @{${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]};
          }
     }
}


write_gb_header($chromosome,$seq_len);


foreach my $gene_name ( keys %Gene ) {
    next if $gene_name =~/^\s*$/;
    printf("\n%4s %-15s %-10s","",'gene',"${$Gene{$gene_name}{'gene'}{$gene_name}}[0]..${$Gene{$gene_name}{'gene'}{$gene_name}}[-1]") if (${$Gene{$gene_name}{'gene'}{$gene_name}}[0] && ${$Gene{$gene_name}{'gene'}{$gene_name}}[-1] && $Gene{$gene_name}{'strand'}{'strand'} eq "+");

    printf("\n%4s %-15s %-10s","",'gene',"complement\(${$Gene{$gene_name}{'gene'}{$gene_name}}[0]..${$Gene{$gene_name}{'gene'}{$gene_name}}[-1]\)") if (${$Gene{$gene_name}{'gene'}{$gene_name}}[0] && ${$Gene{$gene_name}{'gene'}{$gene_name}}[-1] && $Gene{$gene_name}{'strand'}{'strand'} eq "-");

    printf ("\n%4s %-15s %-10s ",'','',"\/gene=\"$gene_name\"")  if ${$Gene{$gene_name}{'gene'}{$gene_name}}[0] && ${$Gene{$gene_name}{'gene'}{$gene_name}}[-1];
    printf ("\n%4s %-15s %-10s ",'','',"\/colour=$colour_list{$gene_name}")  if $colour_list{$gene_name};
    foreach my$mRNA(keys %{$Gene{$gene_name}{'CDS'}}){
      my $cds;

      for(my$coord=0;$coord <= scalar@{$Gene{$gene_name}{'CDS'}{$mRNA}}-1;$coord++){
            $cds.= "${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[0]..${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[1]";
            $cds.= "," if (scalar@{$Gene{$gene_name}{'CDS'}{$mRNA}} > $coord+1);
          }

     if($cds=~/,/){
          my$ocds='join('.$cds.')' if $Gene{$gene_name}{'strand'}{'strand'} eq "+";
          $ocds = 'complement(join('.$cds.'))' if $Gene{$gene_name}{'strand'}{'strand'} eq "-";
          $cds=$ocds;
     }
     else{
          my$ocds=$cds;
          $ocds = 'complement('.$cds.')' if $Gene{$gene_name}{'strand'}{'strand'} eq "-";
          $cds=$ocds;

     }
      printf("\n%4s %-15s %-10s","",'CDS',"$cds") if $cds!~/^\s*$/;
      printf ("\n%4s %-15s %-10s ",'','',"\/colour=$colour_list{$gene_name}")  if $colour_list{$gene_name};




      }

  }

print"\nORIGIN\n";
foreach my$chr(keys %$genome_seq){
  print_gb_seq($$genome_seq{$chr});
}


#write_gb_header($genome_name,$seq_len);
print "//";



######print CDS from coordinates and sequence infor each gene in gff file.

my@cds_file=split(/\//,$ARGV[1]);
$cds_file[-1]=~s/\s+//g;

open my$fh,">$cds_file[-1].cds.fasta";
extract_mRNA($$genome_seq{$chromosome},$fh);



########################################################################################
### subroutines
########################################################################################

sub ReadFasta{ # to read fasta format files into hash. returns hash.

	my $seqfile=shift(@_);
	my$demo_header;
	my$seq_len=0;
	my ($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile";
	#print "Reading Sequences from input file.....Plz wait...\n";
	my%seq_hash=();
	#$seq_hash{'RS_Concatenated'}="";

	$/="\n>";    # Change record seperator to read Fasta
	my$last_N=1;
	while(<FASTA>){
    	chomp;
    	($header,@sequence)=split("\n",$_);

    	$header=~s/>//;						# Remove Leading > from Header
    	$header=~s/\s*$//;					# Remove trailing spaces from header
    	$header=~s/^\s*//;					# Remove Leading spaces from Header
        my@scaffold_name=split(/$opt_delimiter/,$header);
        $header=$scaffold_name[$opt_column - 1];
    	my$sequence= join("",@sequence);
    	$sequence=~ s/\s//g;
    	$sequence=~s/\n//g;
		$seq_len+=length($sequence);

    	if($header=~/^\s*$/){next;}
    	$demo_header=$header;
    	# Make sure no duplicate header exists else they will be overwritten and lost during hashing due to same key.
    	if(!exists $seq_hash{$header}){
    		$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
			#$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
		}
		else {

			# find a uniq header name by adding a number at the end. If header still exists, increase the number by one
			while(exists $seq_hash{$header}){$header=$header.$last_N;$last_N++;}

			$seq_hash{$header}=$sequence;
			#$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;


		}
	}

	my @seq_count=keys (%seq_hash);
	my $seq_count=@seq_count;

	#print "Done....\nNumber of sequences read form input file = $seq_count\n\nExample of one the headers in sequence\n$demo_header\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.

	return(\%seq_hash,$seq_len);

}



sub write_gb_header{
  my$genome_name=shift;
  my$seq_len=shift;

  #### get date


printf("%-16s %-10s\n","LOCUS","$genome_name     $seq_len bp DNA  PLN  13-Nov-2013");
printf("%-16s %-10s\n","ACCESSION", $genome_name);
printf("%-16s %-10s\n","VERSION", 1);
printf("%-16s %-10s\n","KEYWORDS",$genome_name);
printf("%0s %16s %10s",'FEATURES'," ",'Location/Qualifiers');
printf("\n%4s %-15s %-10s"," ",'source',"1..$seq_len");



}


sub write_genbank{

  ###########################

  foreach my $gene_name ( keys %Gene ) {
    next if $gene_name =~/^\s*$/;
    print "\n     gene           ${$Gene{$gene_name}{'gene'}{$gene_name}}[0]..${$Gene{$gene_name}{'gene'}{$gene_name}}[-1]";
    print "                     \/gene=\"$gene_name\"";

    foreach my$mRNA(keys %{$Gene{$gene_name}{'CDS'}}){
      print "\n     mRNA           join(";
      foreach my$coord(@{$Gene{$gene_name}{'CDS'}{$mRNA}}){
        foreach my$ind(0..scalar@$coord-1){
          print "$$coord[$ind].." if $ind%2>0;
          print "$$coord[$ind]," if $ind%2>0;

        }
        print ")";
      }

  }
}
}

sub print_gb_seq{
  my$sequence=shift;
  $sequence=~s/[\W]+//g;
  $sequence=~s/[\s]+//g;
  my@sequence= split(/(.{60})/,$sequence);

  @sequence=map{join(" ", split(/(.{10})/,$_))} @sequence;
  my$count=1;
  #@sequence=map{$_.=$count;$count+=60} @sequence;

  for(my$i=0;$i<=scalar@sequence-1;$i++){
    next if $sequence[$i]=~/^\s*$/;
    $sequence[$i]=sprintf("%10s %2s",$count,"$sequence[$i]\n");
    $count+=60
  }

  $sequence= join("\n", @sequence);
  $sequence=~s/\n+/\n/g;
  print $sequence;

}

sub extract_mRNA{
     my$sequence=shift;
     my$fh=shift;

foreach my $gene_name ( keys %Gene ) {
    next if $gene_name =~/^\s*$/;
    my$complement='false';
    $complement='true' if (${$Gene{$gene_name}{'gene'}{$gene_name}}[0] && ${$Gene{$gene_name}{'gene'}{$gene_name}}[-1] && $Gene{$gene_name}{'strand'}{'strand'} eq "-");
    my$mRNA_seq;
  #  foreach my$mRNA(keys %{$Gene{$gene_name}{'CDS'}}){
  #    for(my$coord=0;$coord <= scalar@{$Gene{$gene_name}{'CDS'}{$mRNA}}-1;$coord++){
  #        $mRNA_seq.=substr($sequence,${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[0]-1, ${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[1]-${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[0]+1);
  #   }
  #   revcomp(\$mRNA_seq) if $Gene{$gene_name}{'strand'}{'strand'} eq "-";
  #   #print "\n******** Reverse complimenting*******\n" if $Gene{$gene_name}{'strand'}{'strand'} eq "-";
  #   print $fh ">$gene_name\n$mRNA_seq\n";
  #}


     foreach my$mRNA(keys %{$Gene{$gene_name}{'CDS'}}){
          my@starts;
          my@ends;
      for(my$coord=0;$coord <= scalar@{$Gene{$gene_name}{'CDS'}{$mRNA}}-1;$coord++){
          push(@starts,${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[0]<${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[1]?${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[0]:${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[1]);
          push(@ends,${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[0]>${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[1]?${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[0]:${${$Gene{$gene_name}{'CDS'}{$mRNA}}[$coord]}[1]);
     }
      @starts=sort {$a<=>$b} @starts;
      @ends=sort {$a<=>$b} @ends;
      for(my$i=0;$i<scalar@starts;$i++){
          $mRNA_seq.=substr($sequence,$starts[$i]-1, $ends[$i]-$starts[$i]+1);
      }
     revcomp(\$mRNA_seq) if $Gene{$gene_name}{'strand'}{'strand'} eq "-";
     #print "\n******** Reverse complimenting*******\n" if $Gene{$gene_name}{'strand'}{'strand'} eq "-";
     print $fh ">$gene_name\n$mRNA_seq\n";
  }





}
}

sub revcomp{
    my$seq=shift;
    $$seq = reverse $$seq;
    $$seq =~ tr/ATGCNatgcn/TACGNtacgn/;
    return;
    #return($$seq);
}