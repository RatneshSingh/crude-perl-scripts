#!/usr/bin/perl -w
use strict;
#use Bio::Tools::GFF;
#use Getopt::Std;
use Getopt::Long;

##################################################################
# usage description
my$usage="

This script reads the features from GFF3 file and returns the requested information for the genes/repeats listed in list file or provided manually.

usage: perl scriptName -options
Options:
-seq|s		Sequence file containing genomic sequences.
-gff|g     		 GFF file name containing coordinates for genes/mRNA/CDS/exon/Repeats etc.
-list|l		'manual'|List file name. List of gene names to extract information for.
		 Skip this option to use all the sequence header as a list.
-feature|f		type of function to perform for the gene names in the list.[CDS]
		'mask': Mask the region containing the genes in the list. print the masked genome in new file.
		'gene': print the gene sequences in new file.
		'mRNA': extract the mRNA sequences based on the coords defined in gff file.
		'exon':extract and join the exon sequences based on the coords defined in gff file.
		'CDS': extract the CDS sequences based on the coords defined in gff file.
		'genbank':Save sequences and features as genbank file for each gene in the list [Not implemented Yet].
		'promoter': Extract promoter region of the genes in the list. The length will be specified by -us (upstream) or -ds (downstream) or both.
		'extract_repeats': Extract repeat regions from genomic sequences using gff file from RepeatMasker.
		'mask_repeats': Mask repeat region on the genomic sequences using gff file from RepeatMasker.
                'gc_around_repeats': Calcuate GC content around repeat regions on genomic sequences using gff file from RepeatMasker.
		'intron-exon-gc':Calculate several statistics about GC content for exons-introns.following types of results will be saved in
                        seperate files as well as one file containing all stat for further analysis.
                        GC_each_introns, GC_each_exon, Len_each_intron, Len_each_exon, GC per 10% of length, GC per sliding window (sliding
                        window is in %of length and not nt).
-window  percent of CDS length to be used as sliding window for intron-exon-gc [10]
-slide   percent of CDS length to be used as step size for sliding window in intron-exon-gc [5]
-out|o		prefix for the output file name. Please dont add extension.
-minlen minimum length of the feature to be used. Discard features smaller than this value[100]
Incase the gene name is short and genomic name is long. Define here
the delimiter to split and which column to use for matching

-delim|d		'delimiter' to split the genomic headers [' ']
-col|c  		column number to use for match [1]
-upsteam|us  		Length of upstream region to be extracted for genes [3000]. Valid when used with -f promoter or -f extract_repeats
-downstream|ds  		Length of downstream region to be extracted for each gene [1000]. Valid when used with -f promoter or -f extract_repeats
-pattern|p		Name feature name to look for if not features names in gff3 file are not commonly used names.
-min_len_gc|n		Minimum length required for upstream and downstream seq for GC calculations. Work with -f gc_around_repeats[100].
-min_len_rep|e		Minimum length of Repeat element to include in GC content[10].
-run_mode|mode|m	'speed'|'memory'.Speed while parsing GFF file [memory].
			'memory': for memory efficient run. Only picks the GFF lines for gene names in pattern file (Slow but memory efficient).
			'speed': for speed. slurps all the lines in GFF file (fast but consumes memory).
-rep_group	table of repeat name(col 1) and group(col2) they belong, to classify names into order/family.
-h			print help
-v          verbose.
";


##################################################################
# option processing
our($opt_seqfile,$opt_GFF_file,$opt_gene_list_file,$opt_feature_type,$opt_out_file,$opt_delimiter,$opt_column,$opt_upstream,
	$opt_downstream,$opt_run_modein_len_for_GC,$opt_feature_name,$opt_run_mode,$opt_min_rep_len,$opt_repeat_group_list,$help,
	$print_raw,$slide,$window,$verbose,$len_lim); ### percent of length;);
my(@pattern,%Gene);
#$opt_GFF_file='Lotus_Genes_short_names.gff' ;
#$opt_gene_list_file='Nn_MADS_39gene.list';
#$opt_seqfile='16738-chromosome.fasta';
$opt_feature_type='CDS';
$opt_delimiter=' ';
$opt_column=1;
$opt_upstream=1000;
$opt_downstream=1000;
$opt_run_modein_len_for_GC=100;
$opt_run_mode='memory';
$opt_min_rep_len=10;
$slide=5; ### percent of length for step size in intron-exon-gc sliding window
$window=10; ### percent of length for window size in intron-exon-gc sliding window

my$result=GetOptions(
						"sequence|seq|s=s" => \$opt_seqfile,
						"gff|g=s" => \$opt_GFF_file,
						"list|l=s" => \$opt_gene_list_file,
						"feature|f=s" => \$opt_feature_type,
						"out|o=s" => \$opt_out_file,
						"slide|bin=f"=>\$slide,
						"window|w=f"=>\$window,
						"delim|delimiter|d=s" => \$opt_delimiter,
						"column|col|c=i" => \$opt_column,
						"upstream|up|us|u=i" => \$opt_upstream,
						"downstream|down|ds|b=i" => \$opt_downstream,
						"min_len_gc|mlg|n=i" => \$opt_run_modein_len_for_GC,
						"pattern|type|p=s" => \$opt_feature_name,
						"run_mode|mode|m=s" => \$opt_run_mode,
						"min_len_rep|mlr|e=i" => \$opt_min_rep_len,
						"rep_group=s" => \$opt_repeat_group_list,
						"print_raw" => \$print_raw,
						"help|h" => \$help,
                        "verbose|v"=>\$verbose,
                        "len_lim|minlen=i"=>\$len_lim,
);


die "\n$usage\n" if $help;


##################################################################
# correct for possible mispelling in function type.
$opt_feature_type = 'exon'   if lc$opt_feature_type eq 'exons';
$opt_feature_type = 'genbank' if lc$opt_feature_type eq 'genebank';
$opt_feature_type = 'transcript' if lc$opt_feature_type eq 'transcripts';
##################################################################
# check for genomic sequence file
die "$usage\nPlease provide file containing genomic sequence with -s flag.\n" if !$opt_seqfile;
# check for gff file
die "$usage\nPlease provide gff file with -g flag\n" if !$opt_GFF_file;

##################################################################
# Read genome sequence information in a hash and return reference of the hash containing sequence and total length of sequences read.
my ($genome_seq,$seq_len) = ReadFasta($opt_seqfile);

##################################################################
# Read gene names from list and create an array @pattern or ask user for manually enter the list.
if($opt_gene_list_file && lc $opt_gene_list_file ne 'manual'){
    open LIST,$opt_gene_list_file or die "Cannot find list file $opt_gene_list_file";
    while(<LIST>){
        next if /^\s*$/;
        #print "Read pattern $_\n";
        chomp $_;
        push(@pattern,$_);
    }
}
# create list from manually entered names
elsif($opt_gene_list_file && lc $opt_gene_list_file eq 'manual'){
    print "\n\nType the name of the gene to retrieve information for\nMultiple gene names should be seperated by spaces:";
    my$gene_names=<STDIN>;
    @pattern=split(/\s+/,$gene_names);
    foreach my$pat(@pattern){$pat=~s/\s+//g;}
}

# Create list from the name of the sequence file
elsif ( !$opt_gene_list_file ) {
	print "\nList file is not provided. Printing all the genes in gff file.\n";
    #@pattern=keys %{$genome_seq};

}
##################################################################
#  get_pattern_from_gff: gets you the reference to hash containing lines matching pattern of gene name.
#  get_gene_info: gets you the reference to hash containing information of asked gene.
my$num_pattern=@pattern;

#$num_pattern > 0 or die "\nNo element in the list of pattern to look for. Exiting.\n";

my$num_done=0;

if($opt_run_mode eq 'memory'){
	foreach my$gene_name(@pattern){
		$num_done++;

		get_gene_info($opt_GFF_file,$gene_name,$opt_run_mode); # gene info will be stored into global variable %Gene
		 print "Read co-ordinates for $gene_name from gff file($num_done out of $num_pattern)\n";

	}
}
else{
    print "\nReading all the elements in GFF file in memory....Please wait";
    get_gene_info($opt_GFF_file,'no need of any name',$opt_run_mode);
    @pattern = keys %Gene if ( !$opt_gene_list_file );
    print ".........Done\n"
}



##################################################################
# Prepare files to save output
my$fh;


if($opt_out_file && (lc$opt_feature_type eq 'gc_around_repeats' || lc$opt_feature_type eq 'intron-exon-gc')){open $fh,">$opt_out_file.$opt_feature_type.table" or die "Cannot open output file for writing\n"; }
elsif($opt_out_file){open $fh,">$opt_out_file.$opt_feature_type.fasta" or die "Cannot open output file for writing\n"; }
else{$fh=*STDOUT;}

## print headers for some output types
my ($fh_gcperintron,$fh_gcperexon,$fh_lenperintron,$fh_lenperexon,$fh_perlen,$fh_perlensliding);
if (lc$opt_feature_type eq 'intron-exon-gc'){
	if (!$opt_out_file) {$opt_out_file=$opt_seqfile;$opt_out_file=~s/\.[^\.]+//;}
	open $fh_gcperintron,">$opt_out_file.$opt_feature_type.GCperintron.table";
	open $fh_gcperexon,">$opt_out_file.$opt_feature_type.GCperexon.table";
	open $fh_lenperintron,">$opt_out_file.$opt_feature_type.Lenperintron.table";
	open $fh_lenperexon,">$opt_out_file.$opt_feature_type.Lenperexon.table";
	open $fh_perlen,">$opt_out_file.$opt_feature_type.perlen.table";
	open $fh_perlensliding, ">$opt_out_file.$opt_feature_type.perlensliding.table";

	print $fh join "\t","GeneName","NumIntrons","gene_length","CDS_length","Intron_length","\%GC_exons","\%GC_introns","GC_each_exons","GC_each_introns","len_each_exon","len_each_intron","GC_per_10percent_of_CDS";
	print $fh_gcperintron join "\t","GeneName","GC_each_introns";
	print $fh_gcperexon join "\t","GeneName","GC_each_exons";
print $fh_lenperintron join "\t","GeneName","Len_each_introns";
	print $fh_lenperexon join "\t","GeneName","Len_each_exons";
	print $fh_perlen join "\t","GeneName","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%";
	my@iegc_header;
	my$count=0;
	my$frag=int((100-$window)/$slide);
	for(my$i=0;  $i <= 100-$window;  $i=int($i+$slide)){my$next_i=int($i+$window)<100?int($i+$window):100; push(@iegc_header,"$i-$next_i\%")}
	print $fh_perlensliding join "\t","GeneName",@iegc_header;
}
##################################################################
# set default values for pattern based on feature(-f) type, if user did not provide the value for pattern (-p).
if(lc$opt_feature_type eq 'mask' && !$opt_feature_name){$opt_feature_name='mask'}
elsif(lc$opt_feature_type eq 'gene' && !$opt_feature_name){$opt_feature_name='gene'}
elsif(lc$opt_feature_type eq 'transcript' && !$opt_feature_name){$opt_feature_name='transcript'}
elsif(lc$opt_feature_type eq 'mrna' && !$opt_feature_name){$opt_feature_name='mRNA'}
elsif(lc$opt_feature_type eq 'cds' && !$opt_feature_name){$opt_feature_name='CDS'}
elsif(lc$opt_feature_type eq 'exon' && !$opt_feature_name){$opt_feature_name='exon'}
elsif(lc$opt_feature_type eq 'genbank' && !$opt_feature_name){$opt_feature_name='genbank'}
elsif(lc$opt_feature_type eq 'promoter' && !$opt_feature_name){$opt_feature_name='gene'}
elsif(lc$opt_feature_type eq 'gc_around_repeats' && !$opt_feature_name){$opt_feature_name='dispersed_repeat'}
else{$opt_feature_name='CDS'}

##################################################################
# Read gff file and retrieve the feature values
print "\nPerforming requested ($opt_feature_type) action on sequences.....\n";
my@attributes_gcup_gcdown;
foreach my$gene_name(@pattern){
    my ($gene_seq,$coords,$seqlen,$seqid,$array_introns,$array_exons);
	# to Mask the feature
    if(lc $opt_feature_type eq 'mask'){
		($gene_seq,$coords,$seqlen,$seqid)=get_seq($gene_name,'mask',$opt_feature_name);
		print {$fh}">$gene_name Exon masked:$seqid$coords length:$seqlen \n$gene_seq\n" if $gene_seq ne "0";
	}
	# toextract the feature
    elsif(lc$opt_feature_type eq 'gene'||lc$opt_feature_type eq 'mrna'||lc$opt_feature_type eq 'cds'||lc$opt_feature_type eq 'exon'||lc$opt_feature_type eq 'transcript'){
		($gene_seq,$coords,$seqlen,$seqid)=get_seq($gene_name,'extract',$opt_feature_name);
		print {$fh}">$gene_name $opt_feature_name:$seqid$coords length:$seqlen \n$gene_seq\n" if $gene_seq ne "0";

        #($gene_seq,$coords,$seqlen,$seqid)=get_introns($gene_name,'extract',$opt_feature_name);
        #print {$fh}">$gene_name Intron:$seqid$coords length:$seqlen \n$gene_seq\n" if $gene_seq ne "0";

	}
    # to extract the intron and exons information
    elsif(lc$opt_feature_type eq 'intron-exon-gc'){


        my($CDS_Seq,$coords,$CDS_len,$seqid,$array_exons)=get_seq($gene_name,'extract',$opt_feature_name);
        my$num_exons=$array_exons?scalar@$array_exons:0 if $array_exons;
        next if length($CDS_Seq) < 100;

		###return ($gene_seq,"NA",0,$seq_id,[])
        my($intron_seq,$intron_coords,$intron_len,$iseqid,$array_introns)=get_introns($gene_name,'extract',$opt_feature_name);# if $num_exons > 1;
		my$num_introns=$array_introns?scalar@$array_introns:0; # if $array_introns;
		$intron_len||=0;
        if ($array_exons) {
            next if scalar@$array_exons < 1;
        ##print as: "GeneName","NumIntrons","gene_length","CDS_length","Intron_length","\%GC_exons","\%GC_introns","GC_each_exons","GC_each_introns","len_each_exon","len_each_intron","GC_per_10percent_of_CDS"
		### filehandles: $fh_gcperintron,$fh_gcperexon,$fh_lenperintron,$fh_lenperexon,$fh_perlen

            print {$fh} join "\t","\n$gene_name",$num_introns,$CDS_len+$intron_len,$CDS_len,$intron_len,GC_content_percent($CDS_Seq),GC_content_percent($intron_seq) if $CDS_Seq ne "0" ;

			 print {$fh_gcperexon} "\n$gene_name";
			 print {$fh_gcperintron} "\n$gene_name";
			 print {$fh_lenperexon} "\n$gene_name";
			 print {$fh_lenperintron} "\n$gene_name";
			 print {$fh_perlen} "\n$gene_name";
			 print {$fh_perlensliding} "\n$gene_name";





            print {$fh} "\t";
            my$GCexoneach;
            foreach (@$array_exons){print {$fh_gcperexon} "\t",$_?GC_content_percent($_):0}
			foreach (@$array_exons){print {$fh} $_?GC_content_percent($_):0,","}
            print {$fh} "\t";

            if (scalar@$array_introns > 0) {
                foreach (@$array_introns){print {$fh_gcperintron} "\t",$_?GC_content_percent($_):0}
                foreach (@$array_introns){print {$fh} $_?GC_content_percent($_):0,","}

            }else{
                print {$fh_gcperintron} "\t0";
                print {$fh} "0";
            }
            print {$fh} "\t";

            foreach (@$array_exons){print {$fh_lenperexon} "\t",$_?length($_):0}
			foreach (@$array_exons){print {$fh} $_?length($_):0,","}
            print {$fh} "\t";
            if (scalar@$array_introns > 0) {
                foreach (@$array_introns){print {$fh_lenperintron} "\t",$_?length($_):0}
                foreach (@$array_introns){print {$fh} $_?length($_):0,","}
            }else{
                print {$fh_gcperintron} "\t0";
                print {$fh} "0";
            }
            print {$fh} "\t";
            my$count=0;
            my$frag=10;
            for(my$i=0;  $i < $CDS_len;  $i=int($count * ($CDS_len/$frag))){

                my$extend=($count < $frag? int(($count+1) * $CDS_len/$frag) - $i :  $CDS_len-$i );
                die "Check files as extend value could not be 10 or lower. current extend:$extend\nCount:$count\t:CDSLen:$CDS_len StartCoord:$i" if $extend < 10;
                print "\nCount:$count\t:CDSLen:$CDS_len StartCoord:$i\tExtend:$extend\tEndCoord:", $i+$extend if $verbose;

                print {$fh_perlen} "\t",GC_content_percent(substr($CDS_Seq,  $i, $extend));
				print {$fh} GC_content_percent(substr($CDS_Seq,  $i, $extend)),$count<9?",":"";
                $count++;

            }
			### perl length sliding window. Calculations adjusts error due to fraction and tries to keep the window size and slide size similar along the length so that the last window does not get too small.
			$count=0;
			$frag=int((100-$window)/$slide);
			my$slide_len=$slide*$CDS_len/100;
			my$window_len=int($window*$CDS_len/100);
            for(my$i=0;  $i <= $CDS_len - $window_len;  $i=int($count * $slide_len)){

			my$extend=($count < $frag ? int($count * $slide_len + $window_len + 0.5) - $i  :  $CDS_len-$i );
                print "\nFrag:$frag\tCount:$count\t:CDSLen:$CDS_len StartCoord:$i\tExtend:$extend\tEndCoord:", $i+$extend if $verbose;

                print {$fh_perlensliding} "\t",GC_content_percent(substr($CDS_Seq,  $i, $extend));
				print {$fh} GC_content_percent(substr($CDS_Seq,  $i, $extend)),$count<$frag-1?",":"";
                $count++;

            }

        }


	}
	# to Write the gene bank file for sequence
	elsif(lc$opt_feature_type eq 'genbank'){die "\nThis method is not implemented yet\n"; $gene_seq=write_genebank_file($gene_name);print {$fh}"$gene_seq\n";}

	# to Extract the promoter sequence. upstream and downstream is with respect to Transcripton start site(gene start).
	elsif(lc$opt_feature_type eq 'promoter'){
		($gene_seq,$coords,$seqlen,$seqid)=extract_promoter($gene_name,$opt_upstream,$opt_downstream,$opt_feature_name);
		print {$fh}">$gene_name Promoter $opt_upstream upstream and $opt_downstream downstream of TSS $seqid$coords length:$seqlen\n$gene_seq\n" if $gene_seq;
	}
	# to Calculate the GC content of sequence upstream and downstream of repeat element. Need to provide GFF3 file produced by repeat masker or any other program.
	elsif(lc$opt_feature_type eq 'gc_around_repeats'){
		my$gc_content_around_repeats=GC_around_repeats($gene_name,$opt_upstream,$opt_downstream,$opt_run_modein_len_for_GC,$opt_feature_name,$opt_min_rep_len);
		push (@attributes_gcup_gcdown,@{$gc_content_around_repeats}) if $gc_content_around_repeats !~ /^\s*$/;
	}
    elsif(lc$opt_feature_type eq 'intron'){

        ($gene_seq,$coords,$seqlen,$seqid)=get_introns($gene_name,'extract',$opt_feature_name);
        print {$fh}">$gene_name $opt_feature_name:$seqid$coords length:$seqlen \n$gene_seq\n" if $gene_seq ne "0";
	}

	# or die
	else{die "Type of function ($opt_feature_type) is not a valid function type.Please choose function type from\ngene\nCDS\nExon\nmRNA\ntranscripts\nMask\nGenbank\ngc_around_repeats"}

	if($gene_seq && $gene_seq eq "0" && $opt_feature_type ne 'gc_around_repeats'){print "No gene feature found for $gene_name\n"}
    #else{print "\nGen_sequence:$gene_seq\n" };

}
##################################################################
# print the masked genome if mask option is selected
if(lc$opt_feature_type eq 'mask'){

    my@filename=split(/\//,$opt_seqfile);

    print "\nWriting masked sequences in file masked.$filename[-1]\n";
    open MASKOUT,">masked.$filename[-1]" or die "Could not create file to write masked genome\n";
    foreach(keys %$genome_seq){
        print MASKOUT">$_\n$$genome_seq{$_}\n";
    }

close MASKOUT;
}

##################################################################
# Summarize the table if gc_around_repeats option is opted

my%repeat_class;
undef %Gene; # empty memory as we dont need this hash anymore.
if(lc$opt_feature_type eq 'gc_around_repeats' && $opt_repeat_group_list ){

	open REPLIST, "$opt_repeat_group_list";
	while(<REPLIST>){
		my($repeat,$class)=split(/\s+/,$_);
		s/\s+//g foreach ($repeat,$class);
		$repeat_class{$repeat}=$class;
	}
}




if (lc$opt_feature_type eq 'gc_around_repeats'){

	# summarizing GC content data
	print "\nSummarizing GC content data....\n";
	summarize_gc_around_repeat(\@attributes_gcup_gcdown,$fh,\%repeat_class,$seq_len);
	print "\nGC content summary table is written in output file\n";

	# write raw table in output file
	if($print_raw){
		my$gc_around_repeats=join("\n",@attributes_gcup_gcdown) ;
		print "\nWriting GC content table in output file\n";
		print {$fh}"$gc_around_repeats" ;
	}
}



close $fh;


#Columns are seperated by TAB not spaces.
#Column 1: "seqid"
#Column 2: "source"
#Column 3: "type"
#Columns 4 & 5: "start" and "end"
#Column 6: "score"
#Column 7: "strand"
#Column 8: "phase"
#Column 9: "attributes"


##################################################################
# Subroutines used  in this program


###################################################################################################################
# Read fasta file and return hash of sequences
###################################################################################################################

sub ReadFasta{ # to read fasta format files into hash. returns hash.

	my $seqfile=shift(@_);
	my$demo_header;
	my$seq_len=0;
	my ($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "Reading Sequences from input file.....Plz wait...\n";
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

	print "Done....\nNumber of sequences read form input file = $seq_count\n\nExample of one the headers in sequence\n$demo_header\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.

	return(\%seq_hash,$seq_len);

}

##########################################################################################################
# Collect gene information abouthe the gene gene_name and store in global variable
##########################################################################################################
sub get_gene_info{

    my$gff_file=shift; # Get GFF file name
    my$gene_name=shift; # Get gene name to look for
	my$speed=shift;
    print "\nReading GFF file and extracting gene information.....";
    $gene_name=~s/\s+$//; #Remove trailing white spaces
    $gene_name=~s/^\s+//; # Remove heading white spaces
    chomp($gff_file,$gene_name);
    my@selected_gff_lines; # array to store lines related to gene name
    open GFF,$gff_file or die "Cannot find gff file $opt_GFF_file";
	#my$gene_name_new="\t$gene_name\t";
	# Pick all the lines in GFF file related to gene names.
  if($speed eq 'memory'){
		while (<GFF>){
			next if m/^\#/; # exclude comment lines
			next if m/^\s*$/; # exclude empty lines
			s/\n//; # Remove newline (enter) form the matched line
			s/^\s+//;
			if(/$gene_name\D/){push(@selected_gff_lines,$_)}
		}
		close GFF;

		# Create hash of all the elements for gene name identified from GFF
		foreach my$line(@selected_gff_lines){
			$line=~s/^\s+//g;
			#print "\nThis is from subroutine get_gene_info:$line";
			my($seqid,$source,$type,$start,$end,$score,$strand,$phase,$attributes,@attributes)=split(/\t/,$line);
			chomp($seqid,$source,$type,$start,$end,$score,$strand,$phase,$attributes,@attributes);
			s/^\s*//g foreach($seqid,$source,$type,$start,$end,$score,$strand,$phase);
			s/\s*$//g foreach($seqid,$source,$type,$start,$end,$score,$strand,$phase);
			$type=~s/\s+//g;
			push(@{$Gene{$gene_name}{$type}{'startend'}},$start,$end);
			push(@{$Gene{$gene_name}{$type}{'line'}},$line);

			$Gene{$gene_name}{$type}{'strand'}=$strand;
			$Gene{$gene_name}{$type}{'chr'}=$seqid;
			$Gene{$gene_name}{$type}{'gene_name'}=$gene_name;
			$Gene{$gene_name}{$type}{'attributes'}=join(" ",$attributes,@attributes) if ($type eq "mRNA");;

		}

	#sort values in each array.
	@{$Gene{$gene_name}{$_}{'startend'}}= sort { $a <=> $b } @{$Gene{$gene_name}{$_}{'startend'}} for keys %{$Gene{$gene_name}};
	# %Gene is a global variable so no need to use return()


  }
  else{
		my$count;
        my$geneid;
        my$ID;
        my$parent;
		while (<GFF>){


			next if m/^\#/; # exclude comment lines
			next if m/^\s*$/; #exclude empty lines
			s/\n//; # Remove newline (enter) form the matched line
			s/^\s+//;
			my$line=$_;
			#print "\nThis is from subroutine get_gene_info:$line";
			my($seqid,$source,$type,$start,$end,$score,$strand,$phase,$attributes,@attributes)=split(/\t/,$line);
			chomp($seqid,$source,$type,$start,$end,$score,$strand,$phase,$attributes,@attributes);
			s/^\s*//g foreach($seqid,$source,$type,$start,$end,$score,$strand,$phase);
			s/\s*$//g foreach($seqid,$source,$type,$start,$end,$score,$strand,$phase);
			$type=~s/\s+//g;

            #$geneid=$1 if $attributes=~/ID=([^;]+)/ ;#&& (lc$type eq 'gene' || lc$type eq 'mrna');
            #$geneid=$1 if $attributes=~/Parent=([^;]+)/ ;#&& (lc$type eq 'gene' || lc$type eq 'mrna');


            if($attributes=~/Parent=([^;]+)/){
                $parent=$1;
                $geneid=$parent;
            }

            if($attributes=~/ID=([^;]+)/){
                $ID=$1;
                $geneid=$ID if (lc$type eq 'gene' || lc$type eq 'mrna' || lc$type eq 'transcript');
            }




            $geneid=~s/^\s*|\s*$//g;
            ## if element belongs to more than one parents, assign it to each parent
            my@parents=split /,/,$geneid;
            foreach my $par(@parents){
                next if $par =~ /^\s*$/;
    			push(@{$Gene{$par}{$type}{'startend'}},$start,$end);
    			push(@{$Gene{$par}{$type}{'line'}},$line);
    			$Gene{$par}{$type}{'strand'}=$strand;
    			$Gene{$par}{$type}{'chr'}=$seqid;
    			$Gene{$par}{$type}{'gene_name'}=$par;
                $Gene{$par}{$type}{'parent'}=$parent;
    			$Gene{$par}{$type}{'attributes'}=join(" ",$attributes,@attributes) if ($type eq "mRNA");;
            }
			$count++;
		}

		#sort values in each array and remove elements smaller than defined length if asked..
		foreach my$seqid(keys %Gene){
			@{$Gene{$seqid}{$_}{'startend'}}= sort { $a <=> $b } @{$Gene{$seqid}{$_}{'startend'}} for keys %{$Gene{$seqid}};
			# %Gene is a global variable so no need to use return()
            if ($len_lim) {
                foreach my$typ(keys %{$Gene{$seqid}}){
                    my$feature_len=${$Gene{$seqid}{$typ}{'startend'}}[-1] - ${$Gene{$seqid}{$typ}{'startend'}}[0] + 1;
                    delete $Gene{$seqid}{$typ} if ($feature_len < $len_lim);
                }
            }
		}
		print "\nRead coordinates of  $count elements from GFF file\n";


  }


}

##########################################################################################################
# Get the sequence of name gene_name
##########################################################################################################
sub get_seq{

    my $gene_name=shift;
	my $do_what=shift; # either 'mask' or 'extract'
	my $feature=shift; # feature name

	chomp($gene_name,$do_what,$feature);


		# get gene info from %Gene (passed as reference) from global variable
        my $gene_seq=();
        my @exons;
        my $gene_start=${$Gene{$gene_name}{$feature}{'startend'}}[0];
        my $gene_end=${$Gene{$gene_name}{$feature}{'startend'}}[-1];
        my $seq_id=$Gene{$gene_name}{$feature}{'chr'};
        my $gene_strand=$Gene{$gene_name}{$feature}{'strand'};
        my $CDS_coord='(';

	if(!$gene_start || $gene_start<=0){return 0;} # return zero if start position of the feature is not found.


	if(exists $$genome_seq{$seq_id}){
            #if($do_what2 eq 'mask'){
            print "This gene region is being masked in $seq_id\n" if $do_what eq 'mask' ;
            for(my$i=0;$i<scalar@{$Gene{$gene_name}{$feature}{'startend'}} - 1; $i=$i+2){

				my $new_substring=();
				my$length=${$Gene{$gene_name}{$feature}{'startend'}}[$i+1] - ${$Gene{$gene_name}{$feature}{'startend'}}[$i] + 1;
				if($do_what eq 'mask'){$new_substring=substr($$genome_seq{$seq_id}, ${$Gene{$gene_name}{$feature}{'startend'}}[$i] - 1, $length,"N"x$length);$gene_seq .=$new_substring}
				else{$new_substring=substr($$genome_seq{$seq_id},${$Gene{$gene_name}{$feature}{'startend'}}[$i] - 1, $length );$gene_seq .=$new_substring}
				$CDS_coord.=${$Gene{$gene_name}{$feature}{'startend'}}[$i].'..'.${$Gene{$gene_name}{$feature}{'startend'}}[$i+1].',';
                push(@exons,$new_substring);
            }
    }
        else{print "For gene $gene_name on Sequence with name $seq_id does not exists in the genomic file\n"; return 0;}

        $CDS_coord.=')';
        $CDS_coord=~s/,\)/\)/; # remove last comma

        my$gene_length=length$gene_seq;


        #print "Seq id:$seq_id\n";
        #print "Exon Start:$gene_start\n";
        #print "Exon end:$gene_end\n";
        #print "Exon strand:$gene_strand\n";
        #print "Exon Length:$gene_length\n";
        #print "Exon Coordinates:$CDS_coord\n";
        #my$gene_seq=substr($$genome_seq{$seq_id},$gene_start-1,$gene_end-$gene_start);
        $gene_seq=revcomp($gene_seq) if  $gene_strand eq "-";

        return ($gene_seq,$CDS_coord,$gene_length,$seq_id,\@exons);


}

##########################################################################################################
# Extract promoter region of gene_name
##########################################################################################################

sub extract_promoter{

    my$gene_name=shift;
	my$upstream=shift;
	my$downstream = shift;
    my$feature=shift;
    ### next if promoter is asked for non-gene or non-mRNA region.
    return undef if (!$Gene{$gene_name}{'gene'}{'startend'} && !$Gene{$gene_name}{'gene'}{'startend'});
    # get gene info from %Gene (passed as reference) from global variable
    my$gene_seq=();
    my$gene_start=${$Gene{$gene_name}{$feature}{'startend'}}[0];
    my$gene_end=${$Gene{$gene_name}{$feature}{'startend'}}[-1];
    my$seq_id=$Gene{$gene_name}{$feature}{'chr'};
    my$gene_strand=$Gene{$gene_name}{$feature}{'strand'};
    my$seq_length=length($$genome_seq{$seq_id});
    my$CDS_coord='(';
    my$promoter_start=0;
    my$promoter_end=0;
    my$extraction_length=0;

    # Calculate co-ordinates of sequence to be extracted upstream and downstream of TSS
    if($gene_strand eq '-' || lc$gene_strand eq 'minus' ){

	$promoter_end=$gene_end+$upstream; # $Gene{$gene_name}{$feature}{'startend'} is sorted so $gene_end is actually gene start for minus strand.
	$promoter_end=$seq_length if $promoter_end > $seq_length;

	$promoter_start=$gene_end-$downstream;
	$promoter_start=1 if $promoter_start<=0;

	$extraction_length=$promoter_end-$promoter_start+1;
    }
    else{

	print "Error: The strand of Feature is not defined. Current value for strand:$gene_strand\t....Assuming 'positive' strand.\n" if ($gene_strand ne '+' && $gene_strand ne 'plus') ;

	$promoter_start=$gene_start-$upstream;
	$promoter_start=1 if $promoter_start <= 0;

	$promoter_end=$gene_start+$downstream;
	$promoter_end=$seq_length if $promoter_end>$seq_length;

	$extraction_length=$promoter_end-$promoter_start+1;

    }

    # extract promoter region.
    $gene_seq =substr($$genome_seq{$seq_id},$promoter_start - 1, $extraction_length);
    $gene_seq=revcomp($gene_seq) if ($gene_strand eq '-' || lc$gene_seq eq 'minus');


    $CDS_coord.=$promoter_start.'..'.$promoter_end.')';
    $CDS_coord=~s/,\)/\)/; # replace ',)' with ')'

    my$promoter_length=length$gene_seq;

    print "Extracting promoter from $seq_id for gene name:$gene_name\tfor feature:$feature\tPromot_start:$promoter_start \tPromo_end:$promoter_end\tExtracted_length:$promoter_length\n" if $verbose;

    return ($gene_seq,$CDS_coord,$promoter_length,$seq_id);




}

##########################################################################################################
# Calculate GC around repeats. return the attributes and GC content or upstream and downstream region
##########################################################################################################

sub GC_around_repeats{

    my$sequence_name = shift;
	my$upstream = shift;
	my$downstream = shift;
	my$min_len=shift;
	my$feature=shift; # look for this term when reading co-ordinates for features
	my$min_element_len=shift;
	if(!$min_element_len){$min_element_len=10}

	my@GC_content;
	return(" ") if !$Gene{$sequence_name}{$feature}{'line'}; # return empty string if the array referenc is not defined
	my$number=@{$Gene{$sequence_name}{$feature}{'line'}};
	#print "\nFound $number Repeat elements for $sequence_name\n";
	#print "\n\n***********************Called subroutine GC_around_repeats for sequence:$sequence_name\n";
	# sort the array based on the start site of elements on sequence. Assuming no overlapping elements
	sort_gff3_lines(\@{$Gene{$sequence_name}{$feature}{'line'}},4);


	#foreach my$line(@{$Gene{$sequence_name}{$feature}{'line'}}){ # whole line is stored in global variable %Gene{genename}{type}{'line'} for each sequence_name
		#my($seqid,$source,$type,$start,$end,$score,$strand,$phase,@attributes)=split(/\t/,$line);
        #chomp($seqid,$source,$type,$start,$end,$score,$strand,$phase,@attributes);
		#my$attributes=join(" ",@attributes);
		#if only one element in the  the fragment (seqid) extract up and downstream sequence
		if($number==1){
			#print "\n*********I am in the loop to process the only element\n";
			# split the line
			my($seqid,$source,$type,$start,$end,$score,$strand,$phase,@attributes)=split(/\t/,${$Gene{$sequence_name}{$feature}{'line'}}[0]);
			chomp($seqid,$source,$type,$start,$end,$score,$strand,$phase,@attributes);
			my$attributes=join(" ",@attributes);
			#$attributes=~s/\s+[\s\d]+$//g;
			#$attributes=~s/Target=//g;
			clean_comment(\$attributes);
			#extract upstream
			my$upstream_seq=extract_upstream($start,$end,$strand,$upstream,$$genome_seq{$seqid});

			#extract downstream
			my$downstream_seq=extract_downstream($start,$end,$strand,$downstream,$$genome_seq{$seqid});


			# filter for minimum length specified by user. Dont calc GC if length is less than minimum allowed (by $min_len)
			my$GCup='ShortSeq';
			my$GCdown='ShortSeq';

			if (num_ATGC($upstream_seq)>=$min_len){$GCup=GC_content($upstream_seq)}
			if (num_ATGC($downstream_seq)>=$min_len){$GCdown=GC_content($downstream_seq)}
			#extract repeat_element
			my$element_seq=extract_element($start,$end,$strand,$$genome_seq{$seqid});
			my$GCelement='ShortSeq';
			if (num_ATGC($element_seq)>=$min_element_len){$GCelement=GC_content($element_seq)}

			#my$GCcontentLine=join("\t",$attributes,$GCelement,$GCup,$GCdown);
			my$GCcontentLine=join("\t",$attributes,length($element_seq),$GCelement,length($upstream_seq),$GCup,length($downstream_seq),$GCdown);
			push(@GC_content,$GCcontentLine);

			#print "\n**********End of the loop to process only element\n Values at the end of the loop are:\nstrand:$strand\nStart:$start\nEnd:$end\nUpstream:$upstream\nDownStream:$downstream\n";

			#print "\nExtracting upstreamseq with following info\nStart:$start\nEnd:$end\nStrand:$strand\nUpstream:$upstream\nSequence:$$genome_seq{$seqid}\nSequence length:".length($$genome_seq{$seqid})."\nExtracted Sequence:$upstream_seq\n";
			#print "\nExtracting Downstreamseq with following info\nStart:$start\nEnd:$end\nStrand:$strand\nDownstream:$downstream\nSequence:$$genome_seq{$seqid}\nSequence length:".length($$genome_seq{$seqid})."\nExtracted sequence;$downstream_seq\n";




		}
		elsif($number >1){ # if more than one element are in one sequence, make sure they dont get counted in in GC calculation.
			# split the information for each element in each line and store them in individual arrays so they can be retrived when needed.
			my(@seqid,$seqid,@source,$source,@type,$type,@start,$start,@end,$end,@score,$score,@strand,$strand,@phase,$phase,@attribute,@attributes,$seqlength,@upstream,@downstream);
			for(my$i=0;$i<$number;$i++){
				($seqid[$i],$source[$i],$type[$i],$start[$i],$end[$i],$score[$i],$strand[$i],$phase[$i],@attribute)=split(/\t/,${$Gene{$sequence_name}{$feature}{'line'}}[$i]);
				#chomp($seqid,$source,$type,$start,$end,$score,$strand,$phase,@attributes);
				s/^\s*// foreach($seqid[$i],$source[$i],$type[$i],$start[$i],$end[$i],$score[$i],$strand[$i],$phase[$i]);
				s/\s*$// foreach($seqid[$i],$source[$i],$type[$i],$start[$i],$end[$i],$score[$i],$strand[$i],$phase[$i]);
				$upstream[$i]=$upstream; # each should have their own values of upstream and downstream other wise need to change for one will affect everyone.
				$downstream[$i]=$downstream;# each should have their own values of upstream and downstream other wise need to change for one will affect everyone.
				$attributes[$i]=join(" ",@attribute);
				#$attributes[$i]=~s/Target=//g;
				clean_comment(\$attributes[$i]);
				$seqlength=length $$genome_seq{$seqid[$i]};
				#push(@attributes_all,$attributes);
			}




			# calculate new upstream and downstream values based on how far the next element is. Purpose is not to include the next element in GC calculation.
			# First and last element will be processed seperately as they dont have any element in at least ones side. rest will be processed seperately.
			# This algorithm assumes that elements lines are sorted based on their location on the sequence. It will mess up if it is not sorted.
			for(my$i=0;$i<$number;$i++){


				if($i==0){
					#print "\n**********I am in loop to process first element\n";
					#for the first element. as in following example
					# --------start[$i]>>>>>>>>>>>>>>>>>>>>>>end[$i]-----------------start[$i+1]>>>>>>>>>>>>>>>>>>>>>>>>end[$i+1]-------------------------end[$i+2]<<<<<<<
					# for positive strand. requested Upstream will be till the start of the sequence if start of the sequence is smaller than reqested upstream.
					# for positive strand. requested Downstream will be till the start of the next element if distance between is smaller than reqested downstreamstream.
					$upstream[$i]=(min($start[$i],$end[$i])-1) if (($strand[$i] eq '+'||lc$strand[$i] eq 'plus') &&
											 (min($start[$i],$end[$i])<$upstream[$i]));

					$downstream[$i]=(min($start[$i+1],$end[$i+1])-max($start[$i],$end[$i])-1) if (($strand[$i] eq '+'||lc$strand[$i] eq 'plus') &&
											 ((min($start[$i+1],$end[$i+1])-max($start[$i],$end[$i])-1)<$downstream[$i]));


					# for negative strand as in following example
					# --------end[$i]<<<<<<<<<<<<<<<<<<<<<<start[$i]-----------------start[$i+1]>>>>>>>>>>>>>>>>>>>>>>>>end[$i+1]-------------------------end[$i+2]<<<<<<<
					# for negative strand. requested downstream will be till the start of the sequence if start of the sequence is smaller than reqested downstream.
					# for negative strand. requested upstream will be till the start of the next element if distance between is smaller than reqested upstreamstream.
					$downstream[$i]=(min($start[$i],$end[$i]) -1)if (($strand[$i] eq '-'||lc$strand[$i] eq 'minus') &&
											 (min($start[$i],$end[$i]-1)<$downstream[$i]));

					$upstream[$i]=(min($start[$i+1],$end[$i+1])-max($start[$i],$end[$i])-1) if (($strand[$i] eq '-'||lc$strand[$i] eq 'minus') &&
											 ((min($start[$i+1],$end[$i+1])-max($start[$i],$end[$i]-1)<$upstream[$i])));
					#print "\n*****************End of the loop to process First element\nValues at the end of the loop are:\nstrand:$strand[$i]\nStart:$start[$i]\nEnd:$end[$i]\nUpstream:$upstream[$i]\nDownStream:$downstream[$i]\n";

				}
				# for the last element
				elsif($i==$number-1){
					#print "\n***************I am in loop to process last element\n";
					# for positive strand. as in following example
					# <<<<<start[$i-2]------start[$i-1]>>>>>>>>>>>>>>>>>>>>>>end[$i-1]-----------------start[$i]>>>>>>>>>>>>>>>>>>>>>>>>end[$i]------------------
					# for positive strand. requested Downstream will be till the end of the sequence if distance till the end of the sequence is smaller than reqested downstream.
					# for positive strand. requested upstream will be till the end of the previous element if distance between is smaller than reqested upstreamstream.
					$upstream[$i]=(min($start[$i],$end[$i])-max($start[$i-1],$end[$i-1])-1)if (($strand[$i] eq '+'||lc$strand[$i] eq 'plus') &&
																						((min($start[$i],$end[$i])-max($start[$i-1],$end[$i-1])-1)<$upstream[$i]));

					$downstream[$i]=($seqlength-max($start[$i],$end[$i])) if (($strand[$i] eq '+'||lc$strand[$i] eq 'plus') &&
																							(($seqlength-max($start[$i],$end[$i]))<$downstream[$i]));

					#print "\nCalculated Downstream value for plus strand:$downstream[$i]\nSeqlength:$seqlength\nMax of start and end:".max($start[$i],$end[$i])."\nand difference is :".($seqlength-max($start[$i],$end[$i])) ."\n" if ($strand[$i] eq '+'||lc$strand[$i] eq 'plus');
					# for negative strand as in following example
					# <<<<<start[$i-2]------start[$i-1]>>>>>>>>>>>>>>>>>>>>>>end[$i-1]-----------------end[$i]<<<<<<<<<<<<<<<<<<<<<<<<start[$i]------------------
					# for negative strand. requested upstream will be till the end of the sequence if distance till the end of the sequence is smaller than reqested upstream.
					# for negative strand. requested downstream will be till the end of the previous element if distance between is smaller than reqested downstream.
					$downstream[$i]=(min($start[$i],$end[$i])-max($start[$i-1],$end[$i-1])-1)if (($strand[$i] eq '-'||lc$strand[$i] eq 'minus') &&
																						((min($start[$i],$end[$i])-max($start[$i-1],$end[$i-1])-1)<$downstream[$i]));

					$upstream[$i]=$seqlength-max($start[$i],$end[$i]) if (($strand[$i] eq '-'||lc$strand[$i] eq 'minus') &&
																							($seqlength-max($start[$i],$end[$i])<$upstream[$i]));
					#print "\nCalculated Downstream value for minus strand :$downstream[$i]\nSeqlength:$seqlength\nMax of start and end:".max($start[$i],$end[$i])."\nand difference is :".($seqlength-max($start[$i],$end[$i])) ."\n" if ($strand[$i] eq '-'||lc$strand[$i] eq 'minus');
					#print "\n***************End of the loop to process last element\n Values at the end of the loop are:\nstrand:$strand[$i]\nStart:$start[$i]\nEnd:$end[$i]\nUpstream:$upstream[$i]\nDownStream:$downstream[$i]\n";

				}
				# for the rest of the element in between first and last.
				else{
					#print "\n*****************I am in loop to process".($i+1)." element\n";
					# for positive strand. as in following example
					# <<<<<start[$i-1]------start[$i]>>>>>>>>>>>>>>>>>>>>>>end[$i]-----------------start[$i+1]>>>>>>>>>>>>>>>>>>>>>>>>end[$i+1]------------------
					# for positive strand. requested Downstream will be till the start of the nest element if distance between them is smaller than reqested downstream.
					# for positive strand. requested upstream will be till the end of the previous element if distance between them is smaller than reqested upstreamstream.
					$upstream[$i]=(min($start[$i],$end[$i])-max($start[$i-1],$end[$i-1])-1) if (($strand[$i] eq '+'||lc$strand[$i] eq 'plus') &&
																						((min($start[$i],$end[$i])-max($start[$i-1],$end[$i-1])-1)<$upstream[$i]));


					$downstream[$i]=(min($start[$i+1],$end[$i+1])-max($start[$i],$end[$i])-1) if (($strand[$i] eq '+'||lc$strand[$i] eq 'plus') &&
																							((min($start[$i+1],$end[$i+1])-max($start[$i],$end[$i])-1)<$downstream[$i]));


					# for negative strand as in following example
					# -------end[$i-1]<<<<<start[$i-1]------end[$i]<<<<<<<<<<<<start[$i]-----------------end[$i+1]<<<<<<<<<<<<<<<<<<start[$i+1]------------------
					# for negative strand. requested upstream will be till the end of the next element if distance between them is smaller than reqested upstream.
					# for negative strand. requested downstream will be till the end of the previous element if distance between is smaller than reqested downstream.
					$downstream[$i]=(min($start[$i],$end[$i])-max($start[$i-1],$end[$i-1])-1)if (($strand[$i] eq '-'||lc$strand[$i] eq 'minus') &&
																						((min($start[$i],$end[$i])-max($start[$i-1],$end[$i-1])-1)<$downstream[$i]));

					$upstream[$i]=(min($start[$i+1],$end[$i+1])-max($start[$i],$end[$i])-1) if (($strand[$i] eq '-'||lc$strand[$i] eq 'minus') &&
																							((min($start[$i+1],$end[$i+1])-max($start[$i],$end[$i])-1)<$upstream[$i]));


					#print "\n******************End of the loop to process".($i+1)." element\n Values at the end of the loop are:\nstrand:$strand[$i]\nStart:$start[$i]\nEnd:$end[$i]\nUpstream:$upstream[$i]\nDownStream:$downstream[$i]\n";


				}


				##extract upstream
				my$upstream_seq=extract_upstream($start[$i],$end[$i],$strand[$i],$upstream[$i],$$genome_seq{$seqid[$i]});
				#print "\nExtracting upstreamseq with following info\nStart:$start[$i]\nEnd:$end[$i]\nStrand:$strand[$i]\nUpstream:$upstream[$i]\nSequence:$$genome_seq{$seqid[$i]}\nSequence length:$seqlength\nExtracted Sequence:$upstream_seq\n";
				#extract downstream
				my$downstream_seq=extract_downstream($start[$i],$end[$i],$strand[$i],$downstream[$i],$$genome_seq{$seqid[$i]});
				#print "\nExtracting Downstreamseq with following info\nStart:$start[$i]\nEnd:$end[$i]\nStrand:$strand[$i]\nDownstream:$downstream[$i]\nSequence:$$genome_seq{$seqid[$i]}\nSequence length:$seqlength\nExtracted Sequence:$downstream_seq\n";
				#save attributes and GC content info for upstream and downstream in an array.
				#$attributes[$i]=~s/\s+[\s\d]+$//g;
				#$attributes[$i]=~s/Target=//g;
				#$attributes[$i]=~s/"//g;
				#$attributes[$i]=~s/Target Motif://g;
				clean_comment(\$attributes[$i]);

				# filter for minimum length specified by user. Dont calc GC if length is less than minimum allowed (by $opt_run_modein_len_for_GC)
				my$GCup='ShortSeq';
				my$GCdown='ShortSeq';
				if (num_ATGC($upstream_seq)>=$min_len){$GCup=GC_content($upstream_seq)}
				if (num_ATGC($downstream_seq)>=$min_len){$GCdown=GC_content($downstream_seq)}

				#extract repeat_element
				my$element_seq=extract_element($start[$i],$end[$i],$strand[$i],$$genome_seq{$seqid[$i]});
				my$GCelement='ShortSeq';
				if (num_ATGC($element_seq)>=$min_element_len){$GCelement=GC_content($element_seq)}

				#my$GCcontentLine=join("\t",$attributes[$i],$GCelement,$GCup,$GCdown);
				my$GCcontentLine=join("\t",$attributes[$i],length($element_seq),$GCelement,length($upstream_seq),$GCup,length($downstream_seq),$GCdown);

				#my$GCcontentLine=join("\t",$attributes[$i],GC_content($upstream_seq),GC_content($downstream_seq));
				push(@GC_content,$GCcontentLine);
			}
		}
	#}


	return(\@GC_content);

}

##############################################################################################
# sort array containing GFF lines based on the column requested by user or column 4 (start)
##############################################################################################


sub sort_gff3_lines{
	#print "****************************For each reference:\n";
	my$ref_gff3lines_Array=shift;
	my$which_column_to_sort_with=shift;
	$which_column_to_sort_with=4 if !$which_column_to_sort_with;

	my(%stored_elements,@list_of_sorting_values,@sorted_gff_array);

	foreach my$line(@{$ref_gff3lines_Array}){
		my@elements=split(/\t/,$line);
		s/^\s+//g foreach @elements;
		s/\s+$//g foreach @elements;
		$stored_elements{$elements[$which_column_to_sort_with-1]}=$line;
		push(@list_of_sorting_values,$elements[$which_column_to_sort_with-1]);
	}

	@list_of_sorting_values=sort{$a<=>$b} @list_of_sorting_values;

	foreach(@list_of_sorting_values){

		push(@sorted_gff_array,$stored_elements{$_})

	}

	#print"\nBefore sorting\n".join("\n",@{$ref_gff3lines_Array});
	# replace the value of provided array with new sorted array.
	@{$ref_gff3lines_Array}=@sorted_gff_array;



	#print"\nAfter sorting\n" .join("\n",@{$ref_gff3lines_Array});
}

##############################################################################################
# write gene bank file using all the features from gff
##############################################################################################


sub write_genebank_file{
 die "\n\ngenbank option  is not implemented yet.\n\n";

}

###################################################################################################################
# return reverse compliment of sequence. Only ATGCN are reversed. Anything else will not be changed
###################################################################################################################

sub revcomp{
    my$seq=shift;

    $seq = reverse $seq;
    $seq =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq;
}

##################################################################################################################################
# Extract upstream sequence using information about feature start, feature end, $strand, upstream length and sequence
##################################################################################################################################
sub extract_upstream{

	my$feature_start=shift;
	my$feature_end=shift;
	my$feature_strand=shift;
	my$upstream=shift;
	my$sequence=shift;
	$upstream=0 if $upstream <0;
	my($extract_end,$extract_start,$extraction_length);
	my$seq_length=length($sequence);
    if($feature_strand eq '-' || lc$feature_strand eq 'minus' ){

		$extract_end=max($feature_end,$feature_start)+$upstream; # $Gene{$gene_name}{$feature}{'startend'} is sorted so $gene_end is actually gene start.
		$extract_end=$seq_length if $extract_end > $seq_length;
		$extract_start=max($feature_end,$feature_start);
	}
    else{

		print "Error: The strand of Feature is not defined. Current value for strand:$feature_strand\t....Assuming 'positive' strand.\n" if ($feature_strand ne '+' && $feature_strand ne 'plus') ;

		$extract_start=min($feature_end,$feature_start)-$upstream-1;
		$extract_start=0 if $extract_start <= 0;
		$extract_end=min($feature_end,$feature_start)-1;

    }

	$extraction_length=$extract_end-$extract_start;
	$extraction_length=0 if $extraction_length <0;
	my$upstream_seq=substr($sequence,$extract_start,$extraction_length);
	return($upstream_seq);
}

##########################################################################################################################################
# Extract downstream sequence using information about feature start, feature end, $strand, downstream length and sequence
##########################################################################################################################################
sub extract_downstream{ # to extract downstream of a feature. need start , end, downstream,strand and sequence

	my$feature_start=shift;
	my$feature_end=shift;
	my$feature_strand=shift;
	my$downstream=shift;
	my$sequence=shift;
	$downstream=0 if $downstream <0;
	my($extract_end,$extract_start,$extraction_length);
	my$seq_length=length($sequence);
	# Calculate co-ordinates of sequence to be extracted upstream and downstream of TSS
    if($feature_strand eq '-' || lc$feature_strand eq 'minus' ){

		$extract_start=min($feature_end,$feature_start)-$downstream-1; # $Gene{$gene_name}{$feature}{'startend'} is sorted so $gene_end is actually gene start.
		$extract_start=0 if $extract_start < 0;
		$extract_end=min($feature_end,$feature_start)-1;

    }
    else{

		print "Error: The strand of Feature is not defined. Current value for strand:$feature_strand\t....Assuming 'positive' strand.\n" if ($feature_strand ne '+' && $feature_strand ne 'plus') ;

		$extract_start=max($feature_end,$feature_start);
		$extract_end=max($feature_end,$feature_start)+$downstream;
		$extract_end=$seq_length if $extract_end>$seq_length;
	}

	$extraction_length=$extract_end-$extract_start;
	$extraction_length=0 if $extraction_length <0;
	my$downstream_seq=substr($sequence,$extract_start,$extraction_length);

	return($downstream_seq);
}


##########################################################################################################################################
# Extract downstream sequence using information about feature start, feature end, $strand, downstream length and sequence
##########################################################################################################################################
sub extract_element{ # to extract downstream of a feature. need start , end, downstream,strand and sequence

	my$feature_start=shift;
	my$feature_end=shift;
	my$feature_strand=shift;
	my$sequence=shift;
	my$seq_length=length($sequence);
	my$feature_length=max($feature_start,$feature_end)-min($feature_start,$feature_end)+1;
	my$feature_seq=substr($sequence,min($feature_start,$feature_end)-1,$feature_length);
	# Calculate co-ordinates of sequence to be extracted upstream and downstream of TSS
    if($feature_strand eq '-' || lc$feature_strand eq 'minus' ){
		$feature_seq=revcomp($feature_seq);
    }
	return($feature_seq);
}

##############################################################################################
# extract introns on if CDS coordinates are provided.
##############################################################################################
sub get_introns{

    my$gene_name=shift;
	my$do_what=shift; # either 'mask' or 'extract'
	my$feature=shift; # feature name

	chomp($gene_name,$do_what,$feature);


		# get gene info from %Gene (passed as reference) from global variable
        my$gene_seq="";
        my@introns;
        my $gene_start=${$Gene{$gene_name}{$feature}{'startend'}}[0];
        my$gene_end=${$Gene{$gene_name}{$feature}{'startend'}}[-1];
        my$seq_id=$Gene{$gene_name}{$feature}{'chr'};
        my$gene_strand=$Gene{$gene_name}{$feature}{'strand'};
        my$CDS_coord='(';

	if(!$gene_start || $gene_start<=0){return ($gene_seq,"NA",0,$seq_id,[]);} # return zero if start position of the feature is not found.


	if(exists $$genome_seq{$seq_id}){
            #if($do_what2 eq 'mask'){
            print "This gene region is being masked in $seq_id\n" if $do_what eq 'mask' ;
            for(my$i=1;$i<scalar@{$Gene{$gene_name}{$feature}{'startend'}} - 1; $i=$i+2){

				my $new_substring=();
				my$length=${$Gene{$gene_name}{$feature}{'startend'}}[$i+1] - ${$Gene{$gene_name}{$feature}{'startend'}}[$i] + 1;
				if($do_what eq 'mask'){$new_substring=substr($$genome_seq{$seq_id}, ${$Gene{$gene_name}{$feature}{'startend'}}[$i] - 1, $length,"N"x$length);$gene_seq .=$new_substring}
				else{$new_substring=substr($$genome_seq{$seq_id},${$Gene{$gene_name}{$feature}{'startend'}}[$i] - 1, $length );$gene_seq .=$new_substring}
				$CDS_coord.=${$Gene{$gene_name}{$feature}{'startend'}}[$i].'..'.${$Gene{$gene_name}{$feature}{'startend'}}[$i+1].',';
                push(@introns,$new_substring);
            }
    }
        else{print "For gene $gene_name on Sequence with name $seq_id does not exists in the genomic file\n"; return 0;}

        $CDS_coord.=')';
        $CDS_coord=~s/,\)/\)/; # remove last comma

        my$gene_length=length$gene_seq;
        return ($gene_seq,"NA",0,$seq_id,[]) if !$gene_seq;

        #print "Seq id:$seq_id\n";
        #print "Exon Start:$gene_start\n";
        #print "Exon end:$gene_end\n";
        #print "Exon strand:$gene_strand\n";
        #print "Exon Length:$gene_length\n";
        #print "Exon Coordinates:$CDS_coord\n";
        #my$gene_seq=substr($$genome_seq{$seq_id},$gene_start-1,$gene_end-$gene_start);
        $gene_seq=revcomp($gene_seq) if  $gene_strand eq "-";

        return ($gene_seq,$CDS_coord,$gene_length,$seq_id,\@introns);


}



##############################################################################################
# Return summary of GC content around repeats
##############################################################################################
sub summarize_gc_around_repeat{

	use Statistics::Descriptive;
	my$Ref_GC_array=shift;
	my$file_handle_to_print=shift;
	my$ref_repClass_hash=shift;
	my$total_len=shift;

	my%gc_array_summary;
	my%gc_array_group_summary;

	#progress bar initials
	my$count=my$next_value=0;
	my$total_lines=@{$Ref_GC_array};
	#@{$Ref_GC_array} lines contain $attributes,length($element_seq),$GCelement,length($upstream_seq),$GCup,length($downstream_seq),$GCdown);
	print "\n\t->Collecting datapoints and creating groups......\n";
	syswrite STDOUT, "\[";
	foreach my$line(@{$Ref_GC_array}){
		#------------------------------------------------------------------------------------
		# progress bar
		$count++;
		if(($count*100/$total_lines)>$next_value){syswrite STDOUT,"=";$next_value+=1; syswrite STDOUT,"$next_value\%" if $next_value%10==0;}
		#------------------------------------------------------------------------------------
		my($name,$len_element,$gcelement,$len_up,$gcup,$len_down,$gcdown)=split(/\t/,$line);
		s/\s+//g foreach ($name,$len_element,$gcelement,$len_up,$gcup,$len_down,$gcdown);

		push(@{$gc_array_summary{$name}{'gcelement'}},$gcelement) if $gcelement=~/[\d\.]+/;
		push(@{$gc_array_summary{$name}{'gcup'}},$gcup) if $gcup=~/[\d\.]+/;
		push(@{$gc_array_summary{$name}{'gcdown'}},$gcdown) if $gcdown=~/[\d\.]+/;
		push(@{$gc_array_summary{$name}{'gcupdown'}},$gcup) if $gcup=~/[\d\.]+/;
		push(@{$gc_array_summary{$name}{'gcupdown'}},$gcdown) if $gcdown=~/[\d\.]+/;
		# process length
		push(@{$gc_array_summary{$name}{'len_element'}},$len_element) if $len_element=~/[\d\.]+/;
		#push(@{$gc_array_summary{$name}{'len_up'}},$len_up) if $len_up=~/[\d\.]+/;
		#push(@{$gc_array_summary{$name}{'len_down'}},$len_down) if $len_down=~/[\d\.]+/;
		#push(@{$gc_array_summary{$name}{'len_updown'}},$len_up) if $len_up=~/[\d\.]+/;
		#push(@{$gc_array_summary{$name}{'len_updown'}},$len_down) if $len_down=~/[\d\.]+/;

		# if name-class hash element contain key with $name, push the values into an array and connect with the key
		if($$ref_repClass_hash{$name}){
			push(@{$gc_array_group_summary{$$ref_repClass_hash{$name}}{'gcelement'}},$gcelement) if $gcelement=~/[\d\.]+/;
			push(@{$gc_array_group_summary{$$ref_repClass_hash{$name}}{'gcup'}},$gcup) if $gcup=~/[\d\.]+/;
			push(@{$gc_array_group_summary{$$ref_repClass_hash{$name}}{'gcdown'}},$gcdown) if $gcdown=~/[\d\.]+/;
			push(@{$gc_array_group_summary{$$ref_repClass_hash{$name}}{'gcupdown'}},$gcup) if $gcup=~/[\d\.]+/;
			push(@{$gc_array_group_summary{$$ref_repClass_hash{$name}}{'gcupdown'}},$gcdown) if $gcdown=~/[\d\.]+/;

			# process length
			push(@{$gc_array_group_summary{$$ref_repClass_hash{$name}}{'len_element'}},$len_element) if $len_element=~/[\d\.]+/;
			#push(@{$gc_array_group_summary{$$ref_repClass_hash{$name}}{'len_up'}},$len_up) if $len_up=~/[\d\.]+/;
			#push(@{$gc_array_group_summary{$$ref_repClass_hash{$name}}{'len_down'}},$len_down) if $len_down=~/[\d\.]+/;
			#push(@{$gc_array_group_summary{$$ref_repClass_hash{$name}}{'len_updown'}},$len_up) if $len_up=~/[\d\.]+/;
			#push(@{$gc_array_group_summary{$$ref_repClass_hash{$name}}{'len_updown'}},$len_down) if $len_down=~/[\d\.]+/;
		}
		# if has element does not exist, may be names are not exactly same in hash key value pair. use pattern match to find class/group
		elsif(my$key=matches_to_key($name,$ref_repClass_hash)){
			#foreach my$key( keys %{$ref_repClass_hash}){
				#if($name=~/$key/){
					push(@{$gc_array_group_summary{$$ref_repClass_hash{$key}}{'gcelement'}},$gcelement) if $gcelement=~/[\d\.]+/;
					push(@{$gc_array_group_summary{$$ref_repClass_hash{$key}}{'gcup'}},$gcup) if $gcup=~/[\d\.]+/;
					push(@{$gc_array_group_summary{$$ref_repClass_hash{$key}}{'gcdown'}},$gcdown) if $gcdown=~/[\d\.]+/;
					push(@{$gc_array_group_summary{$$ref_repClass_hash{$key}}{'gcupdown'}},$gcup) if $gcup=~/[\d\.]+/;
					push(@{$gc_array_group_summary{$$ref_repClass_hash{$key}}{'gcupdown'}},$gcdown) if $gcdown=~/[\d\.]+/;

					#process length
					push(@{$gc_array_group_summary{$$ref_repClass_hash{$key}}{'len_element'}},$len_element) if $len_element=~/[\d\.]+/;
					#push(@{$gc_array_group_summary{$$ref_repClass_hash{$key}}{'len_up'}},$len_up) if $len_up=~/[\d\.]+/;
					#push(@{$gc_array_group_summary{$$ref_repClass_hash{$key}}{'len_down'}},$len_down) if $len_down=~/[\d\.]+/;
					#push(@{$gc_array_group_summary{$$ref_repClass_hash{$key}}{'len_updown'}},$len_up) if $len_up=~/[\d\.]+/;
					#push(@{$gc_array_group_summary{$$ref_repClass_hash{$key}}{'len_updown'}},$len_down) if $len_down=~/[\d\.]+/;
				#}
			#}
		}
		else{
				push(@{$gc_array_group_summary{'unclassified'}{'gcelement'}},$gcelement) if $gcelement=~/[\d\.]+/;
				push(@{$gc_array_group_summary{'unclassified'}{'gcup'}},$gcup) if $gcup=~/[\d\.]+/;
				push(@{$gc_array_group_summary{'unclassified'}{'gcdown'}},$gcdown) if $gcdown=~/[\d\.]+/;
				push(@{$gc_array_group_summary{'unclassified'}{'gcupdown'}},$gcup) if $gcup=~/[\d\.]+/;
				push(@{$gc_array_group_summary{'unclassified'}{'gcupdown'}},$gcdown) if $gcdown=~/[\d\.]+/;
				#process for length
				push(@{$gc_array_group_summary{'unclassified'}{'len_element'}},$len_element) if $len_element=~/[\d\.]+/;
				#push(@{$gc_array_group_summary{'unclassified'}{'len_up'}},$len_up) if $len_up=~/[\d\.]+/;
				#push(@{$gc_array_group_summary{'unclassified'}{'len_down'}},$len_down) if $len_down=~/[\d\.]+/;
				#push(@{$gc_array_group_summary{'unclassified'}{'len_updown'}},$len_up) if $len_up=~/[\d\.]+/;
				#push(@{$gc_array_group_summary{'unclassified'}{'len_updown'}},$len_down) if $len_down=~/[\d\.]+/;

		}
	}
	# Calculate statistics for the datapoints in each array and print the summary.
	print "\n\t->Calculating statistics for collected datapoints at element name level......\n";
	print {$file_handle_to_print}"\nRepeat_Name\t#Elements\tSum_of_lengths\t\%of_Total_seq\tGC_element_mean\tGC_element_stdev\t\#Up_seq\tGC_up_mean\tGC_up_standard_deviation\t\#DownSeq\tGC_down_mean\tGC_down_standard_deviation\t\#UpDownSeqs\tGC_updown_mean\tGC_updown_standard_deviation\n";

	foreach( keys %gc_array_summary){

		my$stat_element = Statistics::Descriptive::Full->new();
		$stat_element->add_data(@{$gc_array_summary{$_}{'gcelement'}});
		if($stat_element->count()!=0){$gc_array_summary{$_}{'gcelement_mean'}=sprintf("%0.2f",$stat_element->mean())} else {$gc_array_summary{$_}{'gcelement_mean'}='NoSeq'}
		if($stat_element->count()!=0){$gc_array_summary{$_}{'gcelement_stdev'}=sprintf("%0.2f",$stat_element->standard_deviation())}else {$gc_array_summary{$_}{'gcelement_stdev'}='NoSeq'}

		my$stat_up = Statistics::Descriptive::Full->new();
		$stat_up->add_data(@{$gc_array_summary{$_}{'gcup'}});
		if($stat_up->count()!=0){$gc_array_summary{$_}{'gcup_mean'}=sprintf("%0.2f",$stat_up->mean())} else {$gc_array_summary{$_}{'gcup_mean'}='NoSeq'}
		if($stat_up->count()!=0){$gc_array_summary{$_}{'gcup_stdev'}=sprintf("%0.2f",$stat_up->standard_deviation())}else {$gc_array_summary{$_}{'gcup_stdev'}='NoSeq'}

		my$stat_down = Statistics::Descriptive::Full->new();
		$stat_down->add_data(@{$gc_array_summary{$_}{'gcdown'}});
		if($stat_down->count()!=0){$gc_array_summary{$_}{'gcdown_mean'}=sprintf("%0.2f",$stat_down->mean())}else {$gc_array_summary{$_}{'gcdown_mean'}='NoSeq'}
		if($stat_down->count()!=0){$gc_array_summary{$_}{'gcdown_stdev'}=sprintf("%0.2f",$stat_down->standard_deviation())}else {$gc_array_summary{$_}{'gcdown_stdev'}='NoSeq'}

		my$stat_updown = Statistics::Descriptive::Full->new();
		$stat_updown->add_data(@{$gc_array_summary{$_}{'gcupdown'}});
		if($stat_updown->count()!=0){$gc_array_summary{$_}{'gcupdown_mean'}=sprintf("%0.2f",$stat_updown->mean())}else {$gc_array_summary{$_}{'gcupdown_mean'}='NoSeq'}
		if($stat_updown->count()!=0){$gc_array_summary{$_}{'gcupdown_stdev'}=sprintf("%0.2f",$stat_updown->standard_deviation())}else {$gc_array_summary{$_}{'gcupdown_stdev'}='NoSeq'}

		# process stat fro length of element and gc up-downs..
		my$stat_len_element = Statistics::Descriptive::Full->new();
		$stat_len_element->add_data(@{$gc_array_summary{$_}{'len_element'}});
		if($stat_len_element->count()!=0){$gc_array_summary{$_}{'len_element_sum'}=$stat_len_element->sum()} else {$gc_array_summary{$_}{'len_element_sum'}='NoSeq'}
		#if($stat_len_element->count()!=0){$gc_array_summary{$_}{'len_element_stdev'}=sprintf("%0.2f",$stat_len_element->standard_deviation())}else {$gc_array_summary{$_}{'len_element_stdev'}='NoSeq'}

		#my$stat_len_up = Statistics::Descriptive::Full->new();
		#$stat_len_up->add_data(@{$gc_array_summary{$_}{'len_up'}});
		#if($stat_len_up->count()!=0){$gc_array_summary{$_}{'len_up_mean'}=sprintf("%0.2f",$stat_len_up->mean())} else {$gc_array_summary{$_}{'len_up_mean'}='NoSeq'}
		#if($stat_len_up->count()!=0){$gc_array_summary{$_}{'len_up_stdev'}=sprintf("%0.2f",$stat_len_up->standard_deviation())}else {$gc_array_summary{$_}{'len_up_stdev'}='NoSeq'}
		#
		#my$stat_len_down = Statistics::Descriptive::Full->new();
		#$stat_len_down->add_data(@{$gc_array_summary{$_}{'len_down'}});
		#if($stat_len_down->count()!=0){$gc_array_summary{$_}{'len_down_mean'}=sprintf("%0.2f",$stat_len_down->mean())}else {$gc_array_summary{$_}{'len_down_mean'}='NoSeq'}
		#if($stat_len_down->count()!=0){$gc_array_summary{$_}{'len_down_stdev'}=sprintf("%0.2f",$stat_len_down->standard_deviation())}else {$gc_array_summary{$_}{'len_down_stdev'}='NoSeq'}
		#
		#my$stat_len_updown = Statistics::Descriptive::Full->new();
		#$stat_len_updown->add_data(@{$gc_array_summary{$_}{'len_updown'}});
		#if($stat_len_updown->count()!=0){$gc_array_summary{$_}{'len_updown_mean'}=sprintf("%0.2f",$stat_len_updown->mean())}else {$gc_array_summary{$_}{'len_updown_mean'}='NoSeq'}
		#if($stat_len_updown->count()!=0){$gc_array_summary{$_}{'len_updown_stdev'}=sprintf("%0.2f",$stat_len_updown->standard_deviation())}else {$gc_array_summary{$_}{'len_updown_stdev'}='NoSeq'}




		print {$file_handle_to_print}"$_\t",
		$stat_len_element->count(),"\t",
		$gc_array_summary{$_}{'len_element_sum'},"\t",
		sprintf("%0.3f",($gc_array_summary{$_}{'len_element_sum'}*100/$total_len)),"\t",
		$gc_array_summary{$_}{'gcelement_mean'},"\t",
		$gc_array_summary{$_}{'gcelement_stdev'},"\t",
		$stat_up->count(),"\t",
		$gc_array_summary{$_}{'gcup_mean'},"\t",
		$gc_array_summary{$_}{'gcup_stdev'},"\t",
		$stat_down->count(),"\t",
		$gc_array_summary{$_}{'gcdown_mean'},"\t",
		$gc_array_summary{$_}{'gcdown_stdev'},"\t",
		$stat_updown->count(),"\t",
		$gc_array_summary{$_}{'gcupdown_mean'},"\t",
		$gc_array_summary{$_}{'gcupdown_stdev'},"\n";
		# for length
		#$stat_len_element->count(),"\t",
		#$gc_array_summary{$_}{'len_element_mean'},"\t",
		#$gc_array_summary{$_}{'len_element_stdev'},"\t",
		#$stat_len_up->count(),"\t",
		#$gc_array_summary{$_}{'len_up_mean'},"\t",
		#$gc_array_summary{$_}{'len_up_stdev'},"\t",
		#$stat_len_down->count(),"\t",
		#$gc_array_summary{$_}{'len_down_mean'},"\t",
		#$gc_array_summary{$_}{'len_down_stdev'},"\t",
		#$stat_len_updown->count(),"\t",
		#$gc_array_summary{$_}{'len_updown_mean'},"\t",
		#$gc_array_summary{$_}{'len_updown_stdev'},"\n";

	}
	# Calculate group statistics for the datapoints in each array and print the group summary.
	print "\n\t->Calculating statistics for collected datapoints at Order/Family level......\n";
	print {$file_handle_to_print}"\n\nRepeat_Group_Summary\n\nRepeat_Name\t#Elements\tSum_of_lengths\t\%of_Total_seq\tGC_element_mean\tGC_element_stdev\t\#Up_seq\tGC_up_mean\tGC_up_standard_deviation\t\#DownSeq\tGC_down_mean\tGC_down_standard_deviation\t\#UpDownSeqs\tGC_updown_mean\tGC_updown_standard_deviation\n";
	foreach( keys %gc_array_group_summary){

		my$stat_element = Statistics::Descriptive::Full->new();
		$stat_element->add_data(@{$gc_array_group_summary{$_}{'gcelement'}});
		if($stat_element->count()!=0){$gc_array_group_summary{$_}{'gcelement_mean'}=sprintf("%0.2f",$stat_element->mean())} else {$gc_array_group_summary{$_}{'gcelement_mean'}='NoSeq'}
		if($stat_element->count()!=0){$gc_array_group_summary{$_}{'gcelement_stdev'}=sprintf("%0.2f",$stat_element->standard_deviation())}else {$gc_array_group_summary{$_}{'gcelement_stdev'}='NoSeq'}

		my$stat_up = Statistics::Descriptive::Full->new();
		$stat_up->add_data(@{$gc_array_group_summary{$_}{'gcup'}});
		if($stat_up->count()!=0){$gc_array_group_summary{$_}{'gcup_mean'}=sprintf("%0.2f",$stat_up->mean())} else {$gc_array_group_summary{$_}{'gcup_mean'}='NoSeq'}
		if($stat_up->count()!=0){$gc_array_group_summary{$_}{'gcup_stdev'}=sprintf("%0.2f",$stat_up->standard_deviation())}else {$gc_array_group_summary{$_}{'gcup_stdev'}='NoSeq'}

		my$stat_down = Statistics::Descriptive::Full->new();
		$stat_down->add_data(@{$gc_array_group_summary{$_}{'gcdown'}});
		if($stat_down->count()!=0){$gc_array_group_summary{$_}{'gcdown_mean'}=sprintf("%0.2f",$stat_down->mean())}else {$gc_array_group_summary{$_}{'gcdown_mean'}='NoSeq'}
		if($stat_down->count()!=0){$gc_array_group_summary{$_}{'gcdown_stdev'}=sprintf("%0.2f",$stat_down->standard_deviation())}else {$gc_array_group_summary{$_}{'gcdown_stdev'}='NoSeq'}

		my$stat_updown = Statistics::Descriptive::Full->new();
		$stat_updown->add_data(@{$gc_array_group_summary{$_}{'gcupdown'}});
		if($stat_updown->count()!=0){$gc_array_group_summary{$_}{'gcupdown_mean'}=sprintf("%0.2f",$stat_updown->mean())}else {$gc_array_group_summary{$_}{'gcupdown_mean'}='NoSeq'}
		if($stat_updown->count()!=0){$gc_array_group_summary{$_}{'gcupdown_stdev'}=sprintf("%0.2f",$stat_updown->standard_deviation())}else {$gc_array_group_summary{$_}{'gcupdown_stdev'}='NoSeq'}

		# groutp stat for length
		my$stat_len_element = Statistics::Descriptive::Full->new();
		$stat_len_element->add_data(@{$gc_array_group_summary{$_}{'len_element'}});
		if($stat_len_element->count()!=0){$gc_array_group_summary{$_}{'len_element_sum'}=$stat_len_element->sum()} else {$gc_array_group_summary{$_}{'len_element_sum'}='NoSeq'}
		#if($stat_len_element->count()!=0){$gc_array_group_summary{$_}{'len_element_stdev'}=sprintf("%0.2f",$stat_len_element->standard_deviation())}else {$gc_array_group_summary{$_}{'len_element_stdev'}='NoSeq'}

		#my$stat_len_up = Statistics::Descriptive::Full->new();
		#$stat_len_up->add_data(@{$gc_array_group_summary{$_}{'len_up'}});
		#if($stat_len_up->count()!=0){$gc_array_group_summary{$_}{'len_up_mean'}=sprintf("%0.2f",$stat_len_up->mean())} else {$gc_array_group_summary{$_}{'len_up_mean'}='NoSeq'}
		#if($stat_len_up->count()!=0){$gc_array_group_summary{$_}{'len_up_stdev'}=sprintf("%0.2f",$stat_len_up->standard_deviation())}else {$gc_array_group_summary{$_}{'len_up_stdev'}='NoSeq'}
		#
		#my$stat_len_down = Statistics::Descriptive::Full->new();
		#$stat_len_down->add_data(@{$gc_array_group_summary{$_}{'len_down'}});
		#if($stat_len_down->count()!=0){$gc_array_group_summary{$_}{'len_down_mean'}=sprintf("%0.2f",$stat_len_down->mean())}else {$gc_array_group_summary{$_}{'len_down_mean'}='NoSeq'}
		#if($stat_len_down->count()!=0){$gc_array_group_summary{$_}{'len_down_stdev'}=sprintf("%0.2f",$stat_len_down->standard_deviation())}else {$gc_array_group_summary{$_}{'len_down_stdev'}='NoSeq'}
		#
		#my$stat_len_updown = Statistics::Descriptive::Full->new();
		#$stat_len_updown->add_data(@{$gc_array_group_summary{$_}{'len_updown'}});
		#if($stat_len_updown->count()!=0){$gc_array_group_summary{$_}{'len_updown_mean'}=sprintf("%0.2f",$stat_len_updown->mean())}else {$gc_array_group_summary{$_}{'len_updown_mean'}='NoSeq'}
		#if($stat_len_updown->count()!=0){$gc_array_group_summary{$_}{'len_updown_stdev'}=sprintf("%0.2f",$stat_len_updown->standard_deviation())}else {$gc_array_group_summary{$_}{'len_updown_stdev'}='NoSeq'}

		# print result

		print {$file_handle_to_print}"$_\t",
		$stat_len_element->count(),"\t",
		$gc_array_group_summary{$_}{'len_element_sum'},"\t",
		sprintf("%0.3f",($gc_array_group_summary{$_}{'len_element_sum'}*100/$total_len)),"\t",
		$gc_array_group_summary{$_}{'gcelement_mean'},"\t",
		$gc_array_group_summary{$_}{'gcelement_stdev'},"\t",
		$stat_up->count(),"\t",
		$gc_array_group_summary{$_}{'gcup_mean'},"\t",
		$gc_array_group_summary{$_}{'gcup_stdev'},"\t",
		$stat_down->count(),"\t",
		$gc_array_group_summary{$_}{'gcdown_mean'},"\t",
		$gc_array_group_summary{$_}{'gcdown_stdev'},"\t",
		$stat_updown->count(),"\t",
		$gc_array_group_summary{$_}{'gcupdown_mean'},"\t",
		$gc_array_group_summary{$_}{'gcupdown_stdev'},"\n",
		# for length
		#$stat_len_element->count(),"\t",
		#$gc_array_group_summary{$_}{'len_element_mean'},"\t",
		#$gc_array_group_summary{$_}{'len_element_stdev'},"\t",
		#$stat_len_up->count(),"\t",
		#$gc_array_group_summary{$_}{'len_up_mean'},"\t",
		#$gc_array_group_summary{$_}{'len_up_stdev'},"\t",
		#$stat_len_down->count(),"\t",
		#$gc_array_group_summary{$_}{'len_down_mean'},"\t",
		#$gc_array_group_summary{$_}{'len_down_stdev'},"\t",
		#$stat_len_updown->count(),"\t",
		#$gc_array_group_summary{$_}{'len_updown_mean'},"\t",
		#$gc_array_group_summary{$_}{'len_updown_stdev'},"\n";

	}







}






##############################################################################################
# return minimum of all the numbers
##############################################################################################
sub min{
  @_= sort { $a <=> $b }@_;
   return $_[0];

}
##############################################################################################
# return maximum of all the numbers
##############################################################################################
sub max{

 @_ = sort { $a <=> $b }@_;
   return $_[-1];


}

##############################################################################################
# return GC content of a sequence.Non ATGC are not included in calculaions
##############################################################################################
sub GC_content{
	my$sequence=shift;
	return 0 if !$sequence;
	my$GCnumber=$sequence=~tr/GCgc/GCgc/;
	my$length_seq=$sequence=~tr/ATGCatgc/ATGCatgc/;
	#print "\nSequence length is 0\nSequence:$sequence\n" if $length_seq==0;
	my$GCcontent=$GCnumber/$length_seq if $length_seq!=0;
	return(sprintf("%0.2f",$GCcontent))if $length_seq!=0;
	return(0) if $length_seq==0;
}

##############################################################################################
# return GC content of sequence in percent. Non ATGC are not included in calculaions
##############################################################################################
sub GC_content_percent{
	my$sequence=shift;
	return 0 if !$sequence;
	my$GCnumber=$sequence=~tr/GCgc/GCgc/;
	my$length_seq=$sequence=~tr/ATGCatgc/ATGCatgc/;

	#print "\nSequence length is 0\nSequence:$sequence\n" if $length_seq==0;
	my$GCcontent=$GCnumber/$length_seq if $length_seq!=0;
	return(sprintf("%0.2f",$GCcontent*100))if $length_seq!=0;
	return(0) if $length_seq==0;
}

##############################################################################################
# return length (only ATGC)
##############################################################################################
sub num_ATGC{
	my$sequence=shift;
	my$length_seq=$sequence=~tr/ATGCatgc/ATGCatgc/;
	return($length_seq);

}

##############################################################################################
# return key if name matches to any key of the hash else return false
##############################################################################################

sub matches_to_key{
	my$name=shift;
	my$ref_repClass_hash=shift;
	my$return_value;
	my$matchfound='false';
	foreach(keys %{$ref_repClass_hash}){
		$return_value=$_ if $name=~/$_/;
		$matchfound='true';
	}
	if($matchfound){return($return_value)}
	else{return('false')}
}

# progress bar
#
#------------------------------------------------------------------------------------
sub show_progress{
		my$total_lines=shift;
		my$current_count=shift;
		my$next_value=shift;
		$current_count++;
		if(($current_count*100/$total_lines)>$next_value){
			syswrite STDOUT,"=";$next_value+=1;
			syswrite STDOUT,"$next_value\%" if $next_value%10==0;
		}
		return($total_lines,$current_count,$next_value)
}
#################################################################################################
# Clean Description for unnecessary comments e.g Target" or Target= etc.
#################################################################################################
sub clean_comment {

	my$text=shift; # reference to the text to be cleaned
	$$text=~s/\s+[\s\d]+$//g;
	$$text=~s/Target=//g;
	$$text=~s/"//g;
	$$text=~s/Target Motif://g;
}
