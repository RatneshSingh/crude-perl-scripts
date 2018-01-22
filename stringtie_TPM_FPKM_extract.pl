#! /usr/bin/perl
use warnings;
use strict;
use Cwd;
use Data::Dumper;
no autovivification;
our (%tpm,%fpm, %cov, @flist,$tpof,$fpof,$cvof,%id_table,%trcs_as_table,%samp_list,%as_ex,%hist_data);

my$data_dir=$ARGV[0];
my$pat=$ARGV[1];
my $min_exp =$ARGV[2];
my $max_exp = $ARGV[3];

$max_exp=$max_exp?$max_exp:10000000000000000000;
$min_exp=$min_exp?$min_exp:0;
$data_dir=$data_dir?$data_dir:Cwd::cwd();
$pat=$pat?$pat:'out';
opendir DIR, $data_dir;
print "\n Looking for file with pattern $pat in $data_dir\n";


if($pat=~/out/i){
    $tpof="Turf_gene_expression_Diamond_TPM.table"   ;
    $fpof="Turf_gene_expression_Diamond_FPKM.table"   ;
    $cvof="Turf_gene_expression_Diamond_COV.table"   ;
}elsif($pat=~/gtf/i){
    $tpof="Turf_transcript_expression_Diamond_TPM.table"   ;
    $fpof="Turf_transcript_expression_Diamond_FPKM.table"   ;
    $cvof="Turf_transcript_expression_Diamond_COV.table"   ;

}





while(my$file=readdir(DIR)){
    next if $file eq '.';
    next if $file eq  '..';
    print "\n\nAnalysing file $file";
    print " ...no match with $pat  ...Skipping\n" if $file !~ /$pat/;
    next if $file !~ /$pat/;
    print " ...Matched restricted name $1 .... Skipping\n" if $file =~ /(Cavalier|Meyer|4790|4895|4905|4935|4945|501)/i;
    next if $file =~ /Cavalier|Meyer|4790|4895|4905|4935|4945|501/i;
    my$dread=0;
    print "\n*******Collecting TPM, FPKM and COV values from: $file\n";
    my$samp_id=$file;
    $samp_id=~s/HISAT2_map_|_L000__R1R2_P.*|.gtf|.out//g;
    push(@flist,$samp_id);
    $samp_list{$samp_id}=$file;
    open CNT,$file;
    if($pat=~/out/){
        while(<CNT>){
            next if /^Gene.ID/;
            next if /^\s*$/;
            my@line=split /\t/;
            next if $line[0] =~ /^STRG/;

            chomp foreach @line;
            $cov{$line[0]}{$samp_id}=$line[6];
            $fpm{$line[0]}{$samp_id}=$line[7];
            $tpm{$line[0]}{$samp_id}=$line[8];

            $dread++;
        }
    }
    elsif($pat =~/gtf/){
        while(<CNT>){
            next if /^#/;
            next if /^\s*$/;
            my@line=split /\t/;
            chomp foreach @line;
            next if $line[2] !~ /transcript/i;
            my@gtel=split(/;/,$line[8]);
            my($gene_id,$transcript_id,$reference_id,$ref_gene_id,$fpkm,$cov,$tpm);
            foreach my $gtin(@gtel){
                $gene_id=$1 if($gtin=~/gene_id\s*\"([^\"]+)\"/);
                $transcript_id=$1 if($gtin=~/transcript_id\s*\"([^\"]+)\"/);
                $reference_id=$1 if($gtin=~/reference_id\s*\"([^\"]+)\"/);
                $ref_gene_id=$1 if($gtin=~/ref_gene_id\s*\"([^\"]+)\"/);
                $fpkm=$1 if($gtin=~/FPKM\s*\"([^\"]+)\"/);
                $cov=$1 if($gtin=~/cov\s*\"([^\"]+)\"/);
                $tpm=$1 if($gtin=~/TPM\s*\"([^\"]+)\"/);





            }
            next if !$reference_id;

            $fpm{$reference_id}{$samp_id}=$fpkm;
            $cov{$reference_id}{$samp_id}=$cov;
            $tpm{$reference_id}{$samp_id}=$tpm;
            ### create dictionary to assign samples gene/transcript IDs to merged IDs.
            $id_table{$gene_id}=$ref_gene_id;
            $id_table{$transcript_id}=$reference_id;



            $dread++;


        }
    }
    print "*******Read $dread genes/transcripts info for Sample: $samp_id ($file)\n";
}



print "\n\nPrinting FPKM, TPM and COV values in external files";
open TPM, ">",$tpof  ;#if($pat =~/out/);
open FPM, ">",$fpof;
open COV, ">",$cvof;#if($pat =~/gtf/);
print TPM "Gene.ID\t".join("\t",@flist) ;#if($pat =~/out/);
print FPM "Gene.ID\t".join("\t",@flist);
print COV "Transcript.ID\t".join("\t",@flist) ;# if($pat =~/gtf/);
foreach my$gene(sort keys %fpm){
    print TPM "\n$gene" ;# if($pat =~/out/);;
    print FPM "\n$gene";
    print COV "\n$gene" ;#if($pat =~/gtf/);
    foreach my$samp_id(sort @flist){
        if(exists($fpm{$gene}{$samp_id})){
		   print TPM "\t$tpm{$gene}{$samp_id}" ;# if($pat =~/out/);
           print FPM "\t$fpm{$gene}{$samp_id}";
           print COV "\t$cov{$gene}{$samp_id}" ;#if($pat =~/gtf/);
	}else{
           print TPM "\t0"  ;#if($pat =~/out/);
           print FPM "\t0";
           print COV "\t0"   ;#if($pat =~/gtf/);
	}
    }
}

print "\nResults are saved in :$tpof, \t $fpof , and \n $cvof files\n";

my $astad="/home/ratnesh.singh/Turf/Turf_RNASeq_alternatesplicing/astalavista";
opendir ASTA, $astad;
while(my$astaf=readdir(ASTA)){
    next if $astaf !~ /ASTALVISTA_ASI.Out.gtf$/;
    my$samp_id=$astaf;
    $samp_id=~s/HISAT2_map_|_L000__R1R2_P.*//g;
    next if ! exists $samp_list{$samp_id};
    open ASTAF,"<","$astad/$astaf" or die "Unable to open file: $astaf\n";
    print "Processing astalvista data for sample $samp_id ($astaf)\n";
    while (my$astal = <ASTAF>){
        next if $astal=~/^\s*$|^#/;
        
        my($transcript_id,$structure)=process_asta($astal);
        
        my($trcs1,$trcs2)=split(/,/,$transcript_id);                        ## "HISAT2_map_1_11_TAGC.23609.4    ,    HISAT2_map_1_11_TAGC.23609.1/HISAT2_map_1_11_TAGC.23609.2/HISAT2_map_1_11_TAGC.23609.3"; 
        my($strs1,$strs2)=split(/,/,$structure);                            ## "0                               ,    1^2-"
        
        #collect_strs(\%trcs_as_table,$samp_id,$trcs1,$strs1);
        #collect_strs(\%trcs_as_table,$samp_id,$trcs2,$strs2);
        assign_as(\%trcs_as_table,$samp_id,$transcript_id,$structure);
        
    }
    close(ASTAF);
}
my$exp="TPM";
my$asof="Turf_Diamond_altsplice_expression_$exp.table";
open(ASOF,">",$asof) or die "Unable to open $asof\n";
print ASOF join("\t","Transcripts",map{$_."_AS",$_."_$exp"}@flist);

foreach my$mgenes(sort keys $trcs_as_table{merged}){
    print  ASOF "\n$mgenes";
    foreach my$samp_id(sort @flist){
        my$strs=$trcs_as_table{merged}{$mgenes}{$samp_id}?unique_list($trcs_as_table{merged}{$mgenes}{$samp_id}):"NA";
        my$tpm=$tpm{$mgenes}{$samp_id}?$tpm{$mgenes}{$samp_id}:0;
        print ASOF "\t$strs\t$tpm";
        
        ## collect expression values for each category of AS
        if ($tpm >= $min_exp && $tpm < $max_exp) {
            #push(@{$as_ex{$strs}},$tpm) if $strs=~/\,/;
            foreach my$lstrs(split(/,/,$strs)){
                push(@{$as_ex{$samp_id}{$lstrs}},$tpm);
                push(@{$hist_data{$lstrs}},$tpm);
            }
        }
    }
}

foreach my$strs(keys %hist_data){
        draw_histogram(\@{$hist_data{$strs}},$strs);
}



foreach my$samp_id(keys %as_ex){
    my%thash=%{$as_ex{$samp_id}};
    draw_boxplot(\%thash,"$samp_id.boxplot");
}
#draw_boxplot(\%as_ex,"Diamond_AS_EX");

close(ASOF);











###########################################################################################################33
sub process_asta{
    
    my$astal=shift;
    
    my@astel=split(/\t+/,$astal);
        my@astcon=split(/;/,$astel[8]);   ####  transcript_id "HISAT2_map_1_11_TAGC.23609.4,HISAT2_map_1_11_TAGC.23609.1/HISAT2_map_1_11_TAGC.23609.2/HISAT2_map_1_11_TAGC.23609.3"; gene_id "Super_scaffold_124:573468-578384W"; flanks "576986-,577237^"; structure "0,1^2-"; splice_chain ",577051^577143-"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2";
        my($gene_id,$transcript_id,$flanks,$structure,$splice_chain,$degree,$dimension);
        foreach my$gtin(@astcon){
            $gene_id=$1 if($gtin=~/gene_id\s*\"([^\"]+)\"/);                ##  gene_id "Super_scaffold_124:573468-578384W"; 
            $transcript_id=$1 if($gtin=~/transcript_id\s*\"([^\"]+)\"/);    ##  transcript_id "HISAT2_map_1_11_TAGC.23609.4,HISAT2_map_1_11_TAGC.23609.1/HISAT2_map_1_11_TAGC.23609.2/HISAT2_map_1_11_TAGC.23609.3"; 
            $flanks=$1 if($gtin=~/flanks\s*\"([^\"]+)\"/);                  ##  flanks "576986-,577237^"; 
            $structure=$1 if($gtin=~/structure\s*\"([^\"]+)\"/);            ##  structure "0,1^2-";
            $splice_chain=$1 if($gtin=~/splice_chain\s*\"([^\"]+)\"/);      ##  splice_chain ",577051^577143-";
            $degree=$1 if($gtin=~/degree\s*\"([^\"]+)\"/);                  ##  degree "4";
            $dimension=$1 if($gtin=~/dimension\s*\"([^\"]+)\"/);            ##  dimension "2_2";
        }
    
    return($transcript_id,$structure)
    
}
sub collect_strs{
    my$refHash=shift;
    my$samp_id=shift;
    my$trcs1=shift;
    my$strs1=shift;
    my$structure=shift;
    foreach my$trcs(split(/\//,$trcs1)){
            if(exists $$refHash{ind}{$trcs}{$samp_id}){
                $$refHash{ind}{$trcs}{$samp_id}=join(",",$$refHash{ind}{$trcs}{$samp_id},$strs1);
            }else{
                $$refHash{ind}{$trcs}{$samp_id}=$strs1 ;
            }
            if(exists $id_table{$trcs} ){
                #print "\nSample TID:$trcs\t RefTID:$id_table{$trcs}";
                if(exists $$refHash{merged}{$id_table{$trcs}}{$samp_id}){
                    $$refHash{merged}{$id_table{$trcs}}{$samp_id}=join(",",$$refHash{merged}{$id_table{$trcs}}{$samp_id},$strs1);
                }else{
                    $$refHash{merged}{$id_table{$trcs}}{$samp_id}=$strs1;
                }
            }
        }
}


sub assign_as{
    my$refHash=shift;
    my$samp_id=shift;
    my$trcsf   =shift;
    my$strs   =shift;
    my($trcs1,$trcs2)=split(/\,/,$trcsf);
    my$ascode=code_to_as($strs);
    #print "\nSample TID:$trcsf\tStrs:$strs\tAs_code:$ascode";
    foreach my $trcs (split(/\//,$trcs2)){
            if(exists $$refHash{ind}{$trcs}{$samp_id}){
                $$refHash{ind}{$trcs}{$samp_id}=join(",",$$refHash{ind}{$trcs}{$samp_id},$ascode);
            }else{
                $$refHash{ind}{$trcs}{$samp_id}=$ascode ;
            }
            if(exists $id_table{$trcs} ){
                #print "\nSample TID:$trcs\t RefTID:$id_table{$trcs}\t$ascode";
                if(exists $$refHash{merged}{$id_table{$trcs}}{$samp_id}){
                    $$refHash{merged}{$id_table{$trcs}}{$samp_id}=join(",",$$refHash{merged}{$id_table{$trcs}}{$samp_id},$ascode);
                }else{
                    $$refHash{merged}{$id_table{$trcs}}{$samp_id}=$ascode;
                }
            }
    }
    
}


sub unique_list{
    my$list=shift;
    my%thash=map{$_=>1} sort split(/,/,$list);   
    return join(",",keys %thash)
    
    
}


sub code_to_as{
    my$code=shift;
    my$astype='other';
    if ($code =~ /^0\,1\-2\^/){$astype='skipped_exon'}
    elsif ($code =~ /^0\,1\^2\-/){$astype='retained_introns'}
    elsif ($code =~ /^1\^\,2\^/){$astype='alt_donors'}
    elsif ($code =~ /^1\-2\^\,3\-4\^/){ $astype='mutually_exclusive_exons'}
    elsif ($code =~ /^1\-\,2\-/){$astype='alt_acceptors'}
    elsif ($code =~ /^1\[/){$astype='alt_TSS'}
    else {$astype='other'}

    return($astype)
}



sub draw_histogram{
    my$refData=shift;
    my$title=shift;
    
    use GD::Graph::histogram;
    my $hgraph= new GD::Graph::histogram(800,1200);
    
    $hgraph->set(title=>"$title",histogram_bins=>100, histogram_type=>'percentage');
    $hgraph->set_x_label_font( ['verdana', 'arial', 'gdMediumBoldFont'], 28);
    $hgraph->set_y_label_font( ['verdana', 'arial', 'gdMediumBoldFont'], 28);
    my$gd= $hgraph->plot($refData);
    
    open(IMG,'>',"$title.png") or die $hgraph->error;
    binmode IMG;
    print IMG $gd->png;
    close IMG;
}

sub draw_boxplot{
    use GD::Graph::boxplot;
    
    my$refData=shift;
    my$title=shift;
    my $hgraph= new GD::Graph::boxplot(1200,1200);
    
    my(@data);
    $data[0]=[];
    $data[1]=[];
    for my$dlab(keys %{$refData}){
        my$tdata=$$refData{$dlab};
        push($data[0],$dlab);
        push($data[1],$tdata);
    }
    
    
    $hgraph->set(title=>"$title",y_max_value=>200,y_min_value=>-10,x_labels_vertical=>1);
    $hgraph->set_x_label_font( ['verdana', 'arial', 'gdMediumBoldFont'], 28);
    $hgraph->set_y_label_font( ['verdana', 'arial', 'gdMediumBoldFont'], 28);
    #print Dumper(@data);
    
    
    my$gd= $hgraph->plot(\@data);
    
    open(IMG,'>',"$title.boxplot.png") or die $hgraph->error;
    binmode IMG;
    print IMG $gd->png;
    close IMG;
    
}