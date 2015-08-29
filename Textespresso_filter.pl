#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
 use Term::ANSIColor qw(:constants);
  $Term::ANSIColor::AUTORESET = 1;
our (%citation,%get_title,@all_genes,$text_file,$and,$or,@patterns,%print,$search_titles,$nocolor);
my $gene_len=3;
my$options=GetOptions(
  "text|t=s"=>\$text_file,
  "and"=>\$and,
  "or"=>\$or,
  "query|q=s{,}"=>\@patterns,
  "title"=>\$print{'title'},
  "authors"=>\$print{'author'},
  "citation"=>\$print{'citation'},
  "genes"=>\$print{'genes'},
  "abstract"=>\$print{'abstract'},
  "genes_all|genes_all"=>\$print{'all_genes'},
  "printall|print_all|all"=>\$print{'print_all'},
  "len_gene_name|lgn=i"=>\$gene_len,
  "search_titles|titlesearch"=>\$search_titles,
  "no_color|nc|nocolor"=>\$nocolor
  );


my$usage="

-text|t=s    text_file containing search output data in endnote format
-and    Show the entries where all the query words occurs atleast once.
-or    Show the entries where atleast one of the queries occurs once
-query|q    key words to search for. Use mulitple words seperated by space.
-title    print 'title',
-authors    print 'author'
-citation    print 'citation'
-genes    print 'genes',
-abstract    print 'abstract'
-genes_all|genes_all    print 'all_genes' at the end
-printall|print_all|all    print everything
-len_gene_name|lgn=i    minimum number of characters in the gene abbreviation
-search_titles|titlesearch    only search in titles


";




open TEXT, "$text_file" or die "\nError--->Cannot open file: $text_file\n$usage";




# define "and" or "or" not both or none.
if($and && $or){print "\n\"and\" and \"or\" cannot be at the same time. Returning to default \"or\"\n"; undef $and;}
elsif(!$and && !$or){$or="true"}

if($print{'print_all'}){
  $print{'title'}=1;
  $print{'author'}=1;
  $print{'citation'}=1;
  $print{'genes'}=1;
  $print{'abstract'}=1;
  $print{'all_genes'}=1;
  }

print "Queries used:", join("\nQuery:",@patterns);
print "\nDocuments will be selected only if all the queries match [AND]\n" if $and;
print "\nDocuments will be selected if any of the queries match [OR]\n" if $or;
push(@patterns," ") if @patterns ==0 ;
$/="\n\%D";
while (<TEXT>){

    chomp;
    my $text_chunk=$_;

    #Collect tilte, abstract, reference to trace back the reference.
    my $title=$1 if $text_chunk=~/\%T\s+([^\n]+)/g;
    my $abstract=$1 if $text_chunk=~/\%X\s+([\w\W]+)/g;
    #$abstract=~s/\%D\s*//g;

    ## collect gene names and remove common words from such as DNA PCR etc..
    my@gene_names1;
    push(@gene_names1,$1) while($title=~/\b([A-Z]{0,1}[a-z]{0,1}[A-Z\-]{2,}[\dA-Z\s]*)\b/g);
    push(@gene_names1,$1) while($abstract=~/\b([A-Z]{0,1}[a-z]{0,1}[A-Z\-]{2,}[\dA-Z\s]*)\b/g);
    # make the gene names uniq by assigning them as keys of a hash and save in array
    my@gene_names=keys %{{ map{$_=>1}@gene_names1}};

    # remove common names and short names which were picked as gene.
    for(my$i=0;$i<scalar@gene_names;$i++){
      $gene_names[$i]=~s/^\s+//g;
      $gene_names[$i]=~s/\s+$//g;
      $gene_names[$i]=~s/\s+/ /g;
      $gene_names[$i]=~s/^[-]+//g;
      $gene_names[$i]=~s/[-]+$//g;
      my$tgnm=$gene_names[$i];
      $tgnm=~s/\s+//g;
      $tgnm=~s/[\b]+//g;
      if(
           length $tgnm < $gene_len
        || length $tgnm > 40
        || $gene_names[$i]=~/^A\s+/g
        || $tgnm eq 'DNA'
        || $tgnm eq 'PCR'
        || $tgnm eq 'ABA'
        || $tgnm eq 'ATP'
        || $tgnm eq 'RNA'
        || $tgnm eq 'siRNA'
        || $tgnm eq 'rRNA'
        || $tgnm eq 'rDNA'
        || $tgnm eq 'miRNA'
        || $tgnm eq 'mRNA'
        || $tgnm eq 'tRNA'
        || $tgnm eq 'mRNA'
        || $tgnm eq 'CO2'
        || $tgnm eq 'qPCR'
        || $tgnm eq 'qRT-PCR'
        || $tgnm !~ /[A-Z]{2,}/
        || $tgnm =~/^[ATGC]{4,}$/
        ){my$temp=splice(@gene_names,$i,1);$i--;}
    }







    #collect citation
    my $J = $text_chunk=~/\%J([^\n]+)/ ? $1 :"Not Found";
    my $V = $text_chunk=~/\%V([^\n]+)/ ? $1 :"Not Found";
    my $P = $text_chunk=~/\%P([^\n]+)/ ? $1 :"Not Found";
    my $Y;
    if ($text_chunk=~/\%D\s*([\d]+)\s*\n/){$Y = $1;}
    elsif($text_chunk=~/^\s*([\d]+)\s*\n/){$Y = $1 }
    else{$Y=0000}
    my$citation="$J($Y),$V,$P.";
    $citation=~s/^\s+//;
    $citation=~s/\s+$//;
    $citation=~s/\s+/ /;

    # collect authors name
    print "" if $text_chunk =~ /\%A([^\n]+)/g;
    my@authors;
    push(@authors,$1) while $text_chunk =~ /\%A([^\n]+)/g;
    my $author=join(", ",@authors);

      my$pick_this;
      if($and ){
        $pick_this="true";
        foreach my$query(@patterns){
          $pick_this="false" if $text_chunk!~/$query/i && !$search_titles;
          $pick_this="false" if $title!~/$query/i && $search_titles;
          last if $pick_this eq "false";
        }
      }
      elsif($or){
        $pick_this="false";
        foreach my$query(@patterns){
          $pick_this="true" if $text_chunk=~/$query/i && !$search_titles;
          $pick_this="true" if $title=~/$query/i && $search_titles;
        }
      }




    if($pick_this eq "true"){
      $citation{$title}{'title'}=$title;
      $citation{$title}{'author'}=$author;
      $citation{$title}{'citation'}=$citation;
      $citation{$title}{'abstract'}=$abstract;
      $citation{$title}{'year'}=$Y;
      push(@all_genes,@gene_names);
      $citation{$title}{'genes'}=join(",",@gene_names);
    }
    else{next}

}

my $num_selected=0;
foreach my$sel_title(sort { $citation{$a}{'year'} <=> $citation{$b}{'year'}} keys %citation){
  $num_selected++;

  if(!$nocolor){
      print "\n\n" if ($print{'title'} || $print{'author'} || $print{'citation'} ||$print{'genes'}||$print{'abstract'});
      print BOLD RED "\nTitle $num_selected: ",$sel_title if $print{'title'};
      print RED "\nAuthors: ",$citation{$sel_title}{'author'} if $print{'author'};
      print BOLD MAGENTA "\nCitation: ",$citation{$sel_title}{'citation'} if $print{'citation'};
      print BOLD BLUE "\nGenes: ",$citation{$sel_title}{'genes'} if $print{'genes'};
      print "\nAbstract: ",$citation{$sel_title}{'abstract'} if $print{'abstract'};
    }
  elsif($nocolor){
    print "\n\n" if ($print{'title'} || $print{'author'} || $print{'citation'} ||$print{'genes'}||$print{'abstract'});
    print "\nTitle $num_selected: ",$sel_title if $print{'title'};
    print "\nAuthors: ",$citation{$sel_title}{'author'} if $print{'author'};
    print "\nCitation: ",$citation{$sel_title}{'citation'} if $print{'citation'};
    print "\nGenes: ",$citation{$sel_title}{'genes'} if $print{'genes'};
    print "\nAbstract: ",$citation{$sel_title}{'abstract'} if $print{'abstract'};
  }


}

print BOLD GREEN ON_RED "\n\nTotal number of documents passed the filter:$num_selected\n" if !$nocolor;
print "\n\nTotal number of documents passed the filter:$num_selected\n" if $nocolor;
my@all_gene_names=keys %{{ map{$_=>1}@all_genes}};
print "\n\nGenes in selected Abstracts:";

print BOLD BLUE join(", ", sort @all_gene_names) if $print{'all_genes'} && !$nocolor;
print BOLD GREEN ON_RED "\n\nTotal number of Gene Names in filtered Documents:",scalar@all_gene_names,"\n" if !$nocolor;
print join(", ", sort @all_gene_names) if $print{'all_genes'} && $nocolor;
print "\n\nTotal number of Gene Names in filtered Documents:",scalar@all_gene_names,"\n" if $nocolor;

print "\n";