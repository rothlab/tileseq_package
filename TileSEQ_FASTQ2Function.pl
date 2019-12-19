#!/usr/bin/perl -w
use strict;
use Cwd qw();

# LAST edit: 2019.10.11
# LAST edited by: Roujia Li

# this script takes fastq files that are not separated into different lanes
# i.e only L001 for each fastq file

#required files (need to be in the same directory)
##file 1: geneName_seq.txt 
###file 1 description: a txt file containing template sequence for alignment and coding sequence including stop codon, the file name needs to be exactly "geneName_seq.txt" (geneName in captial letter), the file is required to have a specific format containing two lines: the first line is 'geneName_template' followed by a space or tab and then the actual template sequence including the coding sequence and upstream/downstream sequences, the second line is 'geneName_coding' followed by a space or tab and then the actual coding sequence
##file 2: tileseq_sam2mut.pl
##file 3: tileseq_mut2func.pl
##file 4: mut2func_info.csv

###three arguments required for this script###
##argument 1: geneName, e.g. CBS###
##argument 2: Qscore e.g. 20 (default value)###
##argument 3: full path for bowtie2, e.g. /home/rothlab/jweile/bin

#What does this script do###
##job1 use bowtie2-build to build a Bowtie index from the template DNA seqeucne which is provided in the geneName_seq.txt file##
##job2:merge fastq.gz files for the same indexed library, unzip all the merged fastq.gz files in the current directory into a new directory called "fastqfile" under the current directory###
##job3.1:generate .sh files for the bowtie alignment and downstream mutation call analysis into a new directory called "shfile"###
##job3.2:submit alignment jobs and mutation call jobs to SGE, mutation call jobs wait for alignment jobs to finish###
##job3.3:output sam files into a new directory called "samfile" and output all mutation call files into a new directory called "mutationcallfile"###
##job4:submit a job to call perl script tilseq_mut2func.pl to analyze the mutational call files and output results files in the directory resultfile


###preparation
###############################################################################################
###open create log files 
open(STDOUT, '>', 'running.log') or die "Can't open log";
open(STDERR, '>', 'running_error.log') or die "Can't open log";
###get current directory path
my $path = Cwd::cwd();
###creating new directories###
mkdir "bowtieIndex";
mkdir "fastqfile";
mkdir "samfile";
mkdir "shfile";
mkdir "mutationCallfile";
mkdir "resultfile";
mkdir "librarycheck";
###get geneName and user-defined Qscore value###
if(!$ARGV[0]){
  die "No argument 1";
}
my $gene = "$ARGV[0]";
if(!$ARGV[1]){
  die "No argument 2";
}
my $Qscore = "$ARGV[1]";
if($Qscore<0 || $Qscore>40){
  die "Quality score not in correct range";
}
###check if geneName agree with the geneName_seq.txt file name and create a fasta file of template sequence under the directory "bowtieIndex"
my $seqfileName = $gene.'_seq.txt';
open(my $seqfile,"<$path/$seqfileName")or die "Could not open file geneName_seq.txt input geneName doesn't match sequence fileName $!";
while(<$seqfile>){
  my @a = split;
  if(/template/){
    if(length($a[1])<30){
     die "the template sequence is too short, check the sequence file";
    }
    my $fastafileName = $gene.'.fasta';
    open(my $fasta,">$path/bowtieIndex/$fastafileName");
    print $fasta ">$gene\n$a[1]\n";
    close $fasta;
  }
}
close $seqfile;

###job1
##############################################################################################
###get the path for bowtie2###
if(!$ARGV[2]){
  die "No argument 3";
}
my $pathBowtie2 = "$ARGV[2]"; 
###get the path for bowtie index files###
my $pathIndexFiles = "$path/bowtieIndex";
###build the index files 
if(system("$pathBowtie2/bowtie2-build --quiet -f $path/bowtieIndex/$gene.fasta $path/bowtieIndex/$gene")==0){
  print "\nbowtie index files created successfully\n";
}else{
  die "bowtie indexing failed";
}


###job2
###############################################################################################
###merge the .fastq.gz files for all experiments R1/R2 reads and decompress to fastq files
my $source_dir = "$path";
my $dest_dir = "$path/fastqfile";

opendir (SOURCE, $source_dir) or die $!;
my @source_files = grep /\.gz/, readdir SOURCE;
foreach my $each_file (@source_files){
  my @splittedFileName = split(/\_/,$each_file);
  open my $file_h , "< $source_dir/$each_file";
  my @contents = <$file_h>;
  open my $dest_file, ">>$dest_dir/$splittedFileName[0]\_$splittedFileName[3]\.fastq\.gz";
  print $dest_file @contents;
  close $file_h;
  close $dest_file;
}
closedir SOURCE;

if(system("gunzip $dest_dir/*.gz")==0){
   print "\nFASTQ files merged and unzipped successfully\n";
}else{
  die "FASTQ files merging and unzipping failed";
} 
###############################################################################################


###job3
###############################################################################################
###store all experimentID in a hash###
###IMPORTANT: a new line needs to be added at the end of the command line in each .sh file for qsub to work properly!!!!###
opendir (DEST, $dest_dir) or die $!;
my @dest_files = grep /\.fastq/, readdir DEST;
my %experimentID;
foreach my $each_file (@dest_files){
  my @splittedFileName = split(/\_/,$each_file);
  $experimentID{$splittedFileName[0]} = 0;
}

my $hold_mut_qsub = '';
foreach my $key (keys %experimentID){
  my $fileNameMut = 'Mut_'.$key.'.sh';
  my $fileNameAlnR1 = 'AlnR1_'.$key.'.sh';
  my $fileNameAlnR2 = 'AlnR2_'.$key.'.sh';
  open(my $mut, ">$path/shfile/$fileNameMut");
  open(my $alnr1, ">$path/shfile/$fileNameAlnR1");
  open(my $alnr2, ">$path/shfile/$fileNameAlnR2");

  print $mut '/usr/bin/perl ', "$path/", 'tileseq_sam2mut.pl ', "$path/samfile/", $key, '_R1.sam ', "$path/samfile/", $key, '_R2.sam ', "$path/", $gene, '_seq.txt ', $Qscore, ' ', "$path/mutationCallfile", "\n"; 
  close $mut;

  print $alnr1 "$pathBowtie2/bowtie2 --quiet --no-head --no-sq --local -x $pathIndexFiles/$gene -U $path/fastqfile/$key\_R1\.fastq -S $path/samfile/$key\_R1\.sam\n";
  close $alnr1;

  print $alnr2 "$pathBowtie2/bowtie2 --quiet --no-head --no-sq --local -x $pathIndexFiles/$gene -U $path/fastqfile/$key\_R2\.fastq -S $path/samfile/$key\_R2\.sam\n";
  close $alnr2;
    
  system("qsub -N AlnR1$gene$key $path/shfile/$fileNameAlnR1\n");
  system("qsub -N AlnR2$gene$key $path/shfile/$fileNameAlnR2\n");
  system("qsub -hold_jid  AlnR1$gene$key,AlnR2$gene$key -N Mut$gene$key $path/shfile/$fileNameMut\n");
  $hold_mut_qsub .= ','."Mut$gene$key";
}
$hold_mut_qsub = substr $hold_mut_qsub,1;
###############################################################################################

###job4
###############################################################################################
###create mut2func.sh file in directory shfile
my $fileNameMut2func = 'mut2func.sh';
open(my $func,">$path/shfile/$fileNameMut2func");
print $func '/usr/bin/perl ', "$path/", 'tileseq_mut2func.pl ', "$path/", 'mut2func_info.csv ', "$path",' ',"$path/", $gene, '_seq.txt', "\n"; 
system("qsub -hold_jid $hold_mut_qsub -N Mut2Func $path/shfile/$fileNameMut2func");

