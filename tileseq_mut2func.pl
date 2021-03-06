#!/usr/bin/perl -w
use strict;
use POSIX;

# LAST edit: 2019.12.10
# LAST edited by: Roujia Li

###five arguments required for this script###
##argument 1: mut2func_info.csv file. This file contains seven lines. The first line contains the corresponding amino acid positions for each tile. The rest eight lines contains the sequencing sample name for each multiplexed librarys sorted into each of the six experiments: nonselect1, nonselect2, select1, select2, wtNonselect1, wtNonselect2, wtSelect1, wtSelect2
##argument 2: the path for the working directory
##argument 3: geneName_seq.txt file

######get the directory path where the AAchange.txt files are stored
my $dir = "$ARGV[1]";
###enter the minimum requirement for normalized mutant abundance in non-selective condition after substraction
my $n_sub = 0;
###enter the sequencing depth used for normalization
my $depth_norm = 1000000;

###enter the experiment ID for each experiment from the first region to the last region and generate a output file for sequencing depth
my %regions_nonpriming;
my %exp;
open(IF,$ARGV[0]) or die "can't open mut2func_info.csv file $!";
my $exp_det = 0;
open(my $fhdep,">$dir/resultfile/sequencingDepth.csv") or die $!;
print $fhdep 'SampleID,SequencingDepth (million reads)',"\n";
while(<IF>){
  chomp;
  my @a = split(/\,/,);
  my $first = shift(@a);
###enter the nonpriming regions
  if($first eq 'tiles'){
    my $n_tile = 0;
    foreach my $tile (@a){
      $n_tile++;
      my @pos = split(/\-/,$tile);
      $regions_nonpriming{$n_tile} = [($pos[0],$pos[1])];
    }
  }
###assign the experiment ID
  if($first eq 'nonselect1'){
    $exp_det ++;
    foreach my $exp (@a) {
      my $filename = $exp.'report.txt';
      open(ID,"<$dir/mutationCallfile/$filename") or die "can't open $filename $!";
      print $filename."\n";
      while(<ID>){
        if(/depth/){
          my @a = split;
          print $fhdep "SampleID$exp\_NS1\,$a[4]\n";
        }
      }
      close ID;
      push(@{$exp{1}},$exp);
    }
  }
  if($first eq 'nonselect2'){
    $exp_det ++;
    foreach my $exp (@a) {
      my $filename = $exp.'report.txt';
      open(ID,"<$dir/mutationCallfile/$filename") or die "can't open $filename $!";
      while(<ID>){
        if(/depth/){
          my @a = split;
          print $fhdep "SampleID$exp\_NS2\,$a[4]\n";
        }
      }
      close ID;
      push(@{$exp{2}},$exp);
    }
  }
  if($first eq 'select1'){
    $exp_det ++;
    foreach my $exp (@a) {
      my $filename = $exp.'report.txt';
      open(ID,"<$dir/mutationCallfile/$filename") or die "can't open $filename $!";
      while(<ID>){
        if(/depth/){
          my @a = split;
          print $fhdep "SampleID$exp\_S1\,$a[4]\n";
        }
      }
      close ID;
      push(@{$exp{3}},$exp);
    }
  }
  if($first eq 'select2'){
    $exp_det ++;
    foreach my $exp (@a) {
      my $filename = $exp.'report.txt';
      open(ID,"<$dir/mutationCallfile/$filename") or die "can't open $filename $!";
      while(<ID>){
        if(/depth/){
          my @a = split;
          print $fhdep "SampleID$exp\_S2\,$a[4]\n";
        }
      }
      close ID;
      push(@{$exp{4}},$exp);
    }
  }
  if($first eq 'wtNonselect1'){
    $exp_det ++;
    foreach my $exp (@a) {
      my $filename = $exp.'report.txt';
      open(ID,"<$dir/mutationCallfile/$filename") or die "can't open $filename $!";
      while(<ID>){
        if(/depth/){
          my @a = split;
          print $fhdep "SampleID$exp\_NONSELECTWT1\,$a[4]\n";
        }
      }
      close ID;
      push(@{$exp{5}},$exp);
    }
  }
  if($first eq 'wtNonselect2'){
    $exp_det ++;
    foreach my $exp (@a) {
      my $filename = $exp.'report.txt';
      open(ID,"<$dir/mutationCallfile/$filename") or die "can't open $filename $!";
      while(<ID>){
        if(/depth/){
          my @a = split;
          print $fhdep "SampleID$exp\_NONSELECTWT2\,$a[4]\n";
        }
      }
      close ID;
      push(@{$exp{6}},$exp);
    }
  }
  if($first eq 'wtSelect1'){
    $exp_det ++;
    foreach my $exp (@a) {
      my $filename = $exp.'report.txt';
      open(ID,"<$dir/mutationCallfile/$filename") or die "can't open $filename $!";
      while(<ID>){
        if(/depth/){
          my @a = split;
          print $fhdep "SampleID$exp\_SELECTWT1\,$a[4]\n";
        }
      }
      close ID;
      push(@{$exp{7}},$exp);
    }
  }
  if($first eq 'wtSelect2'){
    $exp_det ++;
    foreach my $exp (@a) {
      my $filename = $exp.'report.txt';
      open(ID,"<$dir/mutationCallfile/$filename") or die "can't open $filename $!";
      while(<ID>){
        if(/depth/){
          my @a = split;
          print $fhdep "SampleID$exp\_SELECTWT2\,$a[4]\n";
        }
      }
      close ID;
      push(@{$exp{8}},$exp);
    }
  }
}
close IF;
close $fhdep;

if($exp_det == 8){
  ###assign experiment ID arrays to a hash
  my @exp1 = @{$exp{1}}; ##for nonselect1
  my @exp2 = @{$exp{2}}; ##for nonselect2
  my @exp3 = @{$exp{3}}; ##for select1
  my @exp4 = @{$exp{4}}; ##for select2
  my @exp5 = @{$exp{5}}; ##for nonselectwt1
  my @exp6 = @{$exp{6}}; ##for nonselectwt2
  my @exp7 = @{$exp{7}}; ##for selectwt1
  my @exp8 = @{$exp{8}}; ##for selectwt2
  my $n_regions = @exp1;
  ###generate nonpriming region for each experiment and store in a hash
  my %nonPrimingRegion;
  for(my $i=1; $i<=8; $i++){
    for(my $j=1; $j<=$n_regions; $j++){
      my @expID = @{$exp{$i}};
      my @region = @{$regions_nonpriming{$j}};
      $nonPrimingRegion{$expID[$j-1]} = [@region];
    }
  }

###define empty hashes
my %nonselect1raw;
my %nonselect2raw;
my %nonselect1;
my %nonselect2;
my %average_nonselect;
my %sd_nonselect;

my %select1raw;
my %select2raw;
my %select1;
my %select2;
my %average_select;
my %sd_select;

my %controlNS1raw;
my %controlNS2raw;
my %controlNS1;
my %controlNS2;
my %average_controlNS;
my %sd_controlNS;

my %controlS1raw;
my %controlS2raw;
my %controlS1;
my %controlS2;
my %average_controlS;
my %sd_controlS;

#####################################################
###open directory and read all file names in an array, then read each file and assign six hashes for nonselective, selective and control experiments
opendir(my $dh, "$dir/mutationCallfile") or die $!;
my @filenames = readdir $dh;
closedir $dh;

foreach my $filename (@filenames){
  if($filename =~ /\S+AAchange.txt/){ ## parse AA change file, it has all the codon changes from SNP
   my @expID = split(/AAchange/,$filename);
   if($nonPrimingRegion{$expID[0]}){
    open(IM,"$dir/mutationCallfile/$filename") or die "can't open file $filename $!";
    while(<IM>){
      my @a = split;
      my $count = $a[5];
      my $countraw = $a[5];
      my $pos = $a[1];
      my $mutationType;
      ## assign mutation types to different mutations
      if($a[0] eq $a[2]){
        $mutationType = 'SYN';
      }elsif($a[2] eq '_'){
        $mutationType = 'STOP';
      }else{
        $mutationType = 'NONSYN';
      }

      my $mut = $a[0]."\t".$a[1]."\t".$a[2]."\t".$mutationType;
      my $mut_codon = $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]."\t".$a[4]."\t".$mutationType;
      my @region = @{$nonPrimingRegion{$expID[0]}};

      ###for pre-select replicate 1
      if($expID[0] ~~ @exp1){
        if($pos ~~ [$region[0]..$region[1]]){
          $nonselect1{$mut} += $count;
          $nonselect1raw{$mut_codon} += $countraw;
        }
      }

      ###for pre-select replicate 2
      if($expID[0] ~~ @exp2){
        if($pos ~~ [$region[0]..$region[1]]){
          $nonselect2{$mut} += $count;
          $nonselect2raw{$mut_codon} += $countraw;
        }
      }

      ###for post-select replicate 1
      if($expID[0] ~~ @exp3){
        if($pos ~~ [$region[0]..$region[1]]){
          $select1{$mut} += $count;
          $select1raw{$mut_codon} += $countraw;
        }
      }

      ###for post-select replicate 2
      if($expID[0] ~~ @exp4){
        if($pos ~~ [$region[0]..$region[1]]){
          $select2{$mut} += $count;
          $select2raw{$mut_codon} += $countraw;
        }
      }

      ###for wt-control Nonselect replicate 1
      if($expID[0] ~~ @exp5){
        if($pos ~~ [$region[0]..$region[1]]){
          $controlNS1{$mut} += $count;
          $controlNS1raw{$mut_codon} += $countraw;
        }
      }

      ###for wt-control Nonselect replicate 2
      if($expID[0] ~~ @exp6){
        if($pos ~~ [$region[0]..$region[1]]){
          $controlNS2{$mut} += $count;
          $controlNS2raw{$mut_codon} += $countraw;
        }
      }

      ###for wt-control select replicate 1
      if($expID[0] ~~ @exp7){
        if($pos ~~ [$region[0]..$region[1]]){
          $controlS1{$mut} += $count;
          $controlS1raw{$mut_codon} += $countraw;
        }
      }

      ###for wt-control select replicate 2
      if($expID[0] ~~ @exp8){
        if($pos ~~ [$region[0]..$region[1]]){
          $controlS2{$mut} += $count;
          $controlS2raw{$mut_codon} += $countraw;
        }
      }
    }
    close IM;
   }
  }

}
#####################################################

####################################################
###output the file for correlation analysis between replicates
open(my $fhpre,">$dir/resultfile/correlation_nonselect.txt") or die $!;
print $fhpre "nonselect1\tnonselect2\n";
foreach my $mut (keys %nonselect1){
  if($nonselect2{$mut}){
    print $fhpre "$nonselect1{$mut}\t$nonselect2{$mut}\n";
  }
}
close $fhpre;

open(my $fhpost,">$dir/resultfile/correlation_select.txt") or die $!;
print $fhpost "select1\tselect2\n";
foreach my $mut (keys %select1){
  if($select2{$mut}){
    print $fhpost "$select1{$mut}\t$select2{$mut}\n";
  }
}
close $fhpost;

open(my $fhctrNS,">$dir/resultfile/correlation_controlNS.txt") or die $!;
print $fhctrNS "controlNS1\tcontrolNS2\n";
foreach my $mut (keys %controlNS1){
  if($controlNS2{$mut}){
    print $fhctrNS "$controlNS1{$mut}\t$controlNS2{$mut}\n";
  }
}
close $fhctrNS;

open(my $fhctrS,">$dir/resultfile/correlation_controlS.txt") or die $!;
print $fhctrS "controlS1\tcontrolS2\n";
foreach my $mut (keys %controlS1){
  if($controlS2{$mut}){
    print $fhctrS "$controlS1{$mut}\t$controlS2{$mut}\n";
  }
}
close $fhctrS;
####################################################

####################################################
###output the file for raw data
open(my $fhraw,">$dir/resultfile/rawData.txt") or die $!;
my %hash_all;

foreach my $mut_codon (keys %nonselect1raw){
  $hash_all{$mut_codon} = 0;
}
foreach my $mut_codon (keys %nonselect2raw){
  $hash_all{$mut_codon} = 0;
}
foreach my $mut_codon (keys %select1raw){
  $hash_all{$mut_codon} = 0;
}
foreach my $mut_codon (keys %select2raw){
  $hash_all{$mut_codon} = 0;
}
foreach my $mut_codon (keys %controlNS1raw){
  $hash_all{$mut_codon} = 0;
}
foreach my $mut_codon (keys %controlNS2raw){
  $hash_all{$mut_codon} = 0;
}
foreach my $mut_codon (keys %controlS1raw){
  $hash_all{$mut_codon} = 0;
}
foreach my $mut_codon (keys %controlS2raw){
  $hash_all{$mut_codon} = 0;
}

print $fhraw "wt_aa\tpos\tmut_aa\twt_codon\tmut_codon\tannotation\tnonselect1\tnonselect2\tselect1\tselect2\tcontrolNS1\tcontrolNS2\tcontrolS1\tcontrolS2\n";
foreach my $mut_codon (keys %hash_all){
  print $fhraw "$mut_codon\t";
  if($nonselect1raw{$mut_codon}){
    print $fhraw $nonselect1raw{$mut_codon},"\t";
  }else{
    print $fhraw "0\t";
  }
  if($nonselect2raw{$mut_codon}){
    print $fhraw $nonselect2raw{$mut_codon},"\t";
  }else{
    print $fhraw "0\t";
  }
  if($select1raw{$mut_codon}){
    print $fhraw $select1raw{$mut_codon},"\t";
  }else{
    print $fhraw "0\t";
  }
  if($select2raw{$mut_codon}){
    print $fhraw $select2raw{$mut_codon},"\t";
  }else{
    print $fhraw "0\t";
  }
  if($controlNS1raw{$mut_codon}){
    print $fhraw $controlNS1raw{$mut_codon},"\t";
  }else{
    print $fhraw "0\t";
  }
  if($controlNS2raw{$mut_codon}){
    print $fhraw $controlNS2raw{$mut_codon},"\t";
  }else{
    print $fhraw "0\t";
  }
  if($controlS1raw{$mut_codon}){
    print $fhraw $controlS1raw{$mut_codon},"\t";
  }else{
    print $fhraw "0\t";
  }
  if($controlS2raw{$mut_codon}){
    print $fhraw $controlS2raw{$mut_codon},"\n";
  }else{
    print $fhraw "0\n";
  }
}
close $fhraw;
####################################################

####################################################
###calcualte the average normalized mutational counts
foreach my $mut (keys %nonselect1){
    if($nonselect2{$mut}){
      $average_nonselect{$mut} = ($nonselect1{$mut}+$nonselect2{$mut})/2;
      $sd_nonselect{$mut} = stdev(($nonselect1{$mut},$nonselect2{$mut}));
      if(!$select1{$mut}){$select1{$mut}=0;}
      if(!$select2{$mut}){$select2{$mut}=0;}
      if(!$controlNS1{$mut}){$controlNS1{$mut}=0;}
      if(!$controlNS2{$mut}){$controlNS2{$mut}=0;}
      if(!$controlS1{$mut}){$controlS1{$mut}=0;}
      if(!$controlS2{$mut}){$controlS2{$mut}=0;}
      $average_select{$mut} = ($select1{$mut}+$select2{$mut})/2;
      $sd_select{$mut} = stdev(($select1{$mut},$select2{$mut}));
      $average_controlNS{$mut} = ($controlNS1{$mut}+$controlNS2{$mut})/2;
      $sd_controlNS{$mut} = stdev(($controlNS1{$mut},$controlNS2{$mut}));
      $average_controlS{$mut} = ($controlS1{$mut}+$controlS2{$mut})/2;
      $sd_controlS{$mut} = stdev(($controlS1{$mut},$controlS2{$mut}));
    }
}

#####################################################

#####################################################
###substract the normalized mutant abundance in control
my %nonselect1_sub;
my %nonselect2_sub;
my %nonselect_sub;
my %select1_sub;
my %select2_sub;
my %select_sub;
my %sd_nonselect_sub;
my %sd_select_sub;

foreach my $mut (keys %average_nonselect){

    ###calculate nonselect_sub
    if($average_nonselect{$mut}>$average_controlNS{$mut}+3*$sd_controlNS{$mut}){
      $nonselect_sub{$mut} = $average_nonselect{$mut} - $average_controlNS{$mut};
      $sd_nonselect_sub{$mut} = ($sd_nonselect{$mut}**2 + $sd_controlNS{$mut}**2 )**0.5;
    }

    ####calculate select_sub
    if($average_select{$mut}+3*$sd_select{$mut}>=$average_controlNS{$mut}-3*$sd_controlNS{$mut}){
      $select_sub{$mut} = abs($average_select{$mut} - $average_controlS{$mut});
      $sd_select_sub{$mut} = ($sd_select{$mut}**2 + $sd_controlS{$mut}**2 )**0.5;
    }

}
#####################################################

#####################################################
###calculate fold change
my %foldchange;
my %sd_foldchange;
foreach my $mut (keys %nonselect_sub){
  if($select_sub{$mut}){
    $foldchange{$mut} = $select_sub{$mut}/$nonselect_sub{$mut};
    if($select_sub{$mut}==0){
      $sd_foldchange{$mut} = $sd_nonselect_sub{$mut};
    }else{
      $sd_foldchange{$mut} = $foldchange{$mut}*((($sd_select_sub{$mut}/$select_sub{$mut})**2 + ($sd_nonselect_sub{$mut}/$nonselect_sub{$mut})**2)**0.5);
    }
  }
}

open(my $fh,">$dir/resultfile/foldchange.txt") or die $!;
print $fh "wt\tpos\tmut\tannotation\taverageNonselect\tsdNonselect\tcvNonselect\tsdSelect\tcvSelect\tsdFoldchange\tcvFoldchange\tfoldchange\n";
foreach my $mut (keys %foldchange){
  my $CVnonselect = $sd_nonselect_sub{$mut}/$nonselect_sub{$mut};
  my $CVselect;
  my $CVfoldchange;
  if($select_sub{$mut} != 0){
    $CVselect = $sd_select_sub{$mut}/$select_sub{$mut};
    $CVfoldchange = $sd_foldchange{$mut}/$foldchange{$mut};
  }else{
    $CVselect = 0;
    $CVfoldchange = 0;
  }
  if($nonselect_sub{$mut}>$n_sub){
    print $fh "$mut\t$nonselect_sub{$mut}\t$sd_nonselect_sub{$mut}\t$CVnonselect\t$sd_select_sub{$mut}\t$CVselect\t$sd_foldchange{$mut}\t$CVfoldchange\t$foldchange{$mut}\n";
  }
}
close $fh;
#####################################################

#####################################################
###calculate nonsynonymous frequency
my $totalCounts;

foreach my $mut (keys %nonselect_sub){
  if($average_nonselect{$mut}>$n_sub){
    my @a = split(/\t/,$mut);
    if($a[3] eq 'NONSYN'){
      $totalCounts += $nonselect_sub{$mut};
    }
  }
}

print 'The nonsynonymous frequency is: ',$totalCounts/$depth_norm,"\n";
#####################################################
}elsif($exp_det == 2){
  ###assign experiment ID arrays to a hash
  my @exp1 = @{$exp{1}}; ##for nonselect1
  my @exp2 = @{$exp{2}}; ##for nonselect2
  my $gene_template;
  my $gene_coding;
  my $geneName;
  my $length_up;
  my $length_down;
  ### assign geneName
  my $seqfilefullName = "$ARGV[2]";
  my @seqfileName = split(/\//,$seqfilefullName);
  my @Name = split(/\_/,$seqfileName[$#seqfileName]);
  $geneName = $Name[0];
  ###open the geneName_seq.txt file and assign sequence and upstream length
  open(IN,$ARGV[2]) or die "Could not open file geneName_seq.txt $!";
  while(<IN>){
    my @a = split;
    if(/template/){
      $gene_template = $a[1];
    }
    if(/coding/){
      $gene_coding = $a[1];
    }
  }
  close IN;
  $length_up = index($gene_template,$gene_coding);
  $length_down = index(reverse($gene_template),reverse($gene_coding))+3;
  my $seq_nt = $gene_coding;
###enter the experimental condition and the mutation frequency threshold###
my $condition = 'nonselect';
my $f_threshold = 4; ####the threshold is calculated as (1/sequencingdepth)*$f_threshold####
my $f_POP_threshold = 1; ####the threshold is calculated as (1/sequencingdepth)*$f_POP_threshold####

###create a hash for sequencing depth and enter the sequence depth for each tile
my $n_regions = @exp1;
my %seqdep_nonpriming;
for(my $i=0; $i<$n_regions; $i++){
  my $filename = $exp1[$i].'report.txt';
  open(IR,"<$dir/mutationCallfile/$filename");
  while(<IR>){
    if(/depth/){
      my @a = split;
      $seqdep_nonpriming{$i+1} = $a[4];
    }
  }
  close IR;
}
for(my $j=1; $j<=$n_regions; $j++){
  print 'The sequencing depth of region ', $j, ' is: ', $seqdep_nonpriming{$j},"\n";
}

##############################################################
###generate amino acid sequences for the give gene
my $seq_AA;
for(my $i=0; $i<length($seq_nt)/3; $i++){
  if(codon2aa(substr($seq_nt,3*$i,3)) ne '_'){
    $seq_AA .= codon2aa(substr($seq_nt,3*$i,3));
  }
}
###for the examined region, generate a hash with all the possible AA changes in the format as three-letter string as (aa_wt,position,aa_mut), for example, D5W
my %allAAchange;
my @allAA = ('A','V','L','I','M','F','Y','W','R','H','K','D','E','S','T','N','Q','G','C','P');
my @region_last = @{$regions_nonpriming{$n_regions}};
my @region_first = @{$regions_nonpriming{1}};
for(my $i=$region_first[0]; $i<=$region_last[1]; $i++){
    my $aa_wt = substr($seq_AA,$i-1,1);
    foreach my $aa_sub (@allAA){
        if($aa_wt ne $aa_sub){
            my $aa_change = $aa_wt.$i.$aa_sub;
            $allAAchange{$aa_change} = AAaccessible_sub1nt($seq_nt,$aa_change);
        }
    }
}
my $n_allAAchange = keys(%allAAchange);
my $n_accessible = 0;
my $n_inaccessible = 0;
foreach my $mut (keys %allAAchange){
  if($allAAchange{$mut} eq 'accessible'){
    $n_accessible++;
  }else{
    $n_inaccessible++;
  }
}
##############################################################

###generate nonpriming region for each experiment and store in a hash
my %nonPrimingRegion;
my $n_exp = $exp_det;
for(my $i=1; $i<=$n_exp; $i++){
  for(my $j=1; $j<=$n_regions; $j++){
    my @expID = @{$exp{$i}};
    my @region = @{$regions_nonpriming{$j}};
    $nonPrimingRegion{$expID[$j-1]} = [@region];
  }
}

###define empty hashes for deletions
my %replicate1Del;
my %replicate2Del;
my %average_replicateDel;
my %sd_replicateDel;

###define empty hashes for mutations at AA level
my %replicate1;
my %replicate2;
my %average_replicate;
my %sd_replicate;

###define empty hashes for mutations at nt level
my %replicate1nt;
my %replicate2nt;
my %average_replicatent;
my %sd_replicatent;

###define empty hashes for POPs with >1 nt change at nt level
my %replicate1POP;
my %replicate2POP;
my %average_replicatePOP;
my %sd_replicatePOP;

###define empty hashes for POPs with >1 nt change at nt level
my %replicate1nt_acc; # accessible by 1nt
my %replicate2nt_acc;
my %average_replicatent_acc;
my %sd_replicatent_acc;

###define empty hashes for positional mutation counts for all 1-nt, 2-nt and 3-nt mutations, and deletions###
my %posMut1;
my %posMut2;
my %posMut3;
my %posDel;
my %posNonsyn;
for(my $j=1; $j<=$n_regions; $j++){
  my @region = @{$regions_nonpriming{$j}};
  for(my $i=$region[0]; $i<=$region[1]; $i++){
    $posMut1{$i} = 0;
    $posMut2{$i} = 0;
    $posMut3{$i} = 0;
    $posDel{$i} = 0;
    $posNonsyn{$i} = 0;
  }
}

#####################################################
###open directory and read all file names in an array, then read each file and assign hashes
opendir(my $dh, "$dir/mutationCallfile");
my @filenames = readdir $dh;
closedir $dh;

foreach my $filename (@filenames){
  #####for mutations#####
  if($filename =~ /\S+AAchange.txt/){
   my @expID = split(/AAchange/,$filename);
   if($nonPrimingRegion{$expID[0]}){
    open(IM,"$dir/mutationCallfile/$filename");
    while(<IM>){
      my @a = split;
      my $count = $a[5];
      my $pos = $a[1];
      my $mutationType;
      if($a[0] eq $a[2]){
        $mutationType = 'SYN';
      }elsif($a[2] eq '_'){
        $mutationType = 'STOP';
      }else{
        $mutationType = 'NONSYN';
      }

      my $mut = $a[0]."\t".$a[1]."\t".$a[2]."\t".$mutationType;
      #
      my $mut_codon = $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]."\t".$a[4]."\t".$mutationType;

      my $mut_pop = '';
      my $mut_1nt = '';
      my $det_pop = 0;
      for(my $i=0; $i<3; $i++){
        if(substr($a[3],$i,1) eq substr($a[4],$i,1)){
          $det_pop++;
        }
      }
      if($det_pop <2){ # if less than two nt mapped (means there is >1nt change)
        $mut_pop = $mut_codon;
      } elsif ($det_pop = 2){ # if two nt mapped (means there is 1nt change)
        $mut_1nt = $mut_codon;
      }

      my @region = @{$nonPrimingRegion{$expID[0]}};

      ###for pre-select replicate 1
      if($expID[0] ~~ @exp1){
        if($pos ~~ [$region[0]..$region[1]]){
          $replicate1{$mut} += $count;
          $replicate1nt{$mut_codon} += $count;
          if($mut_pop){
            $replicate1POP{$mut_pop} += $count;
          }
          if($mut_1nt){
            $replicate1nt_acc{$mut_1nt} += $count;
          }
        }
      }

      ###for pre-select reicate 2
      if($expID[0] ~~ @exp2){
        if($pos ~~ [$region[0]..$region[1]]){
          $replicate2{$mut} += $count;
          $replicate2nt{$mut_codon} += $count;
          if($mut_pop){
            $replicate2POP{$mut_pop} += $count;
          }
          if($mut_1nt){
            $replicate2nt_acc{$mut_1nt} += $count;
          }

        }
      }

    }
    close IM;
   }
  }

  #####for deletions#####
  if($filename =~ /\S+deletion.txt/){
   my @expID = split(/deletion/,$filename);
   if($nonPrimingRegion{$expID[0]}){
    open(IM,"$dir/mutationCallfile/$filename");
    while(<IM>){
      my @a = split;
      my $count = $a[1];
      my @pos = split(/del\D/,$a[0]);
      my $pos_max = max(@pos)-$length_up;
      my $pos_min = min(@pos)-$length_up;
      my $mut = $a[0];

      my @region = @{$nonPrimingRegion{$expID[0]}};

      ###for pre-select replicate 1
      if($expID[0] ~~ @exp1){
        if($pos_min ~~ [(3*$region[0]-2)..(3*$region[1])] && $pos_max ~~ [(3*$region[0]-2)..(3*$region[1])]){
          $replicate1Del{$mut} += $count;
        }
      }

      ###for pre-select replicate 2
      if($expID[0] ~~ @exp2){
        if($pos_min ~~ [(3*$region[0]-2)..(3*$region[1])] && $pos_max ~~ [(3*$region[0]-2)..(3*$region[1])]){
          $replicate2Del{$mut} += $count;
        }
      }

    }
    close IM;
   }
  }
}
#####################################################

####################################################
###calcualte the average normalized mutational and deletion counts
foreach my $mut (keys %replicate1){
    if($replicate2{$mut}){
      $average_replicate{$mut} = ($replicate1{$mut}+$replicate2{$mut})/2;
      $sd_replicate{$mut} = stdev(($replicate1{$mut},$replicate2{$mut}));
    }
}

foreach my $mut (keys %replicate1nt){
    if($replicate2nt{$mut}){
      $average_replicatent{$mut} = ($replicate1nt{$mut}+$replicate2nt{$mut})/2;
      $sd_replicatent{$mut} = stdev(($replicate1nt{$mut},$replicate2nt{$mut}));
    }
}

foreach my $mut (keys %replicate1Del){
    if($replicate2Del{$mut}){
      $average_replicateDel{$mut} = ($replicate1Del{$mut}+$replicate2Del{$mut})/2;
      $sd_replicateDel{$mut} = stdev(($replicate1Del{$mut},$replicate2Del{$mut}));
    }
}

foreach my $mut (keys %replicate1POP){
    if($replicate2POP{$mut}){
      $average_replicatePOP{$mut} = ($replicate1POP{$mut}+$replicate2POP{$mut})/2;
      $sd_replicatePOP{$mut} = stdev(($replicate1POP{$mut},$replicate2POP{$mut}));
    }
}

# 1-nt accessible change
foreach my $mut (keys %replicate1nt_acc){
    if($replicate2nt_acc{$mut}){
      $average_replicatent_acc{$mut} = ($replicate1nt_acc{$mut}+$replicate2nt_acc{$mut})/2;
      $sd_replicatent_acc{$mut} = stdev(($replicate1nt_acc{$mut},$replicate2nt_acc{$mut}));
    }
}
#####################################################

#####################################################
###output four files for mutation_aa, mutation_nt, mutation_pop and deletion
open(my $fhAA,">$dir/librarycheck/mutAA_$geneName\_$condition\.txt");
print $fhAA "wt\tpos\tmut\tannotation\taverage\tsd\n";
foreach my $mut (keys %average_replicate){
  print $fhAA "$mut\t$average_replicate{$mut}\t$sd_replicate{$mut}\n";
}
close $fhAA;

open(my $fhNT,">$dir/librarycheck/mutCodon_$geneName\_$condition\.txt");
print $fhNT "wt\tpos\tmut\twtCodon\tmutCodon\tannotation\taverage\tsd\n";
foreach my $mut (keys %average_replicatent){
  print $fhNT "$mut\t$average_replicatent{$mut}\t$sd_replicatent{$mut}\n";
}
close $fhNT;

open(my $fhPOP,">$dir/librarycheck/mutPOP_$geneName\_$condition\.txt");
print $fhPOP "wt\tpos\tmut\twtCodon\tmutCodon\tannotation\taverage\tsd\n";
foreach my $mut (keys %average_replicatePOP){
  print $fhPOP "$mut\t$average_replicatePOP{$mut}\t$sd_replicatePOP{$mut}\n";
}
close $fhPOP;

open(my $fh1nt,">$dir/librarycheck/mut1Nt_$geneName\_$condition\.txt");
print $fh1nt "wt\tpos\tmut\twtCodon\tmutCodon\tannotation\taverage\tsd\n";
foreach my $mut (keys %average_replicatent_acc){
  print $fh1nt "$mut\t$average_replicatent_acc{$mut}\t$sd_replicatent_acc{$mut}\n";
}
close $fh1nt;

open(my $fhDEL,">$dir/librarycheck/newdel_$geneName\_$condition\.txt");
print $fhDEL "deletion\taverage\tsd\n";
foreach my $mut (keys %average_replicateDel){
  print $fhDEL "$mut\t$average_replicateDel{$mut}\t$sd_replicateDel{$mut}\n";
}
close $fhDEL;
#####################################################

#####################################################
###calculate mutation and deletion frequency
open(my $fhfreq,">$dir/librarycheck/reportCard_$geneName\_$condition\.txt");
my $fprime100_nonsyn = 0;
my $fprime50_nonsyn = 0;
my $f_nonsyn = 0;
my $length_aa = 0;
my %nonsyn_f100;
my %nonsyn_f50;
my %nonsyn_f;
my %nonsyn_all;
my $prob_noDel = 1;
my $mutation_total=0;
my $POP_total=0;

my $onent_total=0;
my $f_nonsen=1;
print $fhfreq "Tiles\tMut_freq\tDel_freq\tSequencing_depth (thousands reads)\n";
my @start_aa;
my @end_aa;
for(my $i=1; $i<=$n_regions; $i++){
  my @region = @{$regions_nonpriming{$i}};
  my $seq_dep = $seqdep_nonpriming{$i};
  my $count_threshold = 1/$seq_dep;
  my $mut_counts = 0;
  my $mut_counts_f100 = 0;
  my $del_counts = 0;
  foreach my $mut (keys %average_replicate){
    my @a = split(/\t/,$mut);
    if($a[1] ~~ [$region[0]..$region[1]]){
      $mut_counts += $average_replicate{$mut}/$depth_norm;
      if($average_replicate{$mut}>=100 && $a[3] eq 'NONSYN'){
        $fprime100_nonsyn += $average_replicate{$mut}/$depth_norm;
        $nonsyn_f100{$mut} = $average_replicate{$mut};
      }
      if($average_replicate{$mut}>=50 && $a[3] eq 'NONSYN'){
        $fprime50_nonsyn += $average_replicate{$mut}/$depth_norm;
        $nonsyn_f50{$mut} = $average_replicate{$mut};
      }
    }
  }

  # added for 1nt change (RL)
  foreach my $mut (keys %average_replicatent_acc){
    my @a = split(/\t/,$mut);
    if($a[1] ~~ [$region[0]..$region[1]]){
      # print($a[5],"\n");
      if($average_replicatent_acc{$mut}>=50 && $a[5] eq 'NONSYN'){
        $onent_total += $average_replicatent_acc{$mut};
      }
    }
  }

  foreach my $mut (keys %average_replicatent){
    my @a = split(/\t/,$mut);
    my $mutAA = $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[5];
    my $det_pop = 0;
    for(my $i=0; $i<3; $i++){
      if(substr($a[3],$i,1) eq substr($a[4],$i,1)){
        $det_pop++;
      }
    }
    if($a[1] ~~ [$region[0]..$region[1]]){
      if($a[5] eq 'NONSYN'){
        if($average_replicatent{$mut} > ceil($f_threshold*$count_threshold) || ($det_pop<2 && $average_replicatent{$mut}>ceil($f_POP_threshold*$count_threshold))){
          $f_nonsyn += $average_replicatent{$mut}/$depth_norm;
          $nonsyn_f{$mutAA} += $average_replicatent{$mut};
        }
        $nonsyn_all{$mutAA} += $average_replicatent{$mut};
      }

      # added for nonsense mutations (RL)
      if($a[5] eq 'STOP'){
        if($average_replicatent{$mut}>=50){
          $f_nonsen += $average_replicatent{$mut}/$depth_norm;
        }
      }
    }
  }

  foreach my $mut (keys %average_replicatePOP){
    my @a = split(/\t/,$mut);
    if($a[1] ~~ [$region[0]..$region[1]]){
      $POP_total += $average_replicatePOP{$mut};
    }
  }
  foreach my $mut (keys %average_replicate){
    my @a = split(/\t/,$mut);
    if($a[1] ~~ [$region[0]..$region[1]]){
      $mutation_total += $average_replicate{$mut};
    }
  }

  foreach my $mut (keys %average_replicateDel){
    my @pos = split(/del\D/,$mut);
    my $pos_max = max(@pos)-$length_up;
    my $pos_min = min(@pos)-$length_up;
    if($pos_min ~~ [(3*$region[0]-2)..(3*$region[1])] && $pos_max ~~ [(3*$region[0]-2)..(3*$region[1])]){
      $del_counts += $replicate2Del{$mut}/$depth_norm;
    }
  }
  $mut_counts = sprintf("%.2g",$mut_counts);
  $del_counts = sprintf("%.2g",$del_counts);
  my $seqdep_nonprim = sprintf("%.2g",$seqdep_nonpriming{$i} *1000);
  print $fhfreq $region[0],'~',$region[1],"\t$mut_counts\t$del_counts\t$seqdep_nonprim\n";
  $prob_noDel *= (1-$del_counts);
  push(@start_aa,$region[0]);
  push(@end_aa,$region[1]);

}
for(my $i=0; $i<$n_regions; $i++){
  $length_aa += ($end_aa[$i]-$start_aa[$i]+1);
}

my $num_nonsyn_f100 = keys(%nonsyn_f100);
my $num_nonsyn_f50 = keys(%nonsyn_f50);
my $num_nonsyn_f = keys(%nonsyn_f);
my $num_nonsyn_all = keys(%nonsyn_all);

my $n_nonsyn_f_accessible = 0;
my $n_nonsyn_f_inaccessible = 0;

my $n_nonsyn_all_accessible = 0;
my $n_nonsyn_all_inaccessible = 0;

foreach my $mut (keys %nonsyn_f){
  my @a = split(/\t/,$mut);
  my $aa_change = $a[0].$a[1].$a[2];
  if(AAaccessible_sub1nt($seq_nt,$aa_change) eq 'accessible'){
    $n_nonsyn_f_accessible++;
  }else{
    $n_nonsyn_f_inaccessible++;
  }
}

foreach my $mut (keys %nonsyn_all){
  my @a = split(/\t/,$mut);
  my $aa_change = $a[0].$a[1].$a[2];
  if(AAaccessible_sub1nt($seq_nt,$aa_change) eq 'accessible'){
    $n_nonsyn_all_accessible++;
  }else{
    $n_nonsyn_all_inaccessible++;
  }
}

print $fhfreq "\n";

print $fhfreq "frequency_(mut)  = (Total number of reads with mutation) / (Total number of mapped reads in sample)\n";
print $fhfreq "frequency per million = frequency * 10^6\n\n";
# for(my $i=0; $i<$n_regions; $i++){
#   print $fhfreq $start_aa[$i],'-',$end_aa[$i],',';
# }
# explain the columns
print $fhfreq "mut_freq = Sum of frequencies of all mutations in the given tile\n";
print $fhfreq "del_freq = Sum of frequencies of all deletions in the given tile\n";
print $fhfreq "seq_depth = Total number of reads mapped / 1000\n";

print $fhfreq "\n";
print $fhfreq '1. The total number of possible AA changes is ', $length_aa*19,"\n";
print $fhfreq "   = Total number of AA in the input tiles * 19 \n";
print $fhfreq "     1.a 1-nt accessible AA changes\t", $n_accessible,"\n";
print $fhfreq "         = Total number of possible AA changes that can be made by changing one nt\n";
print $fhfreq "     1.b 1-nt inaccessible AA changes\t", $n_inaccessible,"\n";
print $fhfreq "         = Total number of possible AA changes that can NOT be made by changing one nt\n\n";

print $fhfreq "2. Average AA changes per clone\t",sprintf("%.2g",$f_nonsyn),"\n";
print $fhfreq "   = Sum of frequency of all the well-measured non-syn mutations (nt level) in the input tiles\n\n";

print $fhfreq "3. Probability of a clone carrying deletion\t",sprintf("%.2g",(1-$prob_noDel)),"\n";
print $fhfreq "   = sum of the frequency of all the deletions\n\n";

print $fhfreq "4. Probability that clone carrying a nonsense mutation\t",sprintf("%.2g",$f_nonsen),"\n";
print $fhfreq "   = sum of the frequency of all the well-measured nonsense mutations\n\n";

print $fhfreq "5. Fraction of all possible AA changes that are well measured\t",sprintf("%.2g",$num_nonsyn_f/($length_aa*19)),"\n";
print $fhfreq "   = (Total number of well-measured non-syn changes found) / (Total possible AA changes)\n\n";

print $fhfreq "6. Fraction of all possible 1 nt-accessible AA changes that are well measured\t",sprintf("%.2g",$n_nonsyn_f_accessible/$n_accessible),"\n";
print $fhfreq "   = (Number of well-measured 1-nt accessible AA changes found) / (Total number of possible 1nt-accessible AA changes)\n\n";

print $fhfreq "7. Fraction of all possible 1 nt-inaccessible AA changes that are well measured\t",sprintf("%.2g",$n_nonsyn_f_inaccessible/$n_inaccessible),"\n";
print $fhfreq "   = (Number of well-measured 1-nt inaccessible AA changes found (nt-level)) / (Total number of possible 1nt-inaccessible AA changes)\n\n";

print $fhfreq "8. Fraction of all possible AA change seen at least once.\t",sprintf("%.2g",$num_nonsyn_all/($length_aa*19)),"\n";
print $fhfreq "   = (Total number of Non-syn changes found) / (Total possible AA changes)\n\n";

print $fhfreq "9. Fraction of all possible 1 nt-accessible AA change seen at least once.\t",sprintf("%.2g",$n_nonsyn_all_accessible/$n_accessible),"\n";
print $fhfreq "   = (Total number of 1 nt-accessible AA changes) / (Total number of possible 1nt-accessible AA changes)\n\n";

print $fhfreq "10. Fraction of AA changes (well measured) with 1 nt change \t",sprintf("%.2g",$onent_total/$mutation_total),"\n";
print $fhfreq "   = (sum of frequency per million of all the well-measured aa change with 1nt change) / (sum of frequency per million of all the mutations we found)\n\n";


#print $fhfreq "Average AA changes per clone with f'>100\t",sprintf("%.2g",$fprime100_nonsyn),"\n";
#print $fhfreq "Fraction of measured AA changes with f'>100\t",sprintf("%.2g",$num_nonsyn_f100/($length_aa*19)),"\n";
#print $fhfreq "Average AA changes per clone with f'>50\t",sprintf("%.2g",$fprime50_nonsyn),"\n";
#print $fhfreq "Fraction of measured AA changes with f'>50 is:\t",sprintf("%.2g",$num_nonsyn_f50/($length_aa*19)),"\n";


close $fhfreq;
#####################################################

#####################################################
###fill up the four positonal mutation/deletion counts hashes###
foreach my $mut (keys %average_replicatent){
  my @a = split(/\t/,$mut);

  my $seq_dep;
  for(my $i=1; $i<=$n_regions; $i++){
    my @region = @{$regions_nonpriming{$i}};
    if($a[1] ~~ [$region[0]..$region[1]]){
      $seq_dep = $seqdep_nonpriming{$i};
    }
  }
  my $count_threshold = 1/$seq_dep;

  my $det_pop = 0;
  for(my $i=0; $i<3; $i++){
    if(substr($a[3],$i,1) eq substr($a[4],$i,1)){
      $det_pop++;
    }
  }
  if($det_pop==0 && $average_replicatent{$mut}>ceil($f_POP_threshold*$count_threshold)){
    $posMut3{$a[1]} += $average_replicatent{$mut};
  }
  if($det_pop==1 && $average_replicatent{$mut}>ceil($f_POP_threshold*$count_threshold)){
    $posMut2{$a[1]} += $average_replicatent{$mut};
  }
  if($det_pop==2 && $average_replicatent{$mut}>ceil($f_threshold*$count_threshold)){
    $posMut1{$a[1]} += $average_replicatent{$mut};
  }
}

foreach my $mut (keys %average_replicateDel){
  my @pos = split(/del\D/,$mut);
  foreach my $pos (@pos){
    my $posAA = ceil(($pos-$length_up)/3);
    $posDel{$posAA} += $average_replicateDel{$mut};
  }
}

foreach my $mut (keys %nonsyn_f){
  my @a = split(/\t/,$mut);
  $posNonsyn{$a[1]} += 1;

}
#####################################################

#####################################################
###output file for positional mutation/deletion distribution###
open(my $fhland,">$dir/librarycheck/landscape_$geneName\_$condition\.txt");
print $fhland "pos\tnt1\tnt>1\tnt2\tnt3\tdel\tnumberNONSYN\n";
foreach my $pos (keys %posDel){
  print $fhland $pos,"\t";
  print $fhland $posMut1{$pos},"\t";
  print $fhland $posMut3{$pos}+$posMut2{$pos},"\t";
  print $fhland $posMut2{$pos},"\t";
  print $fhland $posMut3{$pos},"\t";
  print $fhland $posDel{$pos},"\t";
  print $fhland $posNonsyn{$pos},"\n";
}
close $fhland;
#####################################################
}else{
  die "the number of experiments is not correct";
}


##############################################################
###define sub for median, takes in an array
sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
##############################################################
sub average{
        my @data = @_;
        if (not @data) {
               die("Empty array");
        }
        my $total = 0;
        foreach (@data) {
                $total += $_;
        }
        my $average = $total / @data;
        return $average;
}
sub stdev{
        my @data = @_;
        if(@data == 1){
                return 0;
        }
        my $average = average(@data);
        my $sqtotal = 0;
        foreach(@data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / @data) ** 0.5;
        return $std;
}
sub covariance {
    my ($array1ref, $array2ref) = @_;
    if (not @$array1ref) {
        die("Empty array1");
    }
    if (not @$array2ref) {
        die("Empty array2");
    }
    if (@$array2ref != @$array1ref) {
        die("unequal length of arrays");
    }

    my $result;

    for(my $i=0; $i<@$array1ref; $i++){
      $result += ($array1ref->[$i]-average(@$array1ref))*($array2ref->[$i]-average(@$array2ref));
    }
    $result /= @$array1ref;
    return $result;

}
sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
}

##############################################################
###define sub for min function, takes in a list
sub min {
    my ($min, @vars) = @_;
    for (@vars) {
        $min = $_ if $_ < $min;
    }
    return $min;
}
##############################################################

##############################################################
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
exit;
}
}
##############################################################

##############################################################
###define sub for accessible amino acid changes by changing one nucleotide positions
###This subroutine needs to call another subroutine codon2aa
###AAaccessible_sub1nt(string1,string2), string1 is the coding sequence of the gene, string2 is the three-letter aa change (aa_wt,position,aa_mut), for example, D5W
sub AAaccessible_sub1nt{
  ####read in coding sequence and the amino acid change
  my ($seq_coding,$mut) = @_;
  ####get the aa pos
  my $pos = substr($mut,1,length($mut)-2);
  ####get the wt_aa, mut_aa and wt_codon
  my $aa_mut = substr($mut,-1,1);
  my $aa_wt = substr($mut,0,1);
  my $codon_wt = substr($seq_coding,3*($pos-1),3);
  ####mutate each nt the other three nts and check if the input mutation is accessible####
  my $det = 0;
  for(my $i=0; $i<3; $i++){
    my $sub;
    ##replace with A and check
    if(substr($codon_wt,$i,1) ne 'A'){
      my $codon_mutToA = $codon_wt;
      $sub = substr($codon_mutToA,$i,1,'A');
      my $aa_mutToA = &codon2aa($codon_mutToA);
      if($aa_mutToA eq $aa_mut){
        $det++;
      }
    }
    ##replace with T and check
    if(substr($codon_wt,$i,1) ne 'T'){
      my $codon_mutToT = $codon_wt;
      $sub = substr($codon_mutToT,$i,1,'T');
      my $aa_mutToT = &codon2aa($codon_mutToT);
      if($aa_mutToT eq $aa_mut){
        $det++;
      }
    }
    ##replace with C and check
    if(substr($codon_wt,$i,1) ne 'C'){
      my $codon_mutToC = $codon_wt;
      $sub = substr($codon_mutToC,$i,1,'C');
      my $aa_mutToC = &codon2aa($codon_mutToC);
      if($aa_mutToC eq $aa_mut){
        $det++;
      }
    }
    ##replace with G and check
    if(substr($codon_wt,$i,1) ne 'G'){
      my $codon_mutToG = $codon_wt;
      $sub = substr($codon_mutToG,$i,1,'G');
      my $aa_mutToG = &codon2aa($codon_mutToG);
      if($aa_mutToG eq $aa_mut){
        $det++;
      }
    }
  }

  if($det==0){
    return 'inaccessible';
  }else{
    return 'accessible';
  }

}
##############################################################
