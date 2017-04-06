#!/usr/bin/perl -w
use strict;
use POSIX;

###five arguments required for this script###
##argument 1: mut2func_info.csv file. This file contains seven lines. The first line contains the corresponding amino acid positions for each tile. The rest six lines contains the sequencing sample name for each multiplexed librarys sorted into each of the six experiments: nonselect1, nonselect2, select1, select2, wt1, wt2
##argument 2: the path for the working directory 
##argument 3: geneName_seq.txt file

######get the directory path where the AAchange.txt files are stored
my $dir = "$ARGV[1]";
###enter the minimum requirement for normalized mutant abundance in non-selective condition after substraction
my $n_sub = 0;
###enter the sequencing depth used for normalization
my $depth_norm = 1000000;

###enter the experiment ID for each experiment from the first region to the last region
my %regions_nonpriming;
my %exp;
open(IF,$ARGV[0]) or die "can't open mut2func_info.csv file $!";
my $exp_det = 0;
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
      push(@{$exp{1}},$exp); 
    }
  }
  if($first eq 'nonselect2'){
    $exp_det ++;
    foreach my $exp (@a) {
      push(@{$exp{2}},$exp); 
    }
  }
  if($first eq 'select1'){
    $exp_det ++;
    foreach my $exp (@a) {
      push(@{$exp{3}},$exp); 
    }
  }
  if($first eq 'select2'){
    $exp_det ++;
    foreach my $exp (@a) {
      push(@{$exp{4}},$exp); 
    }
  }
  if($first eq 'wt1'){
    $exp_det ++;
    foreach my $exp (@a) {
      push(@{$exp{5}},$exp); 
    }
  }
  if($first eq 'wt2'){
    $exp_det ++;
    foreach my $exp (@a) {
      push(@{$exp{6}},$exp); 
    }
  }
}
close IF;

if($exp_det == 6){
  ###assign experiment ID arrays to a hash
  my @exp1 = @{$exp{1}}; ##for nonselect1
  my @exp2 = @{$exp{2}}; ##for nonselect2
  my @exp3 = @{$exp{3}}; ##for select1
  my @exp4 = @{$exp{4}}; ##for select2
  my @exp5 = @{$exp{5}}; ##for wt1
  my @exp6 = @{$exp{6}}; ##for wt2
  my $n_regions = @exp1;
  ###generate nonpriming region for each experiment and store in a hash
  my %nonPrimingRegion;
  for(my $i=1; $i<=6; $i++){
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

my %control1raw;
my %control2raw;
my %control1;
my %control2;
my %average_control;
my %sd_control;

#####################################################
###open directory and read all file names in an array, then read each file and assign six hashes for nonselective, selective and control experiments
opendir(my $dh, "$dir/mutationCallfile") or die $!;
my @filenames = readdir $dh;
closedir $dh;

foreach my $filename (@filenames){
  if($filename =~ /\S+AAchange.txt/){
   my @expID = split(/AAchange/,$filename);
   if($nonPrimingRegion{$expID[0]}){
    open(IM,"$dir/mutationCallfile/$filename") or die "can't open file $filename $!";
    while(<IM>){
      my @a = split;
      my $count = $a[5];
      my $countraw = $a[5];
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
      
      ###for wt-control replicate 1
      if($expID[0] ~~ @exp5){
        if($pos ~~ [$region[0]..$region[1]]){
          $control1{$mut} += $count;
          $control1raw{$mut_codon} += $countraw;
        }
      }
      
      ###for wt-control replicate 2
      if($expID[0] ~~ @exp6){
        if($pos ~~ [$region[0]..$region[1]]){
          $control2{$mut} += $count;
          $control2raw{$mut_codon} += $countraw;
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

open(my $fhctr,">$dir/resultfile/correlation_control.txt") or die $!;
print $fhctr "control1\tcontrol2\n";
foreach my $mut (keys %control1){
  if($control2{$mut}){
    print $fhctr "$control1{$mut}\t$control2{$mut}\n";
  }
}
close $fhctr;
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
foreach my $mut_codon (keys %control1raw){
  $hash_all{$mut_codon} = 0;
}
foreach my $mut_codon (keys %control2raw){
  $hash_all{$mut_codon} = 0;
}


print $fhraw "wt_aa\tpos\tmut_aa\twt_codon\tmut_codon\tannotation\tnonselect1\tnonselect2\tselect1\tselect2\tcontrol1\tcontrol2\n";
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
  if($control1raw{$mut_codon}){
    print $fhraw $control1raw{$mut_codon},"\t";
  }else{
    print $fhraw "0\t";
  }
  if($control2raw{$mut_codon}){
    print $fhraw $control2raw{$mut_codon},"\n";
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
    }
}

foreach my $mut (keys %select1){
    if($select2{$mut}){
      $average_select{$mut} = ($select1{$mut}+$select2{$mut})/2;
      $sd_select{$mut} = stdev(($select1{$mut},$select2{$mut}));
    }
}

foreach my $mut (keys %control1){
    if($control2{$mut}){
      $average_control{$mut} = ($control1{$mut}+$control2{$mut})/2;
      $sd_control{$mut} = stdev(($control1{$mut},$control2{$mut}));
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
  if($average_select{$mut}){
    if($average_control{$mut}){
      if((3*$sd_control{$mut}+$average_control{$mut}) < $average_nonselect{$mut} && (3*$sd_control{$mut}+$average_control{$mut}) < $average_select{$mut}){
        $nonselect1_sub{$mut} = $nonselect1{$mut} - $average_control{$mut};
        $nonselect2_sub{$mut} = $nonselect2{$mut} - $average_control{$mut};
        $nonselect_sub{$mut} = $average_nonselect{$mut} - $average_control{$mut};
        $select1_sub{$mut} = $select1{$mut} - $average_control{$mut};
        $select2_sub{$mut} = $select2{$mut} - $average_control{$mut};
        $select_sub{$mut} = $average_select{$mut} - $average_control{$mut};
        my @nonselect = ($nonselect1{$mut},$nonselect2{$mut});
        my @select = ($select1{$mut},$select2{$mut});
        my @control = ($control1{$mut},$control2{$mut});
      
        $sd_nonselect_sub{$mut} = ($sd_nonselect{$mut}**2 + $sd_control{$mut}**2 )**0.5;
        $sd_select_sub{$mut} = ($sd_select{$mut}**2 + $sd_control{$mut}**2 )**0.5;
      }
    }else{
      $nonselect_sub{$mut} = $average_nonselect{$mut};
      $nonselect1_sub{$mut} = $nonselect1{$mut};
      $nonselect2_sub{$mut} = $nonselect2{$mut};
      $select_sub{$mut} = $average_select{$mut};
      $select1_sub{$mut} = $select1{$mut};
      $select2_sub{$mut} = $select2{$mut};
      $sd_nonselect_sub{$mut} = $sd_nonselect{$mut};
      $sd_select_sub{$mut} = $sd_select{$mut};
    }
  }
}
#####################################################

#####################################################
###calculate fold change
my %foldchange;
my %sd_foldchange;
foreach my $mut (keys %nonselect_sub){
  $foldchange{$mut} = $select_sub{$mut}/$nonselect_sub{$mut};
  my @nonselect_sub = ($nonselect1_sub{$mut},$nonselect2_sub{$mut});
  my @select_sub = ($select1_sub{$mut},$select2_sub{$mut});
  $sd_foldchange{$mut} = $foldchange{$mut}*((($sd_select_sub{$mut}/$select_sub{$mut})**2 + ($sd_nonselect_sub{$mut}/$nonselect_sub{$mut})**2)**0.5);
}

open(my $fh,">$dir/resultfile/foldchange.txt") or die $!;
print $fh "wt\tpos\tmut\tannotation\taverageNonselect\tsdNonselect\tcvNonselect\tsdSelect\tcvSelect\tsdFoldchange\tcvFoldchange\tfoldchange\n";
foreach my $mut (keys %nonselect_sub){
  my $CVnonselect = $sd_nonselect_sub{$mut}/$nonselect_sub{$mut};
  my $CVselect = $sd_select_sub{$mut}/$select_sub{$mut};
  my $CVfoldchange = $sd_foldchange{$mut}/$foldchange{$mut};
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
for(my $i=1; $i<=$region_last[1]; $i++){
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

###define empty hashes for POPs with >2 nt change at nt level
my %replicate1POP;
my %replicate2POP;
my %average_replicatePOP;
my %sd_replicatePOP;

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
      my $mut_codon = $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]."\t".$a[4]."\t".$mutationType;
      
      my $mut_pop = '';
      my $det_pop = 0;
      for(my $i=0; $i<3; $i++){
        if(substr($a[3],$i,1) eq substr($a[4],$i,1)){
          $det_pop++;
        }
      }
      if($det_pop <2){
        $mut_pop = $mut_codon;
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
        }
      }
      
      ###for pre-select replicate 2
      if($expID[0] ~~ @exp2){
        if($pos ~~ [$region[0]..$region[1]]){
          $replicate2{$mut} += $count;
          $replicate2nt{$mut_codon} += $count;
          if($mut_pop){
            $replicate2POP{$mut_pop} += $count;
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
my $prob_noDel = 1;
my $mutation_total=0;
my $POP_total=0;

print $fhfreq "region\tmut_freq\tdel_freq\tseq_depth (million reads)\n";
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
  my $seqdep_nonprim = sprintf("%.2g",$seqdep_nonpriming{$i});
  print $fhfreq $region[0],'~',$region[1],"\t$mut_counts\t$del_counts\t$seqdep_nonprim\n";
  $prob_noDel *= (1-$del_counts);
  $length_aa += ($region[1]-$region[0]+1);
  
}
my $num_nonsyn_f100 = keys(%nonsyn_f100);
my $num_nonsyn_f50 = keys(%nonsyn_f50);
my $num_nonsyn_f = keys(%nonsyn_f);

my $n_nonsyn_f_accessible = 0;
my $n_nonsyn_f_inaccessible = 0;

foreach my $mut (keys %nonsyn_f){
  my @a = split(/\t/,$mut);
  my $aa_change = $a[0].$a[1].$a[2];
  if(AAaccessible_sub1nt($seq_nt,$aa_change) eq 'accessible'){
    $n_nonsyn_f_accessible++;
  }else{
    $n_nonsyn_f_inaccessible++;
  }
}
print $fhfreq "\n";
print $fhfreq "Average AA changes per clone\t",sprintf("%.2g",$f_nonsyn),"\n";
print $fhfreq "Fraction of measured AA changes\t",sprintf("%.2g",$num_nonsyn_f/($length_aa*19)),"\n";
print $fhfreq "Fraction of AA changes with >1 nt change\t",sprintf("%.2g",$POP_total/$mutation_total),"\n";
print $fhfreq "Probability of a clone carring a frameshift mutation\t",sprintf("%.2g",(1-$prob_noDel)),"\n";

print $fhfreq "Total AA chagnes\t", $n_allAAchange,"\n";
print $fhfreq "1-nt accessible AA chagnes\t", $n_accessible,"\n";
print $fhfreq "1-nt inaccessible AA chagnes\t", $n_inaccessible,"\n";
print $fhfreq "Fraction of measured 1-nt accessible AA changes\t",sprintf("%.2g",$n_nonsyn_f_accessible/$n_accessible),"\n";
print $fhfreq "Fraction of measured 1-nt inaccessible AA changes\t",sprintf("%.2g",$n_nonsyn_f_inaccessible/$n_inaccessible),"\n";

print $fhfreq "Average AA changes per clone with f'>100\t",sprintf("%.2g",$fprime100_nonsyn),"\n";
print $fhfreq "Fraction of measured AA changes with f'>100\t",sprintf("%.2g",$num_nonsyn_f100/($length_aa*19)),"\n";
print $fhfreq "Average AA changes per clone with f'>50\t",sprintf("%.2g",$fprime50_nonsyn),"\n";
print $fhfreq "Fraction of measured AA changes with f'>50 is:\t",sprintf("%.2g",$num_nonsyn_f50/($length_aa*19)),"\n";


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

sub average{
        my @data = @_;
        if (not @data) {
               die("Empty arrayn");
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

