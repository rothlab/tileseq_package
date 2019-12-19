#!/usr/bin/perl -w
use strict;
use Cwd qw();

# LAST edit: 2019.10.11
# LAST edited by: Roujia Li

# min match = 130
# if the targeted region is shorter in lenth
# change this number accordingly
# on LINE 107 AND LINE 280 change
# if($totalMatch < $totalMisMatch || $totalMatch < $min_match)
# to if($totalMatch < $min_match)

###five arguments required for this script###
##argument 1: sam file for read 1
##argument 2: sam file for read 2
##argument 3: a txt file containing template sequence for alignment and coding sequence including stop codon, the file name needs to be exactly "geneName_seq.txt" (geneName in captial letter), the file is required to have a specific format containing two lines: the first line is 'geneName_template' followed by a space or tab and then the actual template sequence including the coding sequence and upstream/downstream sequences, the second line is 'geneName_coding' followed by a space or tab and then the actual coding sequence
##argument 4: phred score threshold 
##argument 5: the path for the directory where the mutation call files are saved

#######enter the template sequence (including upstream and downstream sequence),coding sequence, length of upstream and downstream sequences, desired quality score and the expected aligned fragment, and set the output directory
my $gene_template;
my $gene_coding;
my $gene;
open(IG,$ARGV[2]) or die $!;
my $seqfilefullName = "$ARGV[2]";
my @seqfileName = split(/\//,$seqfilefullName);
my @geneName = split(/\_/,$seqfileName[$#seqfileName]);
if($geneName[1] ne 'seq.txt'){
  die "sequence file name is not correct, needs to be geneName_seq.txt";
}
$gene = $geneName[0];

my $det = 0;
while(<IG>){
  my @a = split;
  if(/template/){
    $gene_template = $a[1];
    $det++;
  }
  if(/coding/){
    $gene_coding = $a[1];
    $det++;
  }
}
close IG;
if($det!=2){
  die "input sequence file format is wrong";
}
 
my $length_up = index($gene_template,$gene_coding);
my $length_down = index(reverse($gene_template),reverse($gene_coding))+3;
my $threshold_qual = "$ARGV[3]"; ##defaul value would be 20
my $min_match = 130;
my $alignedFragment = $gene; ##enter the expected aligned fragment
#####create the destinationi directory for the mutation call output files#####
if(!$ARGV[4]){
  die "No argument 4 output path";
}
my $destDirectory = "$ARGV[4]";

######################################################################################################################
open(IM,$ARGV[0]) or die $!;
my $fullName = "$ARGV[0]";
my @fileName = split(/\//,$fullName);
my @splittedFileName = split(/\_/,$fileName[$#fileName]);
my $experiment = $splittedFileName[0];
my $fileNameReport = $experiment.'report.txt';
open(my $fh, ">$destDirectory/$fileNameReport");

my $n_mapped_R1 = 0;
my $n_unmapped_R1 = 0;
my $n_mappedLow_R1 = 0;
my $n_mapped_R2 = 0;
my $n_unmapped_R2 = 0;
my $n_mappedLow_R2 = 0;

my %R1mut;
while(<IM>){
  my @a = split;
  my $length = @a;
  
  if($a[2] ne $alignedFragment){
    #print $fh 'R1',"\t",$a[0],"\t",'wrongORpoorMapping:',$a[2],"\n";
    $n_unmapped_R1++;
    next;
  }
  ####exclude those reads with more mismatches###
  my @match = split(/\d+D|\d+S|\d+I/,$a[5]);
  my $totalMatch = 0;
  foreach my $mat (@match){
    if($mat){
      if(substr($mat,-1,1) eq 'M'){
        $totalMatch += substr($mat,0,length($mat)-1);
      }
    }
  }
  my @mismatch = split(/\d+D|\d+M|\d+I/,$a[5]);
  my $totalMisMatch = 0;
  foreach my $mismat (@mismatch){
    if($mismat){
      if(substr($mismat,-1,1) eq 'S'){
        $totalMisMatch += substr($mismat,0,length($mismat)-1);
      }
    }
  }
    
  if($totalMatch < $totalMisMatch || $totalMatch < $min_match){
    $n_unmapped_R1++;
    next;
  }
  ######################################
  if($length<19){
    $n_mappedLow_R1++;
    next;
  }
  
  $n_mapped_R1++;
  ##parse cigar and quality score
  my @CIGAR = split(/\d+S/,$a[5]);
  my $clipedCIGAR = $CIGAR[@CIGAR-1];
  
  my @pos_cigar = split(/M|D|I/,$clipedCIGAR);
  my @flag_cigar = split(/\d+/,$clipedCIGAR);
  
  ##modify query sequence and quality score
  my $seq = $a[9];
  my $qual = $a[10];
  
  ###cliping $seq and $qual
  if(substr($a[5],1,1) eq 'S'){
    $seq = substr($seq,substr($a[5],0,1),length($seq)-substr($a[5],0,1));
    $qual = substr($qual,substr($a[5],0,1),length($qual)-substr($a[5],0,1));
  }
  
  if(substr($a[5],length($a[5])-1,1) eq 'S'){
    $seq = substr($seq,0,length($seq)-substr($a[5],length($a[5])-2,1));
    $qual = substr($qual,0,length($qual)-substr($a[5],length($a[5])-2,1));
  }  
  
  ###modify based on insertion and deletion
  my $position = 0;
  my $insert='';
  for(my $k=0; $k<@pos_cigar; $k++){
    if($flag_cigar[$k+1] eq 'M'){
      $position += $pos_cigar[$k];
    }
    if($flag_cigar[$k+1] eq 'D'){
      for(my $i=0; $i<$pos_cigar[$k]; $i++){
        $seq = substr($seq,0,$position).'D'.substr($seq,$position,length($seq)-$position);
        $qual = substr($qual,0,$position).'D'.substr($qual,$position,length($qual)-$position);
        $position ++;
      }  
    }
    if($flag_cigar[$k+1] eq 'I'){
      my $pos_insert = $position+$a[3]-1;
      for(my $l=0; $l<$pos_cigar[$k]; $l++){
        my $phred = ord(substr($qual,$position+$l,1))-33;
        if($phred>=$threshold_qual){
          $insert .= $pos_insert.'insert'.substr($seq,$position+$l,1);
        }
      }
      $seq = substr($seq,0,$position).substr($seq,$position+$pos_cigar[$k],length($seq)-$position-$pos_cigar[$k]);
      $qual = substr($qual,0,$position).substr($qual,$position+$pos_cigar[$k],length($qual)-$position-$pos_cigar[$k]);
    }    
  }  
  
  ##parse MD:Z:
  if($a[$length-2] =~ /MD:Z:/){
    my $seqID = $a[0];
    my @mut;
    ###print $a[0],"\t";

    my @b = split(/\:/,$a[$length-2]);
    my $MDZ = $b[2];

    my @match = split(/A|T|C|G|N|\^\D+/,$MDZ);
    my @mismatch = split(/\d+/,$MDZ);
    my $lengthMatch = @match; ###the first element is "" if the MDZ string starts with a letter
    my $lengthMismatch = @mismatch; ###the first element is "" if the MDZ string starts with a number
    
    if(substr($MDZ,0,1) =~ /\d/){
      my $pos = 0;
      for(my $i=0; $i<$lengthMismatch-1; $i++){
        if($match[$i] =~ /\d/){
          $pos += $match[$i];
        }
        if(length($mismatch[$i+1]) == 1){
          $pos++;
          my $phred = ord(substr($qual,$pos-1,1))-33;
          if($phred >= $threshold_qual){
            push(@mut,$mismatch[$i+1].($pos+$a[3]-1).substr($seq,$pos-1,1));
          }  
        }else{
          for(my $m=1; $m<length($mismatch[$i+1]); $m++){
            $pos++;
            my $phred1 = ord(substr($qual,$pos-2,1))-33;
            my $phred2 = ord(substr($qual,$pos,1))-33;
            if($phred1 >= $threshold_qual && $phred2 >= $threshold_qual){
              push(@mut,($pos+$a[3]-1).'del'.substr($mismatch[$i+1],$m,1));
            }  
          }  
        }
      }
      if($insert){
        push(@mut,$insert);
      }
    }else{
      my $pos = 0;
      for(my $i=0; $i<$lengthMatch-1; $i++){
        if(length($mismatch[$i]) == 1){
          $pos++;
          my $phred = ord(substr($qual,$pos-1,1))-33;
          if($phred >= $threshold_qual){
            push(@mut,$mismatch[$i].($pos+$a[3]-1).substr($seq,$pos-1,1));
          }  
        }else{
          for(my $m=1; $m<length($mismatch[$i]); $m++){
            $pos++;
            my $phred1 = ord(substr($qual,$pos-2,1))-33;
            my $phred2 = ord(substr($qual,$pos,1))-33;
            if($phred1 >= $threshold_qual && $phred2 >= $threshold_qual){             
              push(@mut,($pos+$a[3]-1).'del'.substr($mismatch[$i],$m,1));
            }
          }  
        }
        if($match[$i+1] =~ /\d/){
          $pos += $match[$i+1];
        }
      }
      if($lengthMismatch == $lengthMatch){
        $pos++;
        push(@mut,$mismatch[$lengthMismatch-1].($pos+$a[3]-1).substr($seq,$pos-1,1));
      }
      if($insert){
        push(@mut,$insert);
      }
    }
    
    $R1mut{$seqID} = [@mut];
    
  }else{print $fh 'R1',"\t",$a[0],"\t",'No MD:Z:',"\n";}
}
close IM;


open(IN,$ARGV[1]) or die $!;
my %R2mut;
while(<IN>){
  my @a = split;
  my $length = @a;
  
  if($a[2] ne $alignedFragment){
    $n_unmapped_R2++;
    #print $fh 'R2',"\t",$a[0],"\t",'wrongORpoorMapping:',$a[2],"\n";
    next;
  }
  ####exclude those reads with more mismatches###
  my @match = split(/\d+D|\d+S|\d+I/,$a[5]);
  my $totalMatch = 0;
  foreach my $mat (@match){
    if($mat){
      if(substr($mat,-1,1) eq 'M'){
        $totalMatch += substr($mat,0,length($mat)-1);
      }
    }
  }
  my @mismatch = split(/\d+D|\d+M|\d+I/,$a[5]);
  my $totalMisMatch = 0;
  foreach my $mismat (@mismatch){
    if($mismat){
      if(substr($mismat,-1,1) eq 'S'){
        $totalMisMatch += substr($mismat,0,length($mismat)-1);
      }
    }
  }
    
  if($totalMatch < $totalMisMatch || $totalMatch < $min_match){
    $n_unmapped_R2++;
    next;
  }
  ######################################
  if($length<19){
    $n_mappedLow_R2++;
    next;
  }
  
  $n_mapped_R2++;
  ##parse cigar and quality score
  my @CIGAR = split(/\d+S/,$a[5]);
  my $clipedCIGAR = $CIGAR[@CIGAR-1];
  
  my @pos_cigar = split(/M|D|I/,$clipedCIGAR);
  my @flag_cigar = split(/\d+/,$clipedCIGAR);
  
  ##modify query sequence and quality score
  my $seq = $a[9];
  my $qual = $a[10];
  
  ###cliping $seq and $qual
  if(substr($a[5],1,1) eq 'S'){
    $seq = substr($seq,substr($a[5],0,1),length($seq)-substr($a[5],0,1));
    $qual = substr($qual,substr($a[5],0,1),length($qual)-substr($a[5],0,1));
  }
  
  if(substr($a[5],length($a[5])-1,1) eq 'S'){
    $seq = substr($seq,0,length($seq)-substr($a[5],length($a[5])-2,1));
    $qual = substr($qual,0,length($qual)-substr($a[5],length($a[5])-2,1));
  }  
  
  ###modify based on insertion and deletion
  my $position = 0;
  my $insert='';
  for(my $k=0; $k<@pos_cigar; $k++){
    if($flag_cigar[$k+1] eq 'M'){
      $position += $pos_cigar[$k];
    }
    if($flag_cigar[$k+1] eq 'D'){
      for(my $i=0; $i<$pos_cigar[$k]; $i++){
        $seq = substr($seq,0,$position).'D'.substr($seq,$position,length($seq)-$position);
        $qual = substr($qual,0,$position).'D'.substr($qual,$position,length($qual)-$position);
        $position ++;
      }  
    }
    if($flag_cigar[$k+1] eq 'I'){
      my $pos_insert = $position+$a[3]-1;
      for(my $l=0; $l<$pos_cigar[$k]; $l++){
        my $phred = ord(substr($qual,$position+$l,1))-33;
        if($phred>=$threshold_qual){
          $insert .= $pos_insert.'insert'.substr($seq,$position+$l,1);
        }
      }
      $seq = substr($seq,0,$position).substr($seq,$position+$pos_cigar[$k],length($seq)-$position-$pos_cigar[$k]);
      $qual = substr($qual,0,$position).substr($qual,$position+$pos_cigar[$k],length($qual)-$position-$pos_cigar[$k]);
    }    
  }  
  
  ##parse MD:Z:
  if($a[$length-2] =~ /MD:Z:/){
    my $seqID = $a[0];
    my @mut;
    ###print $a[0],"\t";

    my @b = split(/\:/,$a[$length-2]);
    my $MDZ = $b[2];

    my @match = split(/A|T|C|G|N|\^\D+/,$MDZ);
    my @mismatch = split(/\d+/,$MDZ);
    my $lengthMatch = @match; ###the first element is "" if the MDZ string starts with a letter
    my $lengthMismatch = @mismatch; ###the first element is "" if the MDZ string starts with a number
    
    if(substr($MDZ,0,1) =~ /\d/){
      my $pos = 0;
      for(my $i=0; $i<$lengthMismatch-1; $i++){
        if($match[$i] =~ /\d/){
          $pos += $match[$i];
        }
        if(length($mismatch[$i+1]) == 1){
          $pos++;
          my $phred = ord(substr($qual,$pos-1,1))-33;
          if($phred >= $threshold_qual){
            push(@mut,$mismatch[$i+1].($pos+$a[3]-1).substr($seq,$pos-1,1));
          }  
        }else{
          for(my $m=1; $m<length($mismatch[$i+1]); $m++){
            $pos++;
            my $phred1 = ord(substr($qual,$pos-2,1))-33;
            my $phred2 = ord(substr($qual,$pos,1))-33;
            if($phred1 >= $threshold_qual && $phred2 >= $threshold_qual){
              push(@mut,($pos+$a[3]-1).'del'.substr($mismatch[$i+1],$m,1));
            }  
          }  
        }
      }
      if($insert){
        push(@mut,$insert);
      }
    }else{
      my $pos = 0;
      for(my $i=0; $i<$lengthMatch-1; $i++){
        if(length($mismatch[$i]) == 1){
          $pos++;
          my $phred = ord(substr($qual,$pos-1,1))-33;
          if($phred >= $threshold_qual){
            push(@mut,$mismatch[$i].($pos+$a[3]-1).substr($seq,$pos-1,1));
          }  
        }else{
          for(my $m=1; $m<length($mismatch[$i]); $m++){
            $pos++;
            my $phred1 = ord(substr($qual,$pos-2,1))-33;
            my $phred2 = ord(substr($qual,$pos,1))-33;
            if($phred1 >= $threshold_qual && $phred2 >= $threshold_qual){             
              push(@mut,($pos+$a[3]-1).'del'.substr($mismatch[$i],$m,1));
            }
          }  
        }
        if($match[$i+1] =~ /\d/){
          $pos += $match[$i+1];
        }
      }
      if($lengthMismatch == $lengthMatch){
        $pos++;
        push(@mut,$mismatch[$lengthMismatch-1].($pos+$a[3]-1).substr($seq,$pos-1,1));
      }
      if($insert){
        push(@mut,$insert);
      }
    }
    
    $R2mut{$seqID} = [@mut];
    
  }else{print $fh 'R2',"\t",$a[0],"\t",'No MD:Z:',"\n";}
}
close IN;



###get the mutations shared by the R1 and R2 reads and interpret the mutation results

my %hash_AAC=();
my %hash_AAC_multiple=();
my %hash_DEL=();
my %hash_INS=();
#my %hash_NON=();

foreach my $key (keys %R1mut){
#######call SNPs
#######################################################################################################################  
  if($R2mut{$key}){
   my @MUT_R1 = @{$R1mut{$key}};
   my @MUT_R2 = @{$R2mut{$key}};
   if(@MUT_R1 && @MUT_R2){
    ##modify the $seq_ref
    my $seq_ref = $gene_coding;
    my $seq_mut1 = $seq_ref;
    my $seq_mut2 = $seq_ref;
    ####modify $seq_ref1
    foreach (@MUT_R1){
      my $det=0;
      if($_ =~ /del/ || $_ =~/insert/){
        $det++;
      }else{
        my $pos = substr($_,1,length($_)-2);
        if($pos<=$length_up || $pos>length($gene_template)-$length_down){
          $det++;
        }
      }
      if($det==0){
        my $pos = substr($_,1,length($_)-2)-$length_up;
        my $nt_ref = substr($_,0,1);
        my $nt_mut = substr($_,length($_)-1,1);
        if($nt_ref eq substr($seq_ref,$pos-1,1)){
          if($pos <= length($seq_ref)){
              $seq_mut1 = substr($seq_mut1,0,$pos-1).$nt_mut.substr($seq_mut1,$pos,length($seq_mut1)-$pos);
          }else{
            print $fh $key,"\t","SNP position is out range","\n";
          }
        }else{
          print $fh $key,"\t","ref nt is wrong","\n";
        }
      }  
    }
    ####modify $seq_ref2
    foreach (@MUT_R2){
      my $det=0;
      if($_ =~ /del/ || $_ =~/insert/){
        $det++;
      }else{
        my $pos = substr($_,1,length($_)-2);
        if($pos<=$length_up || $pos>length($gene_template)-$length_down){
          $det++;
        }
      }
      if($det==0){
        my $pos = substr($_,1,length($_)-2)-$length_up;
        my $nt_ref = substr($_,0,1);
        my $nt_mut = substr($_,length($_)-1,1);
        if($nt_ref eq substr($seq_ref,$pos-1,1)){
          if($pos <= length($seq_ref)){
              $seq_mut2 = substr($seq_mut2,0,$pos-1).$nt_mut.substr($seq_mut2,$pos,length($seq_mut2)-$pos);
          }else{
            print $fh $key,"\t","SNP position is out range","\n";
          }
        }else{
          print $fh $key,"\t","ref nt is wrong","\n";
        }
      }  
    }
        
    ##translate the $seq_ref and $seq_mut1/$seq_mut2, and compare
    if(length($seq_ref) == length($seq_mut1) && length($seq_ref) == length($seq_mut2)){
      my $aarefs;
      my $positions;
      my $aamuts;
      my $codonrefs;
      my $codonmuts;
      my $count=0;
      for(my $i=0; $i<length($seq_ref)/3; $i++){
          my $codon_ref = substr($seq_ref,$i*3,3);
          my $codon_mut1 = substr($seq_mut1,$i*3,3);
          my $codon_mut2 = substr($seq_mut2,$i*3,3);
          my $aa_ref = &codon2aa($codon_ref);
          my $aa_mut1 = &codon2aa($codon_mut1);
          my $aa_mut2 = &codon2aa($codon_mut2);
          if($codon_ref ne $codon_mut1 && $codon_ref ne $codon_mut2 && $codon_mut1 eq $codon_mut2){
              my $aachange = $aa_ref."\t".($i+1)."\t".$aa_mut1."\t".$codon_ref."\t".$codon_mut1;
              $hash_AAC{$aachange} += 1;
              $count++;
              $aarefs .= ($aa_ref.'|');
              $positions .= (($i+1).'|');
              $aamuts .= ($aa_mut1.'|');
              $codonrefs .= ($codon_ref.'|');
              $codonmuts .= ($codon_mut1.'|');  
          }
      }
      if($count>1){
        my $aachangeMultiple = $aarefs."\t".$positions."\t".$aamuts."\t".$codonrefs."\t".$codonmuts;
        $hash_AAC_multiple{$aachangeMultiple} += 1;
      } 
    }else{
      print $fh $key,"\t","ref modification error\n";
    }
   } 
  }
#######################################################################################################################  
 
#######################################################################################################################  
######call deletion and insertion#####
  if($R2mut{$key}){
   my @MUT_R1 = @{$R1mut{$key}};
   my @MUT_R2 = @{$R2mut{$key}};
   if(@MUT_R1 && @MUT_R2){
    my $deletions='';
    foreach my $mut1 (@MUT_R1){
      foreach my $mut2 (@MUT_R2){
        if($mut1 eq $mut2){
          if($mut1 =~ /del/){
            my $pos = substr($mut1,0,length($mut1)-4);
            my $ntDEL = substr($mut1,length($mut1)-1,1);
            if($ntDEL eq substr($gene_template,$pos-1,1)){
              $deletions .= substr($mut1,0,length($mut1)-4).'del'.substr($mut1,length($mut1)-1,1);
              ##$hash_DEL{substr($mut1,0,length($mut1)-4)."\t".'del'."\t".substr($mut1,length($mut1)-1,1)} += 1;
            }
          }elsif($mut1 =~ /insert/){
            $hash_INS{substr($mut1,0,length($mut1)-7)."\t".'ins'."\t".substr($mut1,length($mut1)-1,1)} += 1;
          }
        }
      }
    }
    if($deletions){
      $hash_DEL{$deletions} += 1;
    }
   } 
  }
####################################################################################################################################### 
}

print $fh "The number of mapped R1 reads is: $n_mapped_R1\n";
print $fh "The number of unmapped R1 reads is: $n_unmapped_R1\n";
print $fh "The number of mapped R1 reads with low quality is: $n_mappedLow_R1\n";
print $fh "The number of mapped R2 reads is: $n_mapped_R2\n";
print $fh "The number of unmapped R2 reads is: $n_unmapped_R2\n";
print $fh "The number of mapped R2 reads with low quality is: $n_mappedLow_R2\n";

###calucualte the sequencing depth (million reads) from the number of mapped R1 and R2 reads, which will be used to normalize the mutation counts###
my $depth = 0;
if($n_mapped_R1<$n_mapped_R2){
  $depth = $n_mapped_R1/1000000;
}else{
  $depth = $n_mapped_R2/1000000;
 }
if($depth==0){
  die("sequencing depth is ZERO!!!");
}
print $fh "The sequencing depth is: $depth million reads\n";
close $fh;

###print out the mutations and normalized mutation counts###
my $fileNameAAchange = $experiment.'AAchange.txt';
my $fileNameAAchangeMultiple = $experiment.'MultipleMut.txt';
my $fileNameDel = $experiment.'deletion.txt';
my $fileNameIns = $experiment.'insertion.txt';
#my $fileNameNon = $experiment.'noncoding.txt';

open(my $AAchange, ">$destDirectory/$fileNameAAchange");
open(my $AAchangeMultiple, ">$destDirectory/$fileNameAAchangeMultiple");
open(my $deletion, ">$destDirectory/$fileNameDel");
open(my $insertion, ">$destDirectory/$fileNameIns");
#open(my $noncoding, ">$destDirectory/$fileNameNon");


foreach my $mut (keys %hash_AAC){
  print $AAchange "$mut\t",$hash_AAC{$mut}/$depth,"\n";
}
foreach my $mut (keys %hash_AAC_multiple){
  print $AAchangeMultiple "$mut\t",$hash_AAC_multiple{$mut}/$depth,"\n";
}
foreach my $mut (keys %hash_DEL){
  print $deletion "$mut\t",$hash_DEL{$mut}/$depth,"\n";
}
foreach my $mut (keys %hash_INS){
  print $insertion "$mut\t",$hash_INS{$mut}/$depth,"\n";
}
#foreach my $mut (keys %hash_NON){
  #print $noncoding "$mut\t",$hash_NON{$mut}/$depth,"\n";
#}

close $AAchange;
close $AAchangeMultiple;
close $deletion;
close $insertion;
#close $noncoding;
######################################################################################################################

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
