#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;
use File::Path;
use Cwd 'abs_path';
use FindBin;
use lib $FindBin::Bin;

#========Variables defination===============#
my $bin_dir=dirname($0);
$bin_dir = abs_path($bin_dir) or die $!;
my $dir=substr($bin_dir,0,rindex($bin_dir,"\/"));
my $in_csv=undef;
my $genome="hg19";
my $fasta=$dir."\/fasta\/hg19.fa";
my $out_dir="./"; #default: current_dir
my $out_name="manifest"; # default: new_manifest.csv & outdated_manifest.csv
my $help=0;


#========Get options================#
unless (
    GetOptions(
        "in=s"   => \$in_csv,
        "g=s"    => \$genome,
        "o=s"    => \$out_dir,
        "n=s"    => \$out_name,
        "h|help" =>\$help
    )
  )
{
    print $!, "\n";
    _usage(1);
}

#==========Run Program==============================#
if ($help){
    _usage();
    exit(0);
}

_check();
_main();
END {
    if($?){
        print "Program exist with Error \n";
    }
}

#============Checking=====================#
sub _check{
    #1) in file checking  
    if ( !defined($in_csv) ) {
    	_error("input manifest file (-in) has not been defined yet\n");
    	_usage(1);
    }else{
    	if ( !-e $in_csv ) {
            _error("input manifest file (-in) '$in_csv' does not exists\n");
            _usage(1);
      }
      $in_csv = abs_path($in_csv);
    }

    #2) genome checking
    if(($genome ne "hg19") and ($genome ne "mm10") and ($genome ne "GRCh38") and ($genome ne "mm9")){
  		_error("Please set genome (-g) as GRCh38, hg19, mm9, or mm10!\n");
  	  _usage(1);
  	}
  	if($genome eq "GRCh38"){
		  $fasta=$dir."\/fasta\/GRCh38.fa";
  		$fasta = abs_path($fasta) or die $!;
  	}
  	if($genome eq "hg19"){
		  $fasta=$dir."\/fasta\/hg19.fa";
  		$fasta = abs_path($fasta) or die $!;
  	}
  	if($genome eq "mm9"){
		  $fasta=$dir."\/fasta\/mm9.fa";
  		$fasta = abs_path($fasta) or die $!;
  	}
  	if($genome eq "mm10"){
		  $fasta=$dir."\/fasta\/mm10.fa";
  		$fasta = abs_path($fasta) or die $!;
  	}
    
    #3) out directory checking
    if ( !defined($out_dir) ) {
    	_error("work directory (-o) has not been defined yet\n");
    	_usage(1);
    }else{
    	if ( !-e $out_dir ) {
            _error("work directory (-o) '$out_dir' does not exists\n");
            _usage(1);
      }
      $out_dir = abs_path($out_dir);
      if(!($out_dir=~/\/$/)){$out_dir.="\/";}
    }
}
    

#============Usage=====================#
sub _usage{
    my ($flag) = @_;
    my $usage=<<USAGE;
Usage: perl code_step1_sw.pl [options] -in manifest.csv
e.g: perl code_step1_sw.pl -o test/ -n OncoChip -g hg19 -in OncoChip.csv 
options:
-g [string]             define the genome: GRCh38, hg19, mm9, or mm10 (default: hg19)
-o [string]             output directory (default: current directory)
-in [string]            required,  input manifest file in csv format
-n [string]             name for output
-h                      this help message

USAGE

    print $usage;
    exit(1) if ($flag);
}


sub _error{
    my ( $s, $flag ) = @_;
    print STDERR"[", scalar(localtime), "] [ERROR] $s";
    print STDERR "\n" unless ($flag);
}


#============Main function=====================#
sub _main{
  my $n_chr="";
  my $str="";
  my %seq;
  my %len;
  open IN1, "< $fasta";
  while(<IN1>){
    chomp($_);
    if($_=~/>/){
      if($str ne ""){
        $seq{$n_chr}=uc($str);
        $len{$n_chr}=length($str);
      }
      my $lenflag=index($_,"\s");
      if($_=~/chr/){
        $n_chr=substr($_,4,($lenflag-4));
      }else{
        $n_chr=substr($_,1,($lenflag-4));
      }
      $str="";
      next;
    }
    $str.=$_;
  }
  close IN1;
  
  if($str ne ""){
    $seq{$n_chr}=uc($str);
    $len{$n_chr}=length($str);
  }
  print "Loading fasta done\n";
  
  if((exists $seq{"M"}) and (!(exists $seq{"MT"}))){$seq{"MT"}=$seq{"M"};$len{"MT"}=$len{"M"};}
  if((exists $seq{"MT"}) and (!(exists $seq{"M"}))){$seq{"M"}=$seq{"MT"};$len{"M"}=$len{"MT"};}
  

  my $failed_str="";
  my $outstr="";
  my $lines=0;
  open IN1, "< $in_csv";
  while(<IN1>){
    my @temp=split(/,/,$_);
    if(@temp < 19){
    	next;
    }else{
    	$lines++;
    }
    
    if($lines==1){$outstr.="ID,Name,Strand,SNP,GenomeBuild,Chr,MapInfo,SourceStrand,SourceSNP,TOPSNP,lengthA,lengthB,mis_match,SNP_plus,type\n";next;}
  
    if(!exists $seq{$temp[9]}){$failed_str.=$temp[0].",".$temp[1].",".$temp[3].",".$temp[8].",".$temp[9].",".$temp[10]."\n";next;}
    
    my $pos1=index($temp[17],"[");
    my $pos2=index($temp[17],"]");
    my $pos1t=index($temp[16],"[");
    my $pos2t=index($temp[16],"]");
    
    my $lA=length($temp[5]);
    my $lB=length($temp[7]);

    if(($temp[3]=~/D\/I/) or ($temp[3]=~/I\/D/)){
    	if($temp[10] > $len{$temp[9]}){$failed_str.=$temp[0].",".$temp[1].",".$temp[3].",".$temp[8].",".$temp[9].",".$temp[10]."\n";next;}
      if($temp[10] < 0){$failed_str.=$temp[0].",".$temp[1].",".$temp[3].",".$temp[8].",".$temp[9].",".$temp[10]."\n";next;}
      
      my $strA=uc($temp[5]);
      my $strA_revcom=&reverse_complement($strA);
          
      my $strrefb=uc(substr($seq{$temp[9]},($temp[10]-$lA-1),$lA));
      my $strreff=uc(substr($seq{$temp[9]},$temp[10],$lA));
 
    	my $lstr=$pos2-$pos1-3;#deleted/inserted string length
    	my %refb;
      $refb{"b"}=uc(substr($seq{$temp[9]},($temp[10]-$lA-1),$lA));
      $refb{"b1"}=uc(substr($seq{$temp[9]},($temp[10]-2-$lA),$lA));
      for(my $i=3;$i<=($lstr+1);$i++){
      	$refb{"b".($i-1)}=uc(substr($seq{$temp[9]},($temp[10]-$i-$lA),$lA));
      }
      $refb{"bf1"}=uc(substr($seq{$temp[9]},($temp[10]-$lA),$lA));
      $refb{"bf2"}=uc(substr($seq{$temp[9]},($temp[10]-$lA+1),$lA));
      $refb{"bf3"}=uc(substr($seq{$temp[9]},($temp[10]-$lA+2),$lA));
      for(my $i=3;$i<$lstr;$i++){
      	$refb{"bf".($i+1)}=uc(substr($seq{$temp[9]},($temp[10]-$lA+$i),$lA));
      }
      
      my %reff;  
      $reff{"f"}=uc(substr($seq{$temp[9]},$temp[10],$lA));
      $reff{"f1"}=uc(substr($seq{$temp[9]},($temp[10]+1),$lA));
      for(my $i=2;$i<=$lstr;$i++){
        $reff{"f".$i}=uc(substr($seq{$temp[9]},($temp[10]+$i),$lA));
      }
      $reff{"fb1"}=uc(substr($seq{$temp[9]},($temp[10]-1),$lA));
      $reff{"fb2"}=uc(substr($seq{$temp[9]},($temp[10]-2),$lA));
      for(my $i=3;$i<=$lstr;$i++){
      	$reff{"fb".$i}=uc(substr($seq{$temp[9]},($temp[10]-$i),$lA));
      }
#      if($_=~/chr18_24239216_AAAC_INDEL/){print $strdownr2."\n"; exit; }
      
      my %count;
      foreach my $key (keys %refb){
      	$count{$key}=(( $refb{$key} ^ $strA ) =~ tr/\0//c);
      }
      foreach my $key (keys %reff){
      	$count{$key}=(( $reff{$key} ^ $strA_revcom ) =~ tr/\0//c);
      }
      
      my $min;
      my $type;
      foreach my $key (sort { $count{$a} <=> $count{$b} } keys %count){
      	$min=$count{$key};
      	$type=$key;
      	last;
      }
      
      my $flag_rev=0;
      my $min_sw=$lA;
      if($min>0){
      	my $sw1=(&sw($strA,$strrefb));
      	my $sw2=(&sw($strA_revcom,$strreff));
#        $min_sw=$lA - (((&sw($strA,$strrefb)) >= (&sw($strA_revcom,$strreff))) ? (&sw($strA,$strrefb)) : (&sw($strA_revcom,$strreff)));
        if($sw1<$sw2){
        	$flag_rev=1;
        	$min_sw=$lA -$sw2;
        }else{
        	$min_sw=$lA -$sw1;
        }
      }
      
      my $tempstr=substr($temp[17],($pos1+3),($pos2-$pos1-3));
      my $tempstrt=substr($temp[16],($pos1t+3),($pos2t-$pos1t-3));
      $outstr.=$temp[0].",".$temp[1].",".$temp[2].",".$temp[3].",".$temp[8].",".$temp[9].",".$temp[10].",".$temp[15].",".$tempstrt.",".$tempstr.",".$lA.",".$lB.",";
      if($min_sw<$min){
        $outstr.=$min_sw.",";
        $type="smith-waterman";
      }else{
        $outstr.=$min.",";
      }
      
      $outstr.="-,".$type."\n";
    }elsif($temp[3]=~/N\/A/){
      if($temp[10] > $len{$temp[9]}){$failed_str.=$temp[0].",".$temp[1].",".$temp[3].",".$temp[8].",".$temp[9].",".$temp[10]."\n";next;}
      if($temp[10] < 0) {$failed_str.=$temp[0].",".$temp[1].",".$temp[3].",".$temp[8].",".$temp[9].",".$temp[10]."\n";next;}
      
      my $temp5=substr($temp[5],0,($lA-1));
      my $strA=uc($temp5);
      my $strA_revcom=&reverse_complement($strA);
      my $strAt=uc($temp[5]);
      my $strAt_revcom=&reverse_complement($strAt);
      
      my $strrefb=uc(substr($seq{$temp[9]},($temp[10]-$lA-1),$lA));
      my $strreff=uc(substr($seq{$temp[9]},$temp[10],$lA));

      my %refb;
      $refb{"b"}=uc(substr($seq{$temp[9]},($temp[10]-$lA-1),$lA));
      $refb{"b1"}=uc(substr($seq{$temp[9]},($temp[10]-2-$lA),$lA));
      $refb{"b2"}=uc(substr($seq{$temp[9]},($temp[10]-3-$lA),$lA));
      $refb{"bf1"}=uc(substr($seq{$temp[9]},($temp[10]-$lA),$lA));
      $refb{"bf2"}=uc(substr($seq{$temp[9]},($temp[10]+1-$lA),$lA));
      $refb{"bf3"}=uc(substr($seq{$temp[9]},($temp[10]+2-$lA),$lA));
      
      my %reff;  
      $reff{"f"}=uc(substr($seq{$temp[9]},$temp[10],$lA));
      $reff{"f1"}=uc(substr($seq{$temp[9]},($temp[10]+1),$lA));
      $reff{"f2"}=uc(substr($seq{$temp[9]},($temp[10]+2),$lA));
      $reff{"fb1"}=uc(substr($seq{$temp[9]},($temp[10]-1),$lA));
      $reff{"fb2"}=uc(substr($seq{$temp[9]},($temp[10]-2),$lA));
      $reff{"fb3"}=uc(substr($seq{$temp[9]},($temp[10]-3),$lA));
      

      my %count;
      foreach my $key (keys %refb){
      	$count{$key}=(( $refb{$key} ^ $strAt ) =~ tr/\0//c);
      }
            
      foreach my $key (keys %reff){
      	$count{$key}=(( $reff{$key} ^ $strAt_revcom ) =~ tr/\0//c);
      }
      
      my $min;
      my $type;
      foreach my $key (sort { $count{$a} <=> $count{$b} } keys %count){
      	$min=$count{$key};
      	$type=$key;
      	last;
      }
      
      if(($type eq "b") and ($min > 0)){
      	my $reftemp=uc(substr($seq{$temp[9]},($temp[10]-$lA-1),($lA-1)));
        my $counttemp=(( $reftemp ^ $strA ) =~ tr/\0//c);
        if($counttemp < $min){$min=$counttemp;$type="b-s";}
      }
      
      if(($type eq "fb1") and ($min > 0)){
      	my $reftemp=uc(substr($seq{$temp[9]},$temp[10],($lA-1)));
        my $counttemp=(( $reftemp ^ $strA_revcom ) =~ tr/\0//c);
        if($counttemp < $min){$min=$counttemp;$type="fb1-s";}
      }

      my $refsnp=uc(substr($seq{$temp[9]},($temp[10]-1),1));
      if(($type eq "bf1") and ($min > 0)){
      	my $prosnp=uc(substr($strAt,($lA-2),1));
      	if($prosnp ne $refsnp){$min=($min-1);}
      }
      if(($type eq "bf2") and ($min > 0)){
      	my $prosnp=uc(substr($strAt,($lA-3),1));
      	if($prosnp ne $refsnp){$min=($min-1);}
      }
      if(($type eq "bf3") and ($min > 0)){
      	my $prosnp=uc(substr($strAt,($lA-4),1));
      	if($prosnp ne $refsnp){$min=($min-1);}
      }
      if(($type eq "fb2") and ($min > 0)){
      	my $prosnp=uc(substr($strAt,1,1));
      	if($prosnp ne $refsnp){$min=($min-1);}
      }
      if(($type eq "fb3") and ($min > 0)){
      	my $prosnp=uc(substr($strAt,2,1));
      	if($prosnp ne $refsnp){$min=($min-1);}
      }


      my $flag_rev=0;
      my $min_sw=$lA;
      if($min>0){
      	my $sw1=(&sw($strA,$strrefb));
      	my $sw2=(&sw($strA_revcom,$strreff));
        if($sw1<$sw2){
        	$flag_rev=1;
        	$min_sw=$lA -$sw2;
        }else{
        	$min_sw=$lA -$sw1;
        }
      }
     
      $outstr.=$temp[0].",".$temp[1].",".$temp[2].",".$temp[3].",".$temp[8].",".$temp[9].",".$temp[10].",".$temp[15].","."-".","."-".",".$lA.",".$lB.",";

      if($min_sw<$min){
        $outstr.=$min_sw.",";
        $type="smith-waterman";
      }else{
        $outstr.=$min.",";
      }

      my $top1;
      my $top12;
      my $top2;
      my $top22;
      my $mposi=int(length($temp[17])/2);
      if((length($temp[17])%2)==0){
      	$top1=substr($temp[17],($mposi-$lA),$lA);
      	$top12=substr($temp[17],($mposi-$lA+1),$lA-1);
      	$top2=substr($temp[17],$mposi,$lA);
      	$top22=substr($temp[17],$mposi,$lA-1);
      }else{
      	$top1=substr($temp[17],($mposi-$lA),$lA);
      	$top12=substr($temp[17],($mposi-$lA+1),$lA-1);
      	$top2=substr($temp[17],($mposi+1),$lA);
      	$top22=substr($temp[17],($mposi+1),$lA-1);
      }
      if(($type=~/^f/) or (($type eq "smith-waterman")&($flag_rev==1))){
      	if(($top1 eq $strAt) or ($top12 eq $strA)){
      	  $outstr.="RevTOP,";
      	}else{
      		$outstr.="SameTOP,";
      	}
      }else{
      	if(($top1 eq $strAt) or ($top12 eq $strA)){
      		$outstr.="SameTOP,";
      	}else{
      		$outstr.="RevTOP,";
      	}
      }
      $outstr.=$type."\n";
    }else{
    	if($temp[10] > $len{$temp[9]}){$failed_str.=$temp[0].",".$temp[1].",".$temp[3].",".$temp[8].",".$temp[9].",".$temp[10]."\n";next;}
      if($temp[10] < 0) {$failed_str.=$temp[0].",".$temp[1].",".$temp[3].",".$temp[8].",".$temp[9].",".$temp[10]."\n";next;}
      
      my $temp5=substr($temp[5],0,($lA-1));
      my $strA=uc($temp5);
      my $strA_revcom=&reverse_complement($strA);
      my $strAt=uc($temp[5]);
      my $strAt_revcom=&reverse_complement($strAt);
      
      my $strrefb=uc(substr($seq{$temp[9]},($temp[10]-$lA-1),$lA));
      my $strreff=uc(substr($seq{$temp[9]},$temp[10],$lA));

      my %refb;
      $refb{"b"}=uc(substr($seq{$temp[9]},($temp[10]-$lA-1),$lA));
      $refb{"b1"}=uc(substr($seq{$temp[9]},($temp[10]-2-$lA),$lA));
      $refb{"b2"}=uc(substr($seq{$temp[9]},($temp[10]-3-$lA),$lA));
      $refb{"bf1"}=uc(substr($seq{$temp[9]},($temp[10]-$lA),$lA));
      $refb{"bf2"}=uc(substr($seq{$temp[9]},($temp[10]+1-$lA),$lA));
      $refb{"bf3"}=uc(substr($seq{$temp[9]},($temp[10]+2-$lA),$lA));
      
      my %reff;  
      $reff{"f"}=uc(substr($seq{$temp[9]},$temp[10],$lA));
      $reff{"f1"}=uc(substr($seq{$temp[9]},($temp[10]+1),$lA));
      $reff{"f2"}=uc(substr($seq{$temp[9]},($temp[10]+2),$lA));
      $reff{"fb1"}=uc(substr($seq{$temp[9]},($temp[10]-1),$lA));
      $reff{"fb2"}=uc(substr($seq{$temp[9]},($temp[10]-2),$lA));
      $reff{"fb3"}=uc(substr($seq{$temp[9]},($temp[10]-3),$lA));
      
      my %count;
      foreach my $key (keys %refb){
      	$count{$key}=(( $refb{$key} ^ $strAt ) =~ tr/\0//c);
      }
            
      foreach my $key (keys %reff){
      	$count{$key}=(( $reff{$key} ^ $strAt_revcom ) =~ tr/\0//c);
      }
      
      
      my $min;
      my $type;
      foreach my $key (sort { $count{$a} <=> $count{$b} } keys %count){
      	$min=$count{$key};
      	$type=$key;
      	last;
      }
      
      if(($type eq "b") and ($min > 0)){
      	my $reftemp=uc(substr($seq{$temp[9]},($temp[10]-$lA-1),($lA-1)));
        my $counttemp=(( $reftemp ^ $strA ) =~ tr/\0//c);
        if($counttemp < $min){$min=$counttemp;$type="b-s";}
      }
      
      if(($type eq "fb1") and ($min > 0)){
      	my $reftemp=uc(substr($seq{$temp[9]},$temp[10],($lA-1)));
        my $counttemp=(( $reftemp ^ $strA_revcom ) =~ tr/\0//c);
        if($counttemp < $min){$min=$counttemp;$type="fb1-s";}
      }

      my $refsnp=uc(substr($seq{$temp[9]},($temp[10]-1),1));
      if(($type eq "bf1") and ($min > 0)){
      	my $prosnp=uc(substr($strAt,($lA-2),1));
      	if($prosnp ne $refsnp){$min=($min-1);}
      }
      if(($type eq "bf2") and ($min > 0)){
      	my $prosnp=uc(substr($strAt,($lA-3),1));
      	if($prosnp ne $refsnp){$min=($min-1);}
      }
      if(($type eq "bf3") and ($min > 0)){
      	my $prosnp=uc(substr($strAt,($lA-4),1));
      	if($prosnp ne $refsnp){$min=($min-1);}
      }
      if(($type eq "fb2") and ($min > 0)){
      	my $prosnp=uc(substr($strAt,1,1));
      	if($prosnp ne $refsnp){$min=($min-1);}
      }
      if(($type eq "fb3") and ($min > 0)){
      	my $prosnp=uc(substr($strAt,2,1));
      	if($prosnp ne $refsnp){$min=($min-1);}
      }


      my $flag_rev=0;
      my $min_sw=$lA;
      if($min>0){
      	my $sw1=(&sw($strA,$strrefb));
      	my $sw2=(&sw($strA_revcom,$strreff));
        if($sw1<$sw2){
        	$flag_rev=1;
        	$min_sw=$lA -$sw2;
        }else{
        	$min_sw=$lA -$sw1;
        }
      }
     
      my $tempstr1=substr($temp[17],($pos1+1),1);
      my $tempstr2=substr($temp[17],($pos2-1),1);
      my $tempstr1t=substr($temp[16],($pos1t+1),1);
      my $tempstr2t=substr($temp[16],($pos2t-1),1);
      $outstr.=$temp[0].",".$temp[1].",".$temp[2].",".$temp[3].",".$temp[8].",".$temp[9].",".$temp[10].",".$temp[15].",".$tempstr1t."\/".$tempstr2t.",".$tempstr1."\/".$tempstr2.",".$lA.",".$lB.",";

      if($min_sw<$min){
        $outstr.=$min_sw.",";
        $type="smith-waterman";
      }else{
        $outstr.=$min.",";
      }

      
      my $tempstr=$tempstr1."\/".$tempstr2;
      my $cstr1=&complement($tempstr1);
      my $cstr2=&complement($tempstr2);
      my $refcurpos=substr($seq{$temp[9]},($temp[10]-1),1);
      if(($type=~/^f/) or (($type eq "smith-waterman")&($flag_rev==1))){
      	if(substr($temp[3],1,3) eq $tempstr){
#      	  if(($cstr1 ne $refcurpos) & ($cstr2 ne $refcurpos)){
#      	  	$outstr.=$tempstr1."\/".$tempstr2.",";
#      	  }else{
      	  	$outstr.=$cstr1."\/".$cstr2.",";
#      	  }
      	}else{
#      		if(($tempstr1 ne $refcurpos) & ($tempstr2 ne $refcurpos)){
#      	  	$outstr.=$cstr1."\/".$cstr2.",";
#      	  }else{
      	  	$outstr.=$tempstr1."\/".$tempstr2.",";
#      	  }
      	}
      }else{
      	if(substr($temp[3],1,3) eq $tempstr){
#      	  if(($tempstr1 ne $refcurpos) & ($tempstr2 ne $refcurpos)){
#      	  	$outstr.=$cstr1."\/".$cstr2.",";
#      	  }else{
      	  	$outstr.=$tempstr1."\/".$tempstr2.",";
#      	  }
      	}else{
#      		if(($cstr1 ne $refcurpos) & ($cstr2 ne $refcurpos)){
#      	  	$outstr.=$tempstr1."\/".$tempstr2.",";
#      	  }else{
      	  	$outstr.=$cstr1."\/".$cstr2.",";
#      	  }
      	}
      }
      $outstr.=$type."\n";
    }
  }
  
  my $out_csv=$out_dir."new_".$out_name."\.csv";
  open OUT1, "> $out_csv";
  print OUT1 $outstr;
  close OUT1;

  my $out_fcsv=$out_dir."outdated_".$out_name."\.csv";
  open OUT2, "> $out_fcsv";
  print OUT2 $failed_str;
  close OUT2;
}


sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub complement {
        my $dna = shift;

	# complement the  DNA sequence
        $dna =~ tr/ACGTacgt/TGCAtgca/;
        return $dna;
}

sub sw{
  my ($seq1, $seq2) = @_;
  
  # scoring scheme
  my $MATCH    =  1; # +1 for letters that match
  my $MISMATCH = -1; # -1 for letters that mismatch
  my $GAP      = -1; # -1 for any gap
  
  # initialization
  my @matrix;
  $matrix[0][0]{score}   = 0;
  $matrix[0][0]{pointer} = "none";
  for(my $j = 1; $j <= length($seq1); $j++) {
      $matrix[0][$j]{score}   = 0;
      $matrix[0][$j]{pointer} = "none";
  }
  for (my $i = 1; $i <= length($seq2); $i++) {
      $matrix[$i][0]{score}   = 0;
      $matrix[$i][0]{pointer} = "none";
  }
  
  # fill
  my $max_i     = 0;
  my $max_j     = 0;
  my $max_score = 0;
  
  for(my $i = 1; $i <= length($seq2); $i++) {
      for(my $j = 1; $j <= length($seq1); $j++) {
          my ($diagonal_score, $left_score, $up_score);
          
          # calculate match score
          my $letter1 = substr($seq1, $j-1, 1);
          my $letter2 = substr($seq2, $i-1, 1);       
          if ($letter1 eq $letter2) {
              $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
          }
          else {
              $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
          }
          
          # calculate gap scores
          $up_score   = $matrix[$i-1][$j]{score} + $GAP;
          $left_score = $matrix[$i][$j-1]{score} + $GAP;
          
          if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
              $matrix[$i][$j]{score}   = 0;
              $matrix[$i][$j]{pointer} = "none";
              next; # terminate this iteration of the loop
          }
          
          # choose best score
          if ($diagonal_score >= $up_score) {
              if ($diagonal_score >= $left_score) {
                  $matrix[$i][$j]{score}   = $diagonal_score;
                  $matrix[$i][$j]{pointer} = "diagonal";
              }
              else {
                  $matrix[$i][$j]{score}   = $left_score;
                  $matrix[$i][$j]{pointer} = "left";
              }
          } else {
              if ($up_score >= $left_score) {
                  $matrix[$i][$j]{score}   = $up_score;
                  $matrix[$i][$j]{pointer} = "up";
              }
              else {
                  $matrix[$i][$j]{score}   = $left_score;
                  $matrix[$i][$j]{pointer} = "left";
              }
          }
          
          # set maximum score
          if ($matrix[$i][$j]{score} > $max_score) {
              $max_i     = $i;
              $max_j     = $j;
              $max_score = $matrix[$i][$j]{score};
          }
      }
  }
  
  # trace-back
  
  my $align1 = "";
  my $align2 = "";
  
  my $j = $max_j;
  my $i = $max_i;
  
  while (1) {
      last if $matrix[$i][$j]{pointer} eq "none";
      
      if ($matrix[$i][$j]{pointer} eq "diagonal") {
          $align1 .= substr($seq1, $j-1, 1);
          $align2 .= substr($seq2, $i-1, 1);
          $i--; $j--;
      }
      elsif ($matrix[$i][$j]{pointer} eq "left") {
          $align1 .= substr($seq1, $j-1, 1);
          $align2 .= "-";
          $j--;
      }
      elsif ($matrix[$i][$j]{pointer} eq "up") {
          $align1 .= "-";
          $align2 .= substr($seq2, $i-1, 1);
          $i--;
      }   
  }
  
#  $align1 = reverse $align1;
#  $align2 = reverse $align2;
#  print "$align1\n";
#  print "$align2\n";
  return $max_score;
}
