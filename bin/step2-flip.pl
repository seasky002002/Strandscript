#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;
use File::Path;
use Cwd 'abs_path';
use FindBin;
use lib $FindBin::Bin;

#========Variables defination===============#
#my $bin_dir=dirname($0);
#$bin_dir = abs_path($bin_dir) or die $!;
#my $dir=substr($bin_dir,0,rindex($bin_dir,"\/"));
my $ncsv_file=undef;
my $out_dir="./"; #default: current_dir
my $map_file=undef;
my $ped_file=undef;
my $cutoff=0.2;
my $out_map=undef;
my $out_ped=undef;
my $help=0;

#========Get options================#
unless (
    GetOptions(
        "in=s"    => \$ncsv_file,
        "ped=s"   => \$ped_file,
        "map=s"   => \$map_file,
        "o=s"     => \$out_dir,
        "c=f"     => \$cutoff,
        "h|help"  =>\$help
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
	  #0) in file checking  
    if ( !defined($ncsv_file) ) {
    	_error("input new manifest file (from step1) (-in) has not been defined yet\n");
    	_usage(1);
    }else{
    	if ( !-e $ncsv_file ) {
            _error("input new manifest file (-in) '$ncsv_file' does not exists\n");
            _usage(1);
      }
      $ncsv_file = abs_path($ncsv_file);
    }
    
    #1) ped file checking  
    if ( !defined($ped_file) ) {
    	_error("input ped file (-ped) has not been defined yet.\n");
    	_usage(1);
    }else{
    	if ( !-e $ped_file ) {
            _error("input ped file (-in) '$ped_file' does not exists.\n");
            _usage(1);
      }
      $ped_file = abs_path($ped_file);
    }
    
    #2) map file checking  
    if ( !defined($map_file) ) {
    	_error("input map file (-map) has not been defined yet.\n");
    	_usage(1);
    }else{
    	if ( !-e $map_file ) {
            _error("input map file (-in) '$map_file' does not exists.\n");
            _usage(1);
      }
      $map_file = abs_path($map_file);
    }

    #3) work directory checking
    if ( !defined($out_dir) ) {
    	_error("work directory (-o) has not been defined yet.\n");
    	_usage(1);
    }else{
    	if ( !-e $out_dir ) {
            _error("work directory (-o) '$out_dir' does not exists.\n");
            _usage(1);
      }
      $out_dir = abs_path($out_dir);
      if(!($out_dir=~/\/$/)){$out_dir.="\/";}
    }
    
    #3) cutoff checking
    if ( ($cutoff < 0) or ($cutoff > 1) ) {
    	_error("please set the cutoff (-c) to [0,1].\n");
    	_usage(1);
    }
}
    

#============Usage=====================#
sub _usage{
    my ($flag) = @_;
    my $usage=<<USAGE;
Usage: perl code_step2-TOP.pl [options] -in new_manifest.csv -map plink.map -ped plin.ped -o work_dir/
e.g: perl code_step2-TOP.pl -o test/ -in new_OncoChip.csv -map OncoChip.map -ped OncoChip.map -c 0.2
options:
-in [string]            required, input new manifest file in csv format from step1
-map [string]           required, plink map file for flipping
-ped [string]           required, plink ped file for flipping
-o [string]             output directory (default: current directory)
-c [real]               cutoff for mismatch percentage of Probe sequence, SNPs would be flipped which mismatch<cutoff (default: 0.2)
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
  my $flag=0;
  my %type;
  my %match;
  my %snp;

  open IN1, "< $ncsv_file";
  while(<IN1>){
  	$flag++;
  	if($flag==1){next;}
    chomp($_);
    my @temp=split(/,/,$_);
    
  
      if($temp[3]=~/I|D/){
      	$type{$temp[1]}="INDEL";
    #  	$snp{$temp[1]}=$temp[-1];
      }else{
      	if($temp[12] < ($cutoff*$temp[10])){
        	if($temp[3]=~/N/){
        		$type{$temp[1]}="SNP";
        		
        		if($temp[13] eq "SameTOP"){
        			$match{$temp[1]}="Y";
        		}else{
        			$match{$temp[1]}="C";
        		}
        		next;
        	}
          if($temp[9] eq $temp[13]){
          	$type{$temp[1]}="SNP";
    #      my $temp=substr($temp[3],1,3);
    #      if($temp eq $temp[14]){
            $match{$temp[1]}="Y";
          }else{
          	my $tempcom=&complement($temp[13]);
          	if($temp[9] eq $tempcom){
          		$type{$temp[1]}="SNP";
      	      $match{$temp[1]}="C";
      	    }else{
      	    	$match{$temp[1]}="N";
      	    }
          }
        }
    #    	$snp{$temp[1]}{"f"}=substr($temp[10],0,1);
    #    	$snp{$temp[1]}{"s"}=substr($temp[10],2,1);
      }
  #  }
  }
  close IN1;
  
  
  my $outstr1="";
  my $n=0;
  my %order;
  open IN2, "< $map_file";
  while(<IN2>){
    chomp($_);
    my @temp=split(/\s+/,$_);
    if(exists $type{$temp[1]}){
    	$outstr1.=$_."\n";
    }
#    }elsif(($temp[0] eq "0") or ($temp[0] eq "XY")){
#    	$outstr1.=$_."\n";
#    	$type{$temp[1]}="0XY";
#    }
    
    $order{$n}=$temp[1];
    $n++;
  }
  close IN2;
  
  
  my $outstr2="";
  open IN3, "< $ped_file";
  while(<IN3>){
    chomp($_);
    my @temp=split(/\s+/,$_);
    $outstr2.=join("\t",@temp[0..5]);
    for(my $i=6;$i<@temp;$i+=2){
      my $t=(($i-6)/2);
      if(exists $type{$order{$t}}){
#      	if(($type{$order{$t}} eq "INDEL") or ($type{$order{$t}} eq "0XY")){
      	if($type{$order{$t}} eq "INDEL"){
  #    		$outstr2.="\t".$snp{$order{$t}}; #need change
          $outstr2.="\t".$temp[$i]."\t".$temp[$i+1];
      	}else{
      		if($match{$order{$t}} eq "C"){
      			$outstr2.="\t".(&complement($temp[$i]))."\t".(&complement($temp[$i+1]));
      		}else{
            $outstr2.="\t".$temp[$i]."\t".$temp[$i+1];
          }
        }
      }
    }
    $outstr2.="\n";
  }
  close IN3;
  
  my $temp1=rindex($map_file,"\.map");
  my $temp2=rindex($map_file,"\/");
  my $name=substr($map_file,($temp2+1),($temp1-$temp2-1));
  my $out_map=$out_dir."flipped_".$name."\.map";
  my $out_ped=$out_dir."flipped_".$name."\.ped";
  open OUT1, "> $out_map";
  print OUT1 $outstr1;
  close OUT1;
  
  open OUT2, "> $out_ped";
  print OUT2 $outstr2;
  close OUT2;
}
  
sub complement {
        my $dna = shift;

	# complement the  DNA sequence
        $dna =~ tr/ACGTacgt/TGCAtgca/;
        return $dna;
}
