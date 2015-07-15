#!/usr/bin/perl
use strict;
# Script to get the distances between all phage sequences and all hosts kmer frequencies
# Argument 0 : directory of host metrics
# Argument 1 : output directory (BEWARE : THE CONTENT OF THIS DIRECTORY WILL BE ALL REMOVED!!)
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[0])))
{
	print "# Script to get the distances between all phage sequences and all hosts kmer frequencies
# Argument 0 : directory of host metrics
# Argument 1 : output directory (BEWARE : THE CONTENT OF THIS DIRECTORY WILL BE ALL REMOVED!!)\n";
	die "\n";
}

my $dir_hosts=$ARGV[0]; # Directory with the host files listing their kmer frequency, generated with Jellyfish (http://www.cbcb.umd.edu/software/jellyfish/)
my $out_dir=$ARGV[1];
`rm $out_dir/*`;
my $dir_phages="Viral_metrics/Kmer_frequencies/";

print "We get the list of refs ... ";
my $liste_ref="list_refs";
open(LI,"<$liste_ref") || die ("pblm opening file $liste_ref\n");
my %check;
my @tab_s;
my @tab_di;
my @tab_tri;
my @tab_tetra;
while(<LI>){
	chomp($_);
	my @tab=split(",",$_);
	foreach(@tab){
		$check{$_}=1;
		if (length($_)==1){push(@tab_s,$_);}
		if (length($_)==2){push(@tab_di,$_);}
		if (length($_)==3){push(@tab_tri,$_);}
		if (length($_)==4){push(@tab_tetra,$_);}
	}
}
close LI;
print "done \n";
print "We get all values for the host ...";
my @tab_file_hosts=<$dir_hosts*_kmer_*.tab>;
my %store_vec_host;
my %check_host;
foreach(@tab_file_hosts){
	# We read the output file from Jellyfish, with the length of the mker indicated in the file name
	my $file=$_;
	$file=~/.*\/(.*)_kmer_(\d)\.tab/;
	my $id=$1;
	my $num=$2;
	# $num will be the kmer length
	$check_host{$id}=1;
	my %store;
	open(TAB,"<$file") || die ("pblm opening $file\n");
	while(<TAB>){
		chomp($_);
		my @tab=split(" ",$_);
		$store{"total"}+=$tab[1];
		if ($check{$tab[0]}==1){
			$store{$tab[0]}=$tab[1];
		}
		else{
			my $r=&revcomp($tab[0]);
			if ($check{$r}==1){
				$store{$r}=$tab[1];
			}
			else{
				print "!!!! $tab[0] / $r ??????\n";
			}
		}
	}
	close TAB;
	if ($store{"total"}==0){
		die("\n!!!!!!!!!!! Pblm, 0 counts in $file\n");
	}
	my @tab_temp;
	if ($num==1){
		foreach(@tab_s){
			push(@tab_temp,$store{$_}/$store{"total"});
		}
	}
	if ($num==2){
		foreach(@tab_di){
			push(@tab_temp,$store{$_}/$store{"total"});
		}
	}
	if ($num==3){
		foreach(@tab_tri){
			push(@tab_temp,$store{$_}/$store{"total"});
		}
	}
	if ($num==4){
		foreach(@tab_tetra){
			push(@tab_temp,$store{$_}/$store{"total"});
		}
	}
	$store_vec_host{$id}{$num}=join(",",@tab_temp);
}
print "done \n";
my @tab_hosts=keys %check_host;
my @tab_phages=<$dir_phages*_kmer_*.tab>;
foreach(@tab_phages){
	# We read the output file from Jellyfish, with the length of the mker indicated in the file name
	my $file=$_;
	$file=~/.*\/(.*)_kmer_(\d)\.tab/;
	my $id=$1;
	my $kmer=$2;
	print "Processing $id / $kmer-mer... ";
	my %store;
	open(TAB,"<$file") || die ("pblm opening file $file\n");
	while(<TAB>){
		chomp($_);
		my @tab=split(" ",$_);
		$store{"total"}+=$tab[1];
		if ($check{$tab[0]}==1){
			$store{$tab[0]}=$tab[1];
		}
		else{
			my $r=&revcomp($tab[0]);
			if ($check{$r}==1){
				$store{$r}=$tab[1];
			}
			else{
				print "!!!! $tab[0] / $r ??????\n";
			}
		}
	}
	close TAB;
	if ($store{"total"}==0){
		die("\n!!!!!!!!!!! Pblm, 0 counts in $file\n");
	}
	my @tab_temp;
	if ($kmer==1){
		foreach(@tab_s){
			push(@tab_temp,$store{$_}/$store{"total"});
		}
	}
	if ($kmer==2){
		foreach(@tab_di){
			push(@tab_temp,$store{$_}/$store{"total"});
		}
	}
	if ($kmer==3){
		foreach(@tab_tri){
			push(@tab_temp,$store{$_}/$store{"total"});
		}
	}
	if ($kmer==4){
		foreach(@tab_tetra){
			push(@tab_temp,$store{$_}/$store{"total"});
		}
	}
	my $vec_phage=join(",",@tab_temp);
	# We have the viral vector, we compare it with all the host vectors
	print " and now the distances ... ";
	my $out_file=$out_dir.$id."_dist.tab";
	open(S1,">>$out_file") || die ("pblm opening file $out_file\n");
	foreach(@tab_hosts){
		my $host=$_;
		my $dist=&get_distance($vec_phage,$store_vec_host{$host}{$kmer});
		if ($dist eq "-"){
			print "Pblm with $id vs $host ?\n";
		}
		print S1 "$id\t$host\t$kmer\t$dist\n";
	}
	close S1;
	print " done \n";
}


sub revcomp{
	my $seq=$_[0];
	$seq=~tr/atgc/tacg/;
	$seq=~tr/ATCG/TAGC/;
	$seq=reverse($seq);
	return $seq;
}


sub get_distance(){
	my $vec_1=$_[0];
	my $vec_2=$_[1];
	my @tab_1=split(",",$vec_1);
	my @tab_2=split(",",$vec_2);
	if (($#tab_1==-1) || ($#tab_2==-1)){
		print "########### pblm with distance $vec_1 vs $vec_2 \n";
		return "-";
	}
	else{
		my $dist;
		for(my $i=0;$i<=$#tab_1;$i++){
			my $d=abs($tab_1[$i]-$tab_2[$i]);
			$dist+=$d;
		}
		$dist/=($#tab_1+1);
		return $dist;
	}
}
