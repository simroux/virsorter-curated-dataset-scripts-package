#!/usr/bin/perl
use strict;
use Parallel::ForkManager;
# Script to parse the benchmark of MCL on different significativity thresholds and inflations
# Argument 0 : abc directory 
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[0])))
{
	print "# Script to parse the benchmark of MCL on different significativity thresholds and inflations
# Argument 0 : abc directory \n";
	die "\n";
}
my @liste_results=<$ARGV[0]/*.dump>;
my @liste_links=<$ARGV[0]/*.abc>;
my %store_neighbor;
my %check_neighbor;
my $file="$ARGV[0]/VirSorter_all_sig1.abc";
print "Parsing $file\n";
open(ABC,"<$file") || die ("pblm opening file $file\n");
while(<ABC>){
	chomp($_);
	my @tab=split("\t",$_);
	$store_neighbor{$tab[0]}{$tab[1]}=$tab[2];
	$store_neighbor{$tab[1]}{$tab[0]}=$tab[2];
	$check_neighbor{"1"}{$tab[0]}{$tab[1]}=1;
	$check_neighbor{"1"}{$tab[1]}{$tab[0]}=1;
	if ($tab[2]>=5){
		$check_neighbor{"5"}{$tab[0]}{$tab[1]}=1;
		$check_neighbor{"5"}{$tab[1]}{$tab[0]}=1;
	}
	if ($tab[2]>=10){
		$check_neighbor{"10"}{$tab[0]}{$tab[1]}=1;
		$check_neighbor{"10"}{$tab[1]}{$tab[0]}=1;
	}
	if ($tab[2]>=15){
		$check_neighbor{"15"}{$tab[0]}{$tab[1]}=1;
		$check_neighbor{"15"}{$tab[1]}{$tab[0]}=1;
	}
	if ($tab[2]>=20){
		$check_neighbor{"20"}{$tab[0]}{$tab[1]}=1;
		$check_neighbor{"20"}{$tab[1]}{$tab[0]}=1;
	}
	if ($tab[2]>=30){
		$check_neighbor{"30"}{$tab[0]}{$tab[1]}=1;
		$check_neighbor{"30"}{$tab[1]}{$tab[0]}=1;
	}
	if ($tab[2]>=50){
		$check_neighbor{"50"}{$tab[0]}{$tab[1]}=1;
		$check_neighbor{"50"}{$tab[1]}{$tab[0]}=1;
	}
}
close ABC;

# Reading the taxonomy of RefSeq file from the cytoscape file that list all RefSeq sequences and their Genus affiliation (as column 5)
my $taxo_file="PhageSorter_all_sig1_taxo-cyto.txt";
my %infos;
my %check;
open(TAX,"<$taxo_file") || die ("pblm opening file $taxo_file\n");
while(<TAX>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[5] ne "unclassified"){
		$infos{$tab[0]}{"genus"}=$tab[5];
		$check{$tab[5]}++;
	}
}
close TAX;
my @tab_genus=keys %check;

my $n_processes = 4;
my $pm = Parallel::ForkManager->new( $n_processes );

foreach(@liste_results){
	my $file=$_;
	$pm->start and next;
	$file=~/PhageSorter_all_sig(\d+)_I(.+)\.dump/;
	my $out_file=$file."parse_results";
	my $sig=$1;
	my $inflation=$2;
	print "Parsing $file -> $sig and $inflation\n";
	my %counts;
	my %counts_bis;
	my $n_cluster=0;
	my $i=0;
	my $n_icc=0;
	my $store_icc=0;
	my %cluster;
	open(DMP,"<$file") || die ("pblm opening file $file\n");
	while(<DMP>){
		chomp($_);
		my @tab=split("\t",$_);
		my $cluster_id="Cluster_".$i;
	# 	print "Cluster $cluster_id\n";
		if ($#tab>0){
			$n_cluster++;
			for (my $j=0;$j<=$#tab;$j++){
				$n_icc++;
				$cluster{$tab[$j]}=$cluster_id;
				if ($tab[$j]=~/.*_gi_.*/){}
				elsif($tab[$j]=~/NC_.*/){
					if (defined($infos{$tab[$j]})){
						$counts{$cluster_id}{$infos{$tab[$j]}{"genus"}}++;
						$counts_bis{$infos{$tab[$j]}{"genus"}}{$cluster_id}++;	
					}
				}
			}		
			# For the ICC
			if ($#tab>2){
				for (my $j=0;$j<=$#tab;$j++){
					my $n_neigh=0;
					my $n_actual=0;
					my %relevant_neigh;
					for (my $k=0;$k<=$#tab;$k++){
						if ($k!=$j){
							if (defined($store_neighbor{$tab[$j]}{$tab[$k]}) && $store_neighbor{$tab[$j]}{$tab[$k]}>=$sig){
								my @tab_rel=keys %relevant_neigh;
								foreach(keys %relevant_neigh){
									$n_actual+=$check_neighbor{$sig}{$tab[$k]}{$_};
								}
								$relevant_neigh{$tab[$k]}=1;
								$n_neigh++;
							}
						}
					}
					if ($n_neigh<=1){
					}
					else{
						my @tab_rel=keys %relevant_neigh;
						my $n_rel=$#tab_rel+1;
						my $n_pot=$n_rel*($n_rel-1)/2;
						my $coeff=$n_actual/$n_pot;
						$store_icc+=$coeff;
						print "cluster $cluster_id, node $tab[$j] -> $n_rel relevant neighbors, with $n_actual connections on $n_pot potential -> $coeff coefficient\n";
					}
				}
			}
			
		}
		$i++;
	}
	close DMP;
	my $avg=0;
	my $n=0;
	foreach(@tab_genus){
		my @k=keys %{$counts_bis{$_}};
		$avg+=$#k+1;
		$n++;
	}
	$avg/=$n;
# 	print "\tAverage -> $sig $inflation $avg ($n)\n";
	my $avg_2=0;
	my $n_2=0;
	foreach(keys %counts){
		my @k=keys %{$counts{$_}};
		$avg_2+=$#k+1;
		$n_2++;
	}
	$avg_2/=$n_2;
# 	print "\tAverage -> $sig $inflation $avg_2 ($n_2)\n";
	my $f_measure=2*($avg*$avg_2)/($avg+$avg_2);
	$store_icc/=$n_icc;
	open(S1,">$out_file") || die ("pblm opening file $out_file\n");
	print S1 "Significativity th\tInflation th\tNb clusters\tF Measure Refseq Genus\tICC\n";
	print S1 "$sig\t$inflation\t$n_cluster\t$f_measure\t$store_icc\n";
	close S1;
	$pm->finish;
}
$pm->wait_all_children;