#!/usr/bin/perl
use strict;
# Script to predict host based on tetranucleotide distance
# Argument 0 : level of host similar to the actual host to exclude from the database
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[0])))
{
	print "# Script to predict host based on tetranucleotide distance
# Argument 0 : level of host similar to the actual host to exclude from the database (empty, species or genus)\n";
	die "\n";
}

my $out_file="host_predictions.tab";
my $distance_dir="Kmer_distance/"; # Directory with one file per viral genome, listing all the distances to all host genomes in the database in the form Viral Genome / Host Genome / Kmer length / Distance (tab-separated)

my $lvl_to_exclude="";
if ($ARGV[0] eq "species"){$lvl_to_exclude="species";}
if ($ARGV[0] eq "genus"){$lvl_to_exclude="genus";}

my @liste_files=<$distance_dir*.tab>; # listing the distance files (one per virus)

# Take the information on the actual host from the sequence summary file
my $summmary_file="Summary.tab";
my %host;
my @tab_taxo_host=("superkingdom","phylum","class","order","family","genus");
open(TAB,"<$summmary_file") || die ("pblm opening file $summmary_file\n");
while(<TAB>){
	chomp($_);
	my @tab=split("\t",$_);
	$tab[5]=&abrege($tab[5]);
	$host{$tab[0]}=$tab[5];
}
close TAB;


# Take taxonomy information on the host genomes from RefSeq
my $microbes="Refseq_14-03-28_v2.lst";
my %taxo_host;
open(TAB,"<$microbes") || die ("pblm opening file $microbes\n");
while(<TAB>){
	chomp($_);
	my @tab=split("\t",$_);
	$tab[13]=~s/\s/_/g;
	$tab[13]=~s/,/_/g;
	$tab[13]=~s/-/_/g;
	$tab[13]=~s/\(/_/g;
	$tab[13]=~s/\)/_/g;
	$tab[13]=~s/\[/_/g;
	$tab[13]=~s/\]/_/g;
	$tab[13]=~s/\//_/g;
	$tab[13]=~s/\\/_/g;
	$tab[13]=~s/:/_/g;
	$tab[13]=~s/;/_/g;
	$tab[13]=~s/=/_/g;
	$tab[13]=~s/'//g;
	$taxo_host{$tab[13]}{"class"}=$tab[6];
	$taxo_host{$tab[13]}{"order"}=$tab[7];
	$taxo_host{$tab[13]}{"family"}=$tab[8];
	$taxo_host{$tab[13]}{"genus"}=$tab[9];
	$taxo_host{$tab[13]}{"species"}=$tab[10];
	$taxo_host{$tab[13]}{"subspecies"}=$tab[11];
}
close TAB;


my %metrics;
my %count_by_dist;
my @tab_taxo_host=("class","order","family","genus");
open(S1,">$out_file") || die ("pblm opening file $out_file\n");
print S1 "Accession\tHost.class\tQuality_pred_class\tQuality_pred_order\tQuality_pred_family\tQuality_pred_genus\tDist\n";
my $n_panic=0;
foreach(@liste_files){
	my $dist_file=$_;
	$dist_file=~/.*\/([^\/]*)_dist\.tab/;
	my $accession=$1;
	my @t2=split("_",$accession);
	my $class=$t2[0];
	my $line=$accession."\t".$class;
	# Checking which host genome has the minimal distance with the viral genome
	my $min=999999;
	my $pred_host="unknown";
	my %store_score;
	my %count_family;
	my $pred_family="";
	# Reading all the distance to hosts for this virus
	open(TAB,"<$dist_file") || die ("pblm opening file $dist_file\n");
	while(<TAB>){
		chomp($_);
		my @tab=split("\t",$_);
		# If we are on tetranucleotide distances
		if ($tab[2]==4){
			# If we are below the minimum and this is not a completely unknown host
			if ($tab[3]<$min && $taxo_host{$tab[1]}{"family"} ne "unknown"){
				# If this is not a host that we chose not to include in the database
				if ($taxo_host{$tab[1]}{$lvl_excluded} ne $taxo_host{$true_host}{$lvl_excluded}){
					# Then we change the minimum and the prediction, as this is our new predicted host
					$min=$tab[3];
					$pred_host=$tab[1];
				}
			}
		}
	}
	close TAB;
	print "$accession\t$taxo_host{$host{$accession}}{family}\t$taxo_host{$pred_host}{family}\t$min\n";
	my $dist_category;
	# Sorting the prediction in different categories based on the distance
	if ($min<0.0004){$dist_category="1";}
	elsif ($min<0.001){$dist_category="2";}
	elsif ($min<0.01){$dist_category="3";}
	else{$dist_category="4";}
	foreach(@tab_taxo_host){
		if ($taxo_host{$pred_host}{$_} eq $taxo_host{$host{$accession}}{$_}){
			$metrics{$_}{"total"}++;
			$metrics{$_}{"good"}++;
			$line.="\tgood";
			$count_by_dist{$dist_category}{$_}{"total"}++;
			$count_by_dist{$dist_category}{$_}{"good"}++;
		}
		elsif($taxo_host{$host{$accession}}{$_} ne "unknown" && $taxo_host{$host{$accession}}{$_} ne ""){
			$metrics{$_}{"total"}++;
			$count_by_dist{$dist_category}{$_}{"total"}++;
			$line.="\tbad";
		}
		else{
			$line.="\tunconclusive";
		}
	}
	$line.="\t".$min;
	print S1 "$line\n";
}
close S1;

foreach(@tab_taxo_host){
	print "$metrics{$_}{good} good family on $metrics{$_}{total} $_ predictable\n";
}

foreach(sort keys %count_by_dist){
	print "class\t$_\t$count_by_dist{$_}{class}{total}\t$count_by_dist{$_}{class}{good}\t$count_by_dist{$_}{class}{bad}\n";
	print "order\t$_\t$count_by_dist{$_}{order}{total}\t$count_by_dist{$_}{order}{good}\t$count_by_dist{$_}{order}{bad}\n";
	print "family\t$_\t$count_by_dist{$_}{family}{total}\t$count_by_dist{$_}{family}{good}\t$count_by_dist{$_}{family}{bad}\n";
	print "genus\t$_\t$count_by_dist{$_}{genus}{total}\t$count_by_dist{$_}{genus}{good}\t$count_by_dist{$_}{genus}{bad}\n";
}


sub abrege(){
	my $word=$_[0];
	$word=~s/\s/_/g;
	$word=~s/,/_/g;
	$word=~s/-/_/g;
	$word=~s/\(/_/g;
	$word=~s/\)/_/g;
	$word=~s/\[/_/g;
	$word=~s/\]/_/g;
	$word=~s/\//_/g;
	$word=~s/\\/_/g;
	$word=~s/:/_/g;
	$word=~s/;/_/g;
	$word=~s/=/_/g;
	$word=~s/'//g;
	return($word);
}
