#!/usr/bin/perl
use strict;
# Script to affiliate the contigs from PhageSorter
# Argument 0 : Fasta file of the proteins
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[3])))
{
	print "# Script to affiliate the contigs from PhageSorter
# Argument 0 : Fasta file of the proteins
# Argument 1 : Affiliation file of the protein clusters
# Argument 2 : Affiliation file of the contigs
# Argument 3 : out file\n";
	die "\n";
}


my $prot_file=$ARGV[0];
my $cluster_affi=$ARGV[1]; # file Affiliation_Phage_protein_clusters-with-virome.tab
my $affi_file=$ARGV[2]; # (output from VirSorter)
my $out_file=$ARGV[3];

# Loading the affiliation of the protein clusters
my %affi;
open(AFI,"<$cluster_affi") || die ("pblm opening file $cluster_affi\n");
while(<AFI>){
	chomp($_);
	my @tab=split("\t",$_);
	my @tab_1=split(";",$tab[1]);
	foreach(@tab_1){
		my @t=split(":",$_);
		my @t2=split(",",$t[1]);
		if ($#t2==0){
			my @t3=split(/\|/,$t2[0]);
			$affi{$tab[0]}{$t[0]}=$t3[0];
		}
		else{
			$affi{$tab[0]}{$t[0]}="unclassified";
			my $total=0;
			foreach(@t2){
				my @t3=split(/\|/,$_);
				$total+=$t3[1];
			}
			my $th=$total/2; # If more than one affiliation, majority rules
			foreach(@t2){
				my @t3=split(/\|/,$_);
				if ($t3[1]>$th){
					$affi{$tab[0]}{$t[0]}=$t3[0];
				}
			}
		}
	}
}
close AFI;


# Listing all genes predicted for each genome
my %check;
my %count;
my %cross_prot_pred;
open(PROT,"<$prot_file") || die ("pblm opening file $prot_file");
while(<PROT>){
	chomp($_);
	if ($_=~/^>(.*)-(gene_\d+)/){
		my $id_genome=$1;
		my $id_gene=$2;
		my $id_prot=$id_genome."-".$id_gene;
		$check{$id_prot}=1;
		$count{$id_genome}++;
		$cross_prot_pred{$id_prot}{$id_genome}=1;
	}
	elsif($-=~/^>/){
		print "$_ ??????\n";
	}
}
close PROT;

my @tab_taxo=("superkingdom","viralkingdom","order","family","genus","species");
my $th_score=50;
my %hit;

# Reading the affiliation file (output from VirSorter)
my $affi_file=$_;
my $c_seq=$1;
open(AFI,"<$affi_file") || die ("pblm opening file $affi_file\n");
while(<AFI>){
	chomp($_);
	if ($_=~/^>.*/){}
	else{
		my @tab=split(/\|/,$_);
		if ($tab[6] ne "-"){
			my $seq=$tab[0];
			if ($check{$seq}==1){
				if ($tab[6]>$th_score){
					my $cluster=$tab[5];
					foreach(keys %{$cross_prot_pred{$seq}}){
						print "hit for $_\n";
						$hit{$_}{$cluster}=$tab[6];
					}
				}
			}
		}
	}
}
close AFI;


my @tab_taxo=("superkingdom","viralkingdom","order","family","genus","species");
my $n_hits=3; ## We perform a LCA on the 3 best hits if they are above 50 and if they are within 75 % of the best BLAST hit score
open(S1,">$out_file") || die ("pblm opening file $out_file\n");
foreach(keys %count){
	my $id_sequence=$_;
	print "$id_sequence\n";
	if (!defined($hit{$id_sequence})){
		print "No hits for $id_sequence\n";
		print S1 "$id_sequence\t$count{$id_sequence}\tunclassified\n";
	}
	else{
		my @tab=sort { $hit{$id_sequence}{$b} <=> $hit{$id_sequence}{$a} } keys %{$hit{$id_sequence}};
		my %merge;
		print "\t Hit 0 on $tab[0] -> score $hit{$id_sequence}{$tab[0]} => Affi ";
		foreach(@tab_taxo){
			my $lvl=$_;
			$merge{$lvl}=$affi{$tab[0]}{$lvl};
			print "$lvl : $affi{$tab[0]}{$lvl} ";
		}
		print "\n";
		my $th_score=0.75*$hit{$id_sequence}{$tab[0]};
		for (my $i=1;$i<=$#tab;$i++){
			if (defined($tab[$i]) && ($hit{$id_sequence}{$tab[$i]}>$th_score)){
				print "\t Hit $i on $tab[$i] -> score $hit{$id_sequence}{$tab[$i]} => Affi ";
				foreach(@tab_taxo){
					print "$_ : $affi{$tab[$i]}{$_} ";
					if ($affi{$tab[$i]}{$_} eq $affi{$tab[0]}{$_}){}
					else{$merge{$_}="unclassified";print " ### ";}
				}
				print "\n";
			}
		}
		my $line="";
		foreach(@tab_taxo){
			$line.=$merge{$_}.",";
		}
		chop($line);
		print S1 "$id_sequence\t$count{$id_sequence}\t$line\n";
	}
}
close S1;