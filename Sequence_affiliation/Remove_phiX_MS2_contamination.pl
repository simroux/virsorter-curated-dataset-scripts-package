#!/usr/bin/perl
use strict;
# Script to check the detection of PhiX and MS2 in Refseq and WGS
# Argument 0 : fasta phage
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[0])))
{
	print "# Script to clean a prot file for Phage Sorter
# Argument 0 : input fasta file 
# Argument 1 : output of blastn against the potential contaminants\n";
	die "\n";
}


my $in_file=$ARGV[0];
my $blast_file=$ARGV[1];

my $out_fasta_file=$in_file;
$out_fasta_file=~s/\.f/_nett.f/;
if ($in_file eq $out_fasta_file){die("pblm with input fasta file $in_file\n");}

my $th_id=95;
my %check;
my $n=0;
# Reading the BLAST
open(BL,"<$out_file") || die ("pblm opening file $out_file\n");
while (<BL>){
	chomp($_);
	my @tab=split("\t",$_);
	if ($tab[2]>$th_id){
		print "We have a likely contamination -> $tab[0] / $tab[1] -> $tab[2]\n";
		print "\t$_\n";
		if (!defined($check{$tab[0]})){$n++;}
		$check{$tab[0]}=1;
	}
}
close BL;

if (!(-e $out_fasta_file)){
	print "Printing out file $out_fasta_file\n";
	open(FA,"<$in_file") || die ("pblm opening file $in_file\n");
	open(S1,">$out_fasta_file") || die ("pblm opening file $out_fasta_file\n");
	my $tag=0;
	while (<FA>){
		chomp($_);
		if ($_=~/^>(\S*)/){
			$tag=0;
			if ($check{$1}==1){print "we removed $1\n";}
			else{$tag=1;}
		}
		if($tag==1){print S1 "$_\n";}
	}
	close FA;
	close S1;
	print "$n sequences removed\n";
}