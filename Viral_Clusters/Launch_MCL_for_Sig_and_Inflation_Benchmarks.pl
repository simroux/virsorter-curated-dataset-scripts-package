#!/usr/bin/perl
use strict;
# Script to perform a mcl on multiple files with different significativity thresholds (i.e. different abc files)
# Argument 0 : ABC file directory
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[0])))
{
	print "# Script to perform a mcl on multiple files with different significativity thresholds (i.e. different abc files)
# Argument 0 : ABC file directory\n";
	die "\n";
}

my @abc_files=<$ARGV[0]/*.abc>;
foreach(@abc_files){
	my $abc_file=$_;
	my $path_to_mcl_bin="../Utiles/Mcl_12_135/mcl-12-135/bin/";
	my $base_file=$abc_file;
	$base_file=~s/\.abc//g;
	my $out_mci=$base_file.".mci";
	my $out_tab=$base_file.".tab";
	my $cmd_mcxload="$path_to_mcl_bin/mcxload -abc $abc_file --stream-mirror -o $out_mci -write-tab $out_tab";
	print "$cmd_mcxload\n";
	my $out=`$cmd_mcxload`;
	print "Mxc Load : $out\n";
	for(my $i=1.5;$i<=5;$i+=0.25){
		my $dump_file=$base_file."_I".$i.".dump";
		my $cmd_mcl="$path_to_mcl_bin/mcl $out_mci -I $i -use-tab $out_tab -o $dump_file";
		print "$cmd_mcl\n";
		$out=`$cmd_mcl`;
		print "Mcl : $out\n";
	}
}