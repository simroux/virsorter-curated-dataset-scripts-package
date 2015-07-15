#!/usr/bin/perl
use strict;
use Parallel::ForkManager;
# Script to get gene sequences from phage sorter contigs and nett fasta file
# Argument 0 : toto
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[1])))
{
	print "# Script to get gene sequences from phage sorter contigs and nett fasta file
# Argument 0 : Directory of reference codon usage files (generated via cusp: http://emboss.bioinformatics.nl/cgi-bin/emboss/cusp)
# Argument 1 : Number of cpus to use
# Argument 2 : Output directory
# Argument 3 (not mandatory) : Number of phage to start from 
# Argument 4 (not mandatory) : Number of phage to stop at \n";
	die "\n";
}

my $path_to_cai="cai";
my $dir_in="Viral_metrics/Codon_usage/";
my $dir_host_ref=$ARGV[0];
my $n_processes = $ARGV[1];
my $out_dir=$ARGV[2];
my $starting_number=0;
if (defined($ARGV[3])){
	$starting_number=$ARGV[3];
}
my $stopping_number=-1;
if (defined($ARGV[4])){
	$stopping_number=$ARGV[4];
}
my $tmp_file_base=$out_dir."tmp_file_to_be_removed";
print "Listing references ..";
my @liste=<$dir_host_ref*_cusp.tab>;
my %check_host;
foreach(@liste){
	my $file=$_;
	$file=~/.*\/(.*)_cusp.tab/;
	$check_host{$1}=$file;
}
print "ok, references listed\n";
my @tab_host=sort keys %check_host;
my $n_hosts=$#tab_host+1;
print "$n_hosts hosts\n";
my @list_phage=<$dir_in*genes.fna>;
my $n_phage=$#list_phage+1;
if ($stopping_number==-1){$stopping_number=$n_phage+10;}
my $out="";

# Processing the files in parallel
my $pm = Parallel::ForkManager->new( $n_processes );
my $i=0;
foreach(@list_phage){
	$i++;
	my $file=$_;
	$file=~/.*\/(.*)_genes.fna/;
	my $id_phage=$1;
	my $out_file=$out_dir.$id_phage."_cai.tab";
	my $tmp_file=$tmp_file_base."_".$i.".tmp";
	if ($i<$starting_number){
		print "We are at $i and will not start until $starting_number -> we skip\n";
	}
	elsif($i>$stopping_number){
		print "We are at $i and stopped at $stopping_number -> we skip\n";
	}
	else{
		$pm->start and next;
		my $tag_done=0;
		if (-e $out_file){
			print "$out_file already here\n";
			open(TAB,"<$out_file") || die ("pblm opening file $out_file\n");
			my $n=0;
			while(<TAB>){
				chomp($_);
				if ($_=~/^>/){$n++;}
			}
			if ($n==$n_hosts){
				print "\t$out_file already done, the number of hosts $n_hosts match\n";
				$tag_done=1;
			}
			else{
				print "\t$out_file already here but looks incomplete so we redo\n";
				`rm $out_file`;
			}
				
		}
		if($tag_done==0){
			print "Processing $file -> $id_phage, which results will go in $out_file\n";
			for(my $i=0;$i<=$#tab_host;$i++){
				my $id_host=$tab_host[$i];
				my $file_host=$check_host{$tab_host[$i]};
				my $header='\>'.$id_host;
				$out=`echo $header >> $out_file`;
				my $cmd_line=$path_to_cai.' -seqall '.$file.' -cfile '.$file_host.' -outfile >(cat >> '.$out_file.') 2> /dev/null';
				$out=`bash -c '$cmd_line'`;
			}
		}
		$pm->finish;
	}
}
$pm->wait_all_children;


