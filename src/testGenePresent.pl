#!/usr/bin/perl

use strict;
use warnings;

die "Input directory (bc/ov/ucec/luad/gbm) required!" unless (@ARGV==1);

my $wd = $ARGV[0];

#open IN, "<BIOGRID.txt" or die $!;
open IN, "</home/duma/TCGA/data/STRING.txt" or die $!;
#open GENES, "<bc_mRNA_genes_top1000.txt" or die $!;
#open OUT, ">bc_mRNA_genes_top866.txt" or die $!;
open GENES, "</home/duma/TCGA/data/$wd/${wd}_genes.txt" or die $!;
open OUT, ">/home/duma/TCGA/data/$wd/${wd}_genes_net.txt" or die $!;

my %uniques = ();
while (<IN>) {
	chomp;
	my @fields = split('\t');
	die "Not exactly two fields in @fields!\n" unless (@fields==2);
	$uniques{$fields[0]} = $fields[0];
	$uniques{$fields[1]} = $fields[1];
}
print "Number of unique genes " . scalar(keys %uniques) . "\n";

my $total = 0; 
my $count = 0;
while(<GENES>) {
	chomp;
	$total++;
	if (exists($uniques{$_})) {
		$count++;
		print OUT "$_\n";
	}
}
print "Number of genes found is $count out of $total\n";

close IN;
close GENES;
close OUT;
