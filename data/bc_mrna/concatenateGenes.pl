#!/usr/bin/perl

use strict;
use warnings;

open IN, "<bc_genes.txt" or die $!;
open OUT, ">bc_genes_header.txt" or die $!;

my @in=<IN>;
my $first = shift(@in);
chomp($first);
print OUT "$first";

foreach my $i(@in) {
	chomp($i);
	print OUT " $i";
}

close IN;
close OUT;
