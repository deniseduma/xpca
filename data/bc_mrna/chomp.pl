#!/usr/bin/perl

use strict;
use warnings;

open IN, "<bc_mRNA_tum_small_net_firehose.txt.back" or die $!;
open OUT, ">bc_mRNA_tum_small_net_firehose.txt" or die $!;
while (<IN>) {
	chomp($_);
	print OUT "$_\n";
}

close IN;
close OUT;
