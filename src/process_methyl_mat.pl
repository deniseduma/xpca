#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);

#open IN, "<DNA.methyl.tum.noy.txt" or die $!;
open IN, "<bc_methyl_original_tum.txt" or die $!;
my @segs = <IN>;
print "size of segs is " . @segs . "\n";
close IN;

my $header = shift(@segs);
my @hfields = split('\t', $header);
shift(@hfields); shift(@hfields); shift(@hfields); shift(@hfields);
my @new_hfields = ();
foreach (@hfields) {
	my @subfields = split('-', $_);
	my $id = join('-', $subfields[0], $subfields[1], $subfields[2]);
	push(@new_hfields, $id);
}
my $new_header = join(' ', @new_hfields);

my %seg2var = ();

#map segments to their methylation variance accross patients
foreach my $line(@segs) {
	my @fields = split('\t', $line);
	my @aux = @fields[4..$#fields];
	my $key = join(' ', $fields[0], @aux);
	$seg2var{$key} = computeVar(\@aux);
	#print $line;
}

#sort map decreasingly by variance
#and retain top 2000 msot variable segements
open OUT, ">bc_methyl.txt" or die $!;
print OUT "$new_header\n";
my $counter = 0;
foreach my $seg(sort { $seg2var{$b} <=> $seg2var{$a} } keys %seg2var) {
	print OUT $seg;
	$counter++;
	if ($counter == 1000) {
		last;
	}
}
close OUT;

sub computeVar {
   my $data = $_[0];
   #if (@$data ==1) {
   #	return 0;
   #}
   my $mean = sum(@$data) / @$data;
   
   my $sqtotal = 0;
   foreach (@$data) {
	$sqtotal += ($_ - $mean) ** 2;
   }
   my $var = $sqtotal / @$data;
   
   return $var;
}

