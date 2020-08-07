#!/usr/bin/env perl
use strict;
use Getopt::Long;
my @classes;
GetOptions(
	"class=s" => \@classes,
);
my %classes = map {$_ => 1;} @classes; #("Epidermis (H)(4)"=> 1, "Epidermis (H)(8)"=>1);
my $header = <>;
chomp $header;
my @headers = split /\t/, $header;
my @indices;
for (my $i=0; $i < @headers; $i++) {
	if ($classes{$headers[$i]}) {
		push @indices, $i;
	}
}
while (<>) {
	chomp;
	my @data = split /\t/;
	my $fix_up=0;
	my $fix_down=0;
	foreach my $i (@indices) {
		if (! ($data[$i] eq "NA")) {
			if ($data[$i] > 0.5) {
				$fix_up = 1;
			}
			elsif ($data[$i] < 0.5) {
				$fix_down = 1;
			}
		}
	}
	if ($fix_up && $fix_down) {
		print $_, "\n";
	}
}
