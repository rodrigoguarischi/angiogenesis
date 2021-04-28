#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my $basedir		= '/work/angiogenesis/data/RNA-seq/mRNA/poliA+';
my @batches		= (4..5,7..17);
my %sample_counts;

foreach my $c_batches ( @batches ){
	my $filepath = $basedir . "/batch_" . $c_batches;
	foreach my $c_file ( glob $filepath . "/*.gz" ){
		my $filename = basename($c_file);
		my ($sample, $barcode) = ($1,$2) if ( $filename =~ /^(.*)_([ACGT]{6})/ );
		my $cmd = "ln -s " . $c_file . " ./" . $sample . "-" . ++$sample_counts{$sample} . "_" . $barcode . ".fastq.gz";
		`$cmd`;
	}
}



