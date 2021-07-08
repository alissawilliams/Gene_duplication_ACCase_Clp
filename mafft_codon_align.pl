#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: $0 input_fasta output_fasta\n\n";

my $input = shift or die ($usage);
my $output = shift or die ($usage);

align_fasta_with_mafft_codon($input, $output);