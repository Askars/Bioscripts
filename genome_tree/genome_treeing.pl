#!/usr/bin/perl
use strict;
use warnings;
use threads;

use Getopt::Long;
use File::Basename;
use Data::Dumper;

use Bio::SeqIO;
use Bio::Seq;

my $options = check_params();

if (! -d $options->{'d'}) {
    die "Unable to find directory: " . $options->{'d'}. "\n";
}

my $markers_by_genome = {};
my $markers_by_marker = {};

my $marker_lengths = {};

opendir(my $dh, $options->{'d'});
while (my $genome_dir = readdir($dh)) {
    my $align_dir = $options->{'d'}. "/$genome_dir/alignDir";
    if (-d $align_dir) {
        opendir(my $align_dh, $align_dir);
        while (my $fasta_file = readdir($align_dh)) {
            if ($fasta_file =~ /^(PMPROK.*).fasta$/ ) {
                my $prefix = $1;
                my $fasta_fh = Bio::SeqIO->new(-format => 'fasta',
                                               -file => "$align_dir/$fasta_file");
                my $seq_obj = $fasta_fh->next_seq();
                if (! $seq_obj) {
                    next;
                }
                $markers_by_genome->{$genome_dir}->{$prefix} = $seq_obj->seq();
                $markers_by_marker->{$prefix}->{$genome_dir} = $seq_obj->seq();
                if ($fasta_fh->next_seq()) {
                    warn "$genome_dir has more than 1 candidate for $prefix, ",
                    "ignoring all but the first\n";
                }
            }
        }
    }
}

# Check that all the markers are the same length.
foreach my $marker_name (keys %{$markers_by_marker}) {
    my @genome_names = keys %{$markers_by_marker->{$marker_name}};
    my $first = pop @genome_names;
    my $first_length = length($markers_by_marker->{$marker_name}->{$first});
    while (my $current = pop @genome_names) {
        my $current_length =  length($markers_by_marker->{$marker_name}->{$current});
        if ($first_length != $current_length) {
            die "Aligned sequence length for $first and $current are discordant ".
                "for marker $marker_name\n. Got $first_length and $current_length respectively.\n" 
        }
    }
    $marker_lengths->{$marker_name} = $first_length;
}


my $fasta_out = Bio::SeqIO->new(-format => 'fasta');
foreach my $genome_name (sort {$a cmp $b} keys %{$markers_by_genome}) {
    my $concat_seq = '';
    foreach my $marker_name (sort {$a cmp $b} keys %{$markers_by_marker}) {
        my $seq = $markers_by_genome->{$genome_name}->{$marker_name};
        if (defined($seq)) {
            $concat_seq .= $seq
        } else {
            $concat_seq .= "-" x $marker_lengths->{$marker_name};
        }
    }
    my $seq_obj = Bio::Seq->new(-alphabet => 'protein',
                                -id => $genome_name,
                                -seq => $concat_seq);
    $fasta_out->write_seq($seq_obj);
}

################################################################################
# Subroutine: check_params()
# Handles command args via Getopt::Long and returns a reference to a hash of
# options.
################################################################################

sub check_params {
    my @standard_options = ( "help+", "man+");
    my %options;
    GetOptions( \%options, @standard_options, "d:s");
    exec("pod2usage $0") if $options{'help'};
    exec("perldoc $0")   if $options{'man'};
    exec("pod2usage $0") if (!( $options{'d'}));
    return \%options;
}

__DATA__

=head1 NAME

    aligned_fasta_subseq.pl
   
=head1 DESCRIPTION

    Given an aligned fasta file (a FASTA file where all entries are aligned
    against each other) and trim points, this script will create an aligned
    FASTA file. All entries will be trimmed based on the absolute positions of
    the sequences in the alignment, or relative to the (unaligned) nucleotide
    positions of a specific entry in the file if an ID is specified.

=head1 SYNOPSIS

    aligned_fasta_subseq.pl -f <aligned_fasta> [-tr int[,int]] [-id <fasta_id>]
        [-help] [-man]

        -f      Input aligned FASTA.
        -tr     Trim points (either 'start' or 'start, stop'). Relative if -id
                is specified, absolute otherwise.
        -id     Identifier for relative position.
=cut
