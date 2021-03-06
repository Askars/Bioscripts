#!/usr/bin/perl
use strict;
use warnings;
use threads;

use Getopt::Long;
use File::Basename;

use Bio::SeqIO;

my $options = check_params();

my $MAX_THREADS = 1; # Default to 1 thread.
my $phylosift_bin = "/srv/whitlam/bio/apps/sw/phylosift/phylosift_20120514/bin/phylosift";


if ($options->{'a'}) {
    $MAX_THREADS = $options->{'a'}
}

if (! -d $options->{'d'}) {
    die "Unable to find directory: " . $options->{'d'}. "\n";
}

# Calling Perl Threads with an open directory handle is asking for trouble. Place
# file names in an array and call from there.
my @files;
opendir(my $dh, $options->{'d'});
while (my $protein_file = readdir($dh)) {
    if (-f $options->{'d'} . "/$protein_file") {
        push @files, $protein_file;
    }
}
closedir($dh);

foreach my $protein_file (@files) {
    # Loop infinitely until there is a spare thread for processing.
    while (1) {
        if (threads->list(threads::all) < $MAX_THREADS) {
            threads->create('runPhylosift', ($options->{'d'} . "/$protein_file", $options->{'o'}));
            print STDERR $protein_file, "\n";
            last;
        } else {
            foreach my $thr (threads->list(threads::joinable)) {
                my $retcode = $thr->join();
            }
            sleep 10;
        }
    }
}

# Wait for all processing threads to complete.
while (threads->list(threads::all) > 0) {
    foreach my $thr (threads->list(threads::joinable)) {
        my $retcode = $thr->join();
    }
    sleep 10;
}

sub runPhylosift {
    my $protein_file = shift;
    my $output_prefix = shift;
    my($filename, $directories, $suffix) = fileparse($protein_file);
    my $output_folder = $filename ."_search";
    if (defined($output_prefix)) {
        $output_folder = "$output_prefix/$output_folder";
    }
    system("$phylosift_bin search --besthit --output=$output_folder $protein_file");
    system("$phylosift_bin align --besthit --output=$output_folder $protein_file");
    return 0;
}

################################################################################
# Subroutine: check_params()
# Handles command args via Getopt::Long and returns a reference to a hash of
# options.
################################################################################


sub check_params {
    my @standard_options = ( "help+", "man+");
    my %options;
    GetOptions( \%options, @standard_options, "d:s", "a:i", "o:s");
    exec("pod2usage $0") if $options{'help'};
    exec("perldoc $0")   if $options{'man'};
    exec("pod2usage $0") if (!( $options{'d'} && $options{'o'}));
    if (($options{'a'}) && ($options{'a'} < 1)) {
        die "Number of threads (-a) must be a positive integer. Got " . $options{'a'} . "\n";
    }
    return \%options;
}

__DATA__

=head1 NAME

    run_genomes.pl
   
=head1 DESCRIPTION

    Calls Phylosift on a directory containing protein fasta files, and produces
    multiple alignments.
    
=head1 SYNOPSIS

    run_genomes.pl -d <input_dir> -o <output_dir> [-a <num_threads>] 
        [-help] [-man]

        -d      Directory of protein fasta files (1 file per genome)
        -o      Directory to output results (defaults to .)
        -a      Number of threads to run.
        
=cut
