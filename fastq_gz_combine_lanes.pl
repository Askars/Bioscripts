#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min max);
use IO::Compress::Gzip qw($GzipError);
use IPC::Open3;
use threads;

my $options = checkParams();

print STDERR Dumper($options);

my %combined;
opendir(my $dh, $options->{'dir'});
while (my $fn = readdir($dh)) {
    if ($fn =~ /(.*)_(L[0-9]+)_(R[12])_[0-9]+\.fastq\.gz/) {
         my $key = $1 . "_" . $3;
         push @{$combined{$key}}, $options->{'dir'} . "/" . $fn;
    }
}
closedir($dh);

my @threads;
foreach my $key (keys %combined) {
    while (scalar @threads > 9) {
        remove_finished_threads(\@threads);
        sleep(10);
    }
    my @sorted_files = sort {$a cmp $b} @{$combined{$key}};
    my $thr = threads->create('create_combined_file', $key, \@sorted_files);
    push @threads, $thr;
}

while (scalar @threads > 0) {
    remove_finished_threads(\@threads);
    sleep(10);
}


sub remove_finished_threads {
    my $threads_arr_ptr = shift;
    my @threads_arr = @{$threads_arr_ptr};
    my @new_threads_arr;
    foreach my $thr (@threads_arr) {
        if ($thr->is_joinable()) {
            $thr->join();
            next;
        }
        push @new_threads_arr, $thr;
    }
    @{$threads_arr_ptr} = @new_threads_arr;
}

sub create_combined_file {
    my $output_fn = shift;    
    my $files_to_combine_ptr = shift;
    my @files_to_combine = @{$files_to_combine_ptr};

    $output_fn .= "_001.fastq.gz";

    print join(" ", @files_to_combine) . " > " . "$output_fn\n";
    my $out_fh = new IO::Compress::Gzip $output_fn or die "gzip failed: $GzipError";
    my $buf;

    # Need to use IPC::Open3 instead of IO::Uncompress::Gzip as Illumina produces gzip files that causes perl to signal EOF early (before actual eof).
    foreach my $in_fn (@files_to_combine) {
        my($in, $in_fh, $err);
        my $pid = open3($in, $in_fh, $err, 'zcat', $in_fn);
        while(sysread($in_fh, $buf, 4096)) {
            syswrite($out_fh, $buf);
        }
        waitpid($pid, 0);
    }
    close($out_fh);
    print "Done: $output_fn\n";
}

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help+", "dir|d:s");
    my %options;

    # Add any other command line options, and the code to handle them
    #
    GetOptions( \%options, @standard_options );

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    if(!exists $options{'dir'} ) { printParamError ("Please specify a directory with the -d (--dir) flag."); }
    #if(!exists $options{'two'} ) { printParamError ("We need a second read to process!"); }
    
    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #
    my ($error) = @_;
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}


__DATA__

=head1 NAME

Combines the multilane files produced by NextSeqs and MiSeqs into single files in the current directory.

=head1 COPYRIGHT

=head1 DESCRIPTION

=head1 SYNOPSIS

Usage: fastq_gz_combine_lanes.pl -d <directory> 

=cut
