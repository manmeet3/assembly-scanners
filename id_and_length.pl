#!/usr/bin/perl

use Bio::SeqIO;
use Getopt::Std;
use vars qw( $opt_f );

my $bio_io = &init();

while( $seq_o = $bio_io->next_seq() ) {
    printf("%s\t%d\n", $seq_o->display_id, $seq_o->length );
}

sub init {
    getopts( 'f:' );
    unless( -f $opt_f ) {
	print( "id_and_length.pl -f <fasta file>\n" );
	print( "Makes an output file with ID and length of each sequence.\n" );
	exit( 0 );
    }

    return Bio::SeqIO->new( '-file' => $opt_f,
			    '-format' => 'fasta' );
}
