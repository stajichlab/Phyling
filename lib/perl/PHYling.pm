package PHYling;

use strict;
use warnings;

use vars qw(@ISA @EXPORT @EXPORT_OK);
require Exporter;
@ISA = qw(Exporter);

our $DEBUG = 0;
our %DECOMPRESS = ( 'gz' => 'zcat',
		   'bz2' => 'bzcat' );

our %COMPRESS = ( 'gz' => 'gzip -c',
		 'bz2' => 'bzip2 -c' );


@EXPORT = qw(&debug %DECOMPRESS %COMPRESS);
@EXPORT_OK = qw(&debug $DEBUG);


# this will need to be fixed for non-linux and/or force users to 
# to keep datafiles uncompressed

sub read_marker_refproteins {
    my ($dir,$marker_name) = @_;
    my $marker_file = File::Spec->catfile($dir,$marker_name);
    for my $ext ( qw(fa fas fasta pep aa seq) ) {
	if( -f $marker_file .".$ext" ) {
	    $marker_file .= ".$ext";
	    last;
	}
    }
    my $seqio = Bio::SeqIO->new(-format => 'fasta',
				-file   => $marker_file);
    my $seqs = [];
    while( my $seq =$seqio->next_seq ) {
	push @{$seqs}, $seq;
    }
    $seqs;
}

sub translate_cdna {
    my ($infile,$outfile,$force) = @_;
    my $rc = 1;
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $infile ) {

	my $in = Bio::SeqIO->new(-format => 'fasta', -file => $infile);
	my $out = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile");
	while( my $s = $in->next_seq ) {
	    my $rseq = $s->revcom;
	    my $id = $s->display_id;
	    my $tseq = $s->translate(-frame=> 0,
				     -terminator => 'X');
	    $out->write_seq($tseq);
	}

    }
    $rc;
}
sub make_6frame {
    my ($infile,$outfile,$force) = @_;
    my $rc = 1;
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $infile ) {

	my $in = Bio::SeqIO->new(-format => 'fasta', -file => $infile);
	my $out = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile");
	while( my $s = $in->next_seq ) {
	    my $rseq = $s->revcom;
	    my $id = $s->display_id;
	    for my $frame ( 0,1,2) {
		my $tseq = $s->translate(-frame=> $frame,
					 -terminator => 'X');
		$tseq->display_id(sprintf("%s_%d",
					  $id,$frame+1));
		$out->write_seq($tseq);
		$tseq = $rseq->translate(-frame=> $frame,
					 -terminator => 'X');
		$tseq->display_id(sprintf("%s_%d",
					  $id,$frame+4));
		$out->write_seq($tseq);	    
	    }
	    $rc++;
	}
    } else {
	debug("6frame translation already created");
    }
    $rc;
}

sub fix_read_ids {
    my ($infile,$outfile,$interleaved) = @_;
    my $rc = 1;
    my $fh;
    my $ext;
    if( $infile =~ /\.(gz|bz2)$/ ) {
	$ext = $1;
	open($fh => "$DECOMPRESS{$ext} $infile |") || die "($ext) cannot $DECOMPRESS{$ext} $infile: $!";
    } else {
	open($fh => $infile) || die "Cannot open $infile: $!\n";
    }
    open( my $ofh => ">$outfile") || die("Cannot open $outfile: $!\n");
    while(<$fh>) {
	if( /^>(\S+)(\s+(\S+))?/ ) {	    
	    my ($id,$pair_bc) = ($1,$3);
	    $id =~ s/[-:\/#|]/_/g;
	    if( $interleaved && $pair_bc ) {
		if( $pair_bc =~ /([12]):/ ) {
		    $id .= "_$1";
		} 
	    }
	    $_ = ">$id\n";
	    $rc++;
	} 		
	print $ofh $_;
    } 
    close($fh);
    close($ofh);
    $rc;
}

sub parse_hmmtable {
    my ($infile,$hmmer_cutoff) = @_;
    open(my $fh => $infile) || die "cannot open $infile: $!";
    my $seen;
    while(<$fh>) {
	next if /^\#/;
	my @row = split(/\s+/,$_);
	my $t = $row[0];
	my $q = $row[3];
	# this may be unnecessary
	my $evalue = $row[6];
	next if $evalue > $hmmer_cutoff;
	my $id = $t;
	$id =~ s/_[0-6]$//;
	$seen->{$q}->{$id}++;
    }
    $seen;
}

sub parse_diamond_blastx {
    my $file = shift;
    open(my $fh => $file ) || die "cannot open $file: $!";
    my $seen;
    while(<$fh>) {
	my ($read,$marker, $pid, $mismatch,$gaps, 
	    $gaplen, $qstart,$qend,$tstart,$tend, $evalue,$bits) = split;
	$seen->{$marker}->{$read}++;
    }
    $seen;
}

sub parse_config {
    my $config = shift;
    my $apps = {};
    open(my $fh => $config) || die "cannot open $config: $!";
    while(<$fh>) {
	chomp;
	if(/([^=]+)=(\S+)/ ) {	
	    $apps->{$1} = $2;
	} else {
	    warn("cannot parse line $_\n");
	}
    }
    if( $DEBUG ) {
	while( my ($app,$path) = each %$apps ) {
	    debug("app is $app with path = $path\n");
	}
    }
    $apps;
}


sub seqcount {
    my $file = shift;
    open(my $fh => "grep -c '^>' $file |") || die $!;
    my $n = <$fh>;
    $n =~ s/\s+//g;
    $n;
}

sub seq_lengths {
    my $file = shift;
    my $in = Bio::SeqIO->new(-format => 'fasta', -file => $file);
    my $res = {};
    while( my $s = $in->next_seq ) {
	debug($s->display_id. " ". $s->length,"\n");
	$res->{$s->display_id} = $s->length;
    }
    $res;
}

sub debug {
    my @msg = @_;
    warn(join(" ",@msg)) if $DEBUG;
}

1;
