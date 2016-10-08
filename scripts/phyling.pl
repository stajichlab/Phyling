#!env perl
use strict;
use warnings;

use FindBin qw($Bin);
use lib "$FindBin::Bin/../lib/perl";
use File::Copy qw(move);
use Env qw(PHYLINGHOME);
use File::Spec;
use Getopt::Long;
use File::Temp qw(tempfile);
use PHYling qw(debug);
use IO::String;
use Bio::SeqIO;
use Bio::AlignIO;

my $SLEEP_TIME = 240; # sleep 180 sec which should give the blatserver enough time to start up

my $buffer_end_start = 3; # what is the sloppy overhang allowed when bringing in more seqs

my $exonerate_options = '-m p2g --bestn 1 --joinfilter 1 --verbose 0 --ryo ">%ti (%tab - %tae) score=%s rank=%r\n%tcs\n" --showcigar no --showvulgar no --showalignment no --refine full'; # --exhaustive';

my %uncompress = ('bz2' => 'bzcat',
		  'gz'  => 'zcat');

my @EXPECTED_APPS = qw(FASTQ_TO_FASTA HMMALIGN HMMSEARCH TRANSEQ 
                       FASTA TFASTY HMMEMIT CAP3
                       CDBFASTA CDBYANK PHRAP SREFORMAT 
                       TRIMAL FASTTREE MUSCLE EXONERATE DIAMOND);

$ENV{WISECONFIGDIR} = '/opt/linux/centos/7.x/x86_64/pkgs/genewise/2.4.1/';
my $app_conf;

if( $PHYLINGHOME) {
    $app_conf = File::Spec->catfile($PHYLINGHOME, qw(lib apps.conf));

} elsif ($Bin) {
    $app_conf = File::Spec->catfile($Bin, qw(.. lib apps.conf));
}

my $fix_names = 0;
my $qual_offset = 33;
my $hmmer_cutoff = '1e-10';
my $contig_match_cutoff = '1e-8';
my $contig_transsearch_cutoff = '1e-3';
my $Max_rounds = 10; # max iterations
my $scaffold_separator = 'N'x15; # 5 amino acid break in the scaffolded contigs, codon
my ($hmm2_models,$hmm3_models,$marker_hmm,$marker_fasta_dir);
my $clade = 'AFTOL70';
my $CPUs = 1;
my $cleanup = 0;
my $force = 0; # force re-processing files even if there is cached intermediate
my $tmpdir;
my $prefix;
my $seqprefix;
my $rDNA_hmm;
my $do_MSA = 0;
my $consensus_db = 'JGI_1086.consensus';
my $diamond_db;
my $interleaved = 1;
my $consensus_folder = 'consensus';
my $port = 8001+int rand(10000);
my $outdir = 'phyling_run';
GetOptions('ac|app|appconf:s' => \$app_conf,
	   'v|debug!'         => \$PHYling::DEBUG,
	   'force!'           => \$force,
	   'maxrounds:i'      => \$Max_rounds,
	   'q|qual:i'         => \$qual_offset,
	   't|temp|tmpdir:s'  => \$tmpdir,
	   'cleanup!'         => \$cleanup,
	   'p|prefix:s'       => \$prefix,
	   'sp|prefix:s'      => \$seqprefix,	   
	   'interleaved!'     => \$interleaved,
	   'c|clade:s'        => \$clade,
	   'cpus|cpu:i'           => \$CPUs,
	   'hmm:s'            => \$marker_hmm,
	   'hmm2:s'           => \$hmm2_models,
	   'hmm3:s'           => \$hmm3_models,
	   'md|markerdir:s'   => \$marker_fasta_dir,
	   'rDNA!'            => \$rDNA_hmm,
	   'cons|consensus:s' => \$consensus_folder,
	   'consdb:s'         => \$consensus_db,
	   'dmdb:s'           => \$diamond_db,
	   'msa!'             => \$do_MSA,
	   'port:i'           => \$port,
	   'o|out:s'          => \$outdir,
    );

mkdir($consensus_folder) unless -d $consensus_folder;

my $error = 0;
debug("app conf is $app_conf\n");
if( ! $app_conf ) {
    die("Must provide app config file via $PHYLINGHOME or --app --appconf\n");
}
my $paths = &PHYling::parse_config($app_conf);
for my $p ( @EXPECTED_APPS ) {
    if( ! exists $paths->{$p} || ! -x $paths->{$p} ) {
	debug("cannot find Application $p ",$paths->{$p},"\n");
	$error = 1;
    }
}
die("Error in App config file\n") if $error;

$diamond_db = $consensus_db.".dmnd" unless defined $diamond_db;
if( ! -f $diamond_db ) {
    my $cmd = sprintf("%s blastx --threads %d --in %s --db %s",
		      $paths->{'DIAMOND'},$CPUs, $consensus_db, $diamond_db);
    debug("CMD: $cmd\n");
    `$cmd`;
}

my $in_file = shift @ARGV;

if( ! defined $in_file || ! -f $in_file ) {
    die("must provide an input sequencing file to process");
}
	
my (undef,$dir,$fname) = File::Spec->splitpath($in_file);

my $base = $prefix;
if( ! $prefix && $fname =~ /(\S+)\.(fq|fq|fa|fast\w+|seq)/) {
    $base = $1;
} else {
    $base = $$;
}
if( ! $tmpdir ) {
    $tmpdir = File::Spec->catdir($outdir,$base.".PHYling");
} 
if(  ! -d $tmpdir ) {
    mkdir($tmpdir);
}

my $fasta_file;

my ($data_type,$ext);
if( $fname =~ /(\S+)\.(fastq|fq)(\.(gz|bz2))?/ ) {
    $prefix = $1 if ! defined $prefix;
    $ext = $4;
    $fasta_file = File::Spec->catfile($tmpdir,"$1.fasta");
    if( ! -f $fasta_file || $force ) {
	my $cmd;
	if( $fname =~ /\.(bz2|gz)$/ ) {
	    $cmd = sprintf("%s %s | %s -Q %d -o %s",$PHYling::DECOMPRESS{$1},
			   $in_file,
			   $paths->{'FASTQ_TO_FASTA'},$qual_offset,$fasta_file);
	} else {
	    $cmd = sprintf("%s -Q %d -i %s -o %s",$paths->{'FASTQ_TO_FASTA'},
			   $qual_offset, $in_file,$fasta_file);
	}
	debug("CMD: $cmd\n");
	`$cmd`;
    }
    $data_type = 'FASTQ';
} elsif( $fname =~ /(\S+)\.(seq|fasta|fa|fas)(\.(gz|bz2))?/ ) {
    $prefix = $1 if ! defined $prefix;
    $ext = $4;
    $fasta_file = $in_file;
    $data_type = 'FASTA';
} else {
    debug("unknown extension in file $fname ($in_file)\n");
    exit;
}
$seqprefix ||= $prefix;

# GET RID OF PROBLEMATIC CHARACTERS IN THE READ IDs
# some caching, if we already 

my $tmp_rename = File::Spec->catfile($tmpdir,$prefix.".fasta.fix");
if( $force || ! -f $tmp_rename ) {
    &PHYling::fix_read_ids($fasta_file,$tmp_rename,$interleaved);
} else {
    debug("Read file already renamed\n");
}
$fasta_file = $tmp_rename;

#INDEX THE READS FILE FOR LATER PROCESSING
&index_file($fasta_file);

# make 
warn("using port $port\n");
my $bitfile = &make_2bit_file($fasta_file,$force);
my $blat_ready = &start_gfServer($port,$bitfile);

my $marker_table = sprintf("%s.diamond.tab",
			   File::Spec->catfile($tmpdir,$seqprefix));
&run_diamond_markers($fasta_file, $consensus_db, $marker_table);
my $reads_per_marker = &PHYling::parse_diamond_blastx($marker_table);

my @trim_files;
for my $marker ( keys %$reads_per_marker ) {
    warn("marker is $marker\n");
    my $marker_cons = File::Spec->catfile($consensus_folder,"$marker.cons");
    if( ! -f $marker_cons ) {
	die "cannot process without consesnsi already created for now";
	&make_consensus_HMM(File::Spec->catfile($hmm3_models,$marker.".hmm"),
			    $marker_cons);
    } else {
	debug("consensus HMM for $marker_cons already created\n");
    }
    
    #my $marker_seqs = &read_marker_refproteins($marker_fasta_dir,$marker);
    my $cdnafile = File::Spec->catfile($tmpdir,$prefix.".$marker.candidate.cdna");
    my $pepfile = File::Spec->catfile($tmpdir,$prefix.".$marker.candidate.pep");
    my $scaffoldfile = File::Spec->catfile($tmpdir,$prefix.".$marker.ord_scaf.fa");
    debug("cdnafile is $cdnafile\n");
    if( -f $cdnafile ) {
	debug("cdnafile exists\n");
    }
    next if ( ! $force && -f $cdnafile);
    
    if( ! -f $scaffoldfile || $force ) {
	my @reads = keys %{$reads_per_marker->{$marker}};
	my $reads_file = File::Spec->catfile($tmpdir,
					     $prefix.".$marker.r1.fasta");
	if( $force || ! -f $reads_file || 
	    -M $reads_file > -M $marker_table) {
	    &retrieve_reads($fasta_file,\@reads,$reads_file);
	}
	
#	my $contigsfile = &assemble_reads_phrap($reads_file);
	my $contigsfile = &assemble_reads_cap3($reads_file);
	my $contig_count = &PHYling::seqcount($contigsfile);

	if( $contig_count == 0 ) {
	    warn("No assembled contigs or singlets available\n");
	    next;
	}
	debug("seqcount for $contigsfile is $contig_count\n");
	my $change = 1;
	my $rounds = 0;
	while( $contig_count > 1 && 
	       $change > 0 && $rounds < $Max_rounds) {

	    my $added = &search_and_add($fasta_file,$contigsfile,
					$reads_file);
	    warn("added $added reads to $reads_file\n");
	    last if $added == 0;
	    
	    #$contigsfile = &assemble_reads_phrap($reads_file);
	    $contigsfile = &assemble_reads_cap3($reads_file);
	    my $newcount = &PHYling::seqcount($contigsfile);
	    warn("contig count was $contig_count newcount is $newcount\n");
	    $change = ($contig_count - $newcount);
	    $contig_count = $newcount;
	    $rounds++;
	}
	warn("performed $rounds of additional read gathering\n");
	
	my $updated_contigs = &stitch_order_contigs($marker_cons,$contigsfile);
	# merge the contigs, in their new order, into one scaffold with some Ns between
	my $scaff_seq = join($scaffold_separator, (map { defined $_  ? $_->seq : '' } @$updated_contigs));
	my $scaffold = Bio::Seq->new(-id => "$prefix.$marker.scaffold",
				     -seq => $scaff_seq);
	Bio::SeqIO->new(-format => 'fasta', 
			-file =>">$scaffoldfile")->write_seq($scaffold);
	
	debug("Scaffold file: $scaffoldfile\n");
    }
    
    if( -f $scaffoldfile ) {
	&exonerate_best_model($marker_cons,$scaffoldfile,$cdnafile);
	&PHYling::translate_cdna($cdnafile,$pepfile,$force);
    } else {
	warn("no scaffold to process for $marker\n");
    }
    last if $PHYling::DEBUG;
}

if( $rDNA_hmm ) {
    
}

&stop_gfServer($port);

END {
 &stop_gfServer($port);

}

sub start_gfServer {
    my ($Port,$bitfile) = @_;
    #my $pid = fork();
    #die "fork failed" unless defined $pid;
    #if ($pid == 0) {
	# child process goes here
	warn("Port is $Port bitfile is $bitfile\n");
	my $cmd = sprintf("%s start localhost %d %s -canStop",
			  $paths->{GFSERVER},$Port,$bitfile,$bitfile);
	debug("CMD: $cmd\n");
	system("$cmd &");
   # }
    # force sleeping so gf server has enough time to startup
    sleep($SLEEP_TIME);
    return 1;
}
sub stop_gfServer {
    my ($Port) = @_;
	warn("calling stop on $Port");
    system($paths->{GFSERVER},'stop','localhost',$Port);    
}

sub make_2bit_file {
    my ($fasta) = shift;
    my ($base) = $fasta;
    $base .= ".2bit";
    if( $force ||
	! -f $base ||
	-M $base > -M $fasta ) {
	my $cmd = sprintf("%s %s %s",$paths->{FATOTWOBIT},$fasta,$base);
	debug("CMD: $cmd\n");
	`$cmd`;
    }
    $base;
}

sub trim_aln {
    my ($infile) = shift;
    my $outfile = $infile . '.trim';
    my $rc = 1;
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $infile ) {
	my $in = Bio::AlignIO->new(-format => 'clustalw',
				   -file   => $infile);

	my $out = Bio::AlignIO->new(-format => 'fasta',
				    -file   => ">$infile.2");
	if( my $aln = $in->next_aln ) {
	    $aln->map_chars('\.','-');
	    $aln->set_displayname_flat(1);
	    $out->write_aln($aln);	    
	}
	move("$infile.2",$infile);
	my $cmd = sprintf("%s -in %s -out %s -automated1 -fasta",
			  $paths->{TRIMAL},$infile,$outfile);
	debug("CMD: $cmd\n");
	`$cmd`;
    }
    $outfile;
}

sub hmmalign {
    my ($hmm_model, $inseqs, $outfile) = @_;
    my $rc = 1;
    debug("outfile is $outfile\n");
    if( $force || ! -f $outfile ) {
	my ($ifh,$infile) = tempfile('tmpXXXX',UNLINK=>1);
	my $out = Bio::SeqIO->new(-format => 'fasta', -fh => $ifh);
	$out->write_seq(@$inseqs);	
	close($ifh);
	my $cmd = sprintf("%s --trim --amino %s %s > %s",
			  $paths->{HMMALIGN},
			  $hmm_model, $infile,
			  $outfile.".stk");
	debug("CMD: $cmd\n");
	`$cmd`;
	$cmd = sprintf("%s clustal %s > %s",
		       $paths->{SREFORMAT},$outfile.".stk",
		       $outfile);
	debug("CMD: $cmd\n");
	`$cmd`;
    }
    
}

=head2 exonerate_best_model

 Title   : exonerate_best_model
 Usage   :
 Function: Finds best protein2genome alignment with exonerate
 Example :
 Returns : Status of run. Creates a CDS out file
 Args    :


=cut


sub exonerate_best_model {
    my ($inpepfile,$contigfile,$outfile) = @_;
    my $rc = 0;
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $inpepfile ) {

	my $cmd = sprintf("%s %s %s %s > %s",
			  $paths->{EXONERATE},
			  $inpepfile,$contigfile,
			  $exonerate_options,
			  $outfile);
	debug("CMD: $cmd\n");
	$rc = `$cmd`;
    }
    $rc;
}

sub assemble_reads_phrap {
    my $infile = shift;
    my $contigs = $infile.".contigs";
    if( $force || 
	! -f $contigs ||
	-M $infile < -M $contigs ) {	
	my $cmd = sprintf("%s %s",$paths->{PHRAP},$infile);
	debug("CMD: $cmd\n");
	`$cmd`;
	`cat $infile.singlets >> $contigs`;
    }
    if( ! -f "$contigs.renum" ||
	-M $contigs < -M "$contigs.renum" ) {	

	my $renumber = Bio::SeqIO->new(-format => 'fasta',
				       -file   => $contigs);

	my $out = Bio::SeqIO->new(-format => 'fasta',
				  -file   => ">$contigs.renum");
	my $i = 1;
	while( my $s = $renumber->next_seq ) {
	    $s->display_id(sprintf("ctg%d",$i++));
	    $out->write_seq($s);
	}
	$out->close();
	$out = undef;
    }
    "$contigs.renum";
}

sub assemble_reads_cap3 {
    my $infile = shift;
    my $contigs = $infile.".cap.contigs";
    if( $force ||
        ! -f $contigs ||
        -M $infile < -M $contigs ) {
        my $cmd = sprintf("%s %s",$paths->{CAP3},$infile);
        debug("CMD: $cmd\n");
        `$cmd`;
        #`cat $infile.singlets >> $contigs`;
    }
    $contigs;
}

sub retrieve_reads {
    my ($infile,$reads_ar,$outfile) = @_;
    if( ! -f "$infile.cidx" ) {
	&index_file($infile);
    }
    my $cmd = sprintf("| %s %s.cidx > %s",
		      $paths->{CDBYANK},$infile,$outfile);
    debug("CMD: $cmd\n");
    open(my $to_cdbyank => $cmd) || die "Cannot open $cmd: $!\n";
    my $i = 0;
    for my $read ( @$reads_ar ) {
	print $to_cdbyank $read,"\n";
	$i++;
    }
    $i;
}

sub get_read {
    my ($infile,$read_name) = @_;
    if( ! -f "$infile.cidx" ) {
	&index_file($infile);
    }
    my $cmd = sprintf("%s %s.cidx -a %s |",
		      $paths->{CDBYANK},$infile,$read_name);
    debug("CMD: $cmd\n");
    open(my $readseq => $cmd) || die "Cannot open $cmd: $!\n";
    my $seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $readseq);
    my @seqs;
    while(my $s = $seqio->next_seq ) {
	push @seqs, $s;
    }
    @seqs;
}


sub run_diamond_markers {
    my ($seqfile, $seqdb, $output) = @_;
    if( $force || 
	! -f $output ||
	-M $output > -M $seqdb ) {
	my $cmd = sprintf("%s blastx -d %s -q %s -o %s", $paths->{'DIAMOND'}, $seqdb, $seqfile, $output);
	warn("CMD: $cmd\n");
	my $rc = `$cmd`;
    }
    $output;
}

sub run_hmmsearch_markers {
    my ($seqdb, $markerdb,$hmmfilepref) = @_;
    my $table = $hmmfilepref.".hmmsearch.domtbl";
    my $rpt = $hmmfilepref.".hmmsearch.out";
    my $rc = 1;
    if( $force ||
	! -f $table ||
	-M $table > -M $seqdb ) {
	my $cmd = sprintf("%s -E %s --cpu %d --domtblout %s %s %s > %s ",
			  $paths->{HMMSEARCH}, 
			  $hmmer_cutoff,$CPUs,
			  $table,$markerdb,$seqdb,$rpt);
	warn("CMD: $cmd\n");
	$rc = `$cmd`;
    } else {
      sleep($SLEEP_TIME);
    }
    ($table,$rpt);
}




sub make_6frame_transeq {
    my ($infile,$outfile) = @_;
    my $rc = 1;
    if( $force ||
	! -f $outfile ||
	-M $outfile > -M $infile ) {
	my $cmd = sprintf("%s -trim -clean -frame 6 %s %s",
			  $paths->{TRANSEQ},$infile, $outfile);
	debug("CMD: $cmd\n");
	$rc=`$cmd`;
    } else {
	debug("6frame translation EMBOSS already run\n");
    }
    $rc;
}

sub index_file {
    my $infile = shift;
    my $rc = 1;
    if( $force ||
	! -f $infile.".cidx" || 
	-M $infile.".cidx" > -M $infile) {
	my $cmd = sprintf("%s %s",$paths->{CDBFASTA},$infile);
	debug("CMD: $cmd\n");
	$rc = `$cmd`;
    } else {
	debug("index file already created\n");
    }
    $rc;
}



=head2 stitch_order_contigs

 Title   : stitch_order_contigs
 Usage   : &stitch_order_contigs($markerpep,$contigs);
 Function: Reorder and scaffold contigs based on a protein query sequence from 
 Returns : Update contigs file
 Args    :


=cut

sub stitch_order_contigs {
    my ($marker_cons,$contigfile) = @_;
    my $cmd = sprintf("%s -T %d -m 8c -E %s %s %s",
		      $paths->{TFASTY}, $CPUs,
		      $contig_transsearch_cutoff,
		      $marker_cons,
		      $contigfile);
    debug("running $cmd\n");
    open(my $run => "$cmd |") || die "cannot run: $cmd\n";
    my @results;
    while(<$run>) {
	next if /^\#/;
	chomp;
	my ($q,$h,$pid,$match,$mismatch,$gap,$qstart,$qend,
	    $tstart,$tend,$evalue,$bits) = split(/\t/,$_);
	my ($tstrand) = (1,1);	
	if( $tstart > $tend ) { 
	    ($tend,$tstart) = ($tstart,$tend);
	    $tstrand = -1;
	}
	debug("result is $qstart,$qend,$h,$tstart,$tend,$tstrand,$evalue,$bits,$pid \n");
	push @results, [$qstart,$qend,$h,$tstart,$tend,$tstrand,$evalue,$bits,$pid];
    }
    my %contigs;
    my $read_contigs = Bio::SeqIO->new(-format => 'fasta', -file => $contigfile);
    while(my $s = $read_contigs->next_seq ) {
	debug("seq id is ", $s->display_id,"\n");
        $contigs{$s->display_id} = $s;
    }
    my $new_order;
    # sort by query protein alignment order
    my %seen;
    for my $res ( sort { $a->[0] <=> $b->[0] } @results ) {
	debug (join("\t", @$res),"\n");
	next if $seen{$res->[2]}++;
	next if ! defined $contigs{$res->[2]};
	if( $res->[5] < 0 ) {
	    push @$new_order, $contigs{$res->[2]}->revcom;
	} else {
	    push @$new_order, $contigs{$res->[2]};
	}
    }

    $new_order; #bizzare love triangle
}

sub search_and_add {
    my ($searchdb,$queryfile,$outputreads) = @_;

    warn("$queryfile\n");
    my $qlens = &PHYling::seq_lengths($queryfile);
    my $rlens = &PHYling::seq_lengths($outputreads);

    my (undef,$searchdir,$fname) = File::Spec->splitpath($searchdb);
    $searchdir = '.';
# replace this with BLAT and near-identity?
#    my $cmd = sprintf("%s -T %d -E %s -m 8c %s %s |",
#		      $paths->{FASTA},$CPUs,$contig_match_cutoff,
#		      $queryfile, $searchdb);
    my $cmd = sprintf("%s %s %s %s %s stdout -out=blast8 |",
		      $paths->{GFCLIENT},'localhost',$port,$searchdir,
		      $queryfile);
    
    debug("CMD: $cmd\n");
    open(my $fasta_res => $cmd) || die "cannot run $cmd\n";
    my @results;
    my %readnames;
    while(<$fasta_res>) {
	next if /^\#/;
	debug($_);
	chomp;
	my ($q,$h,$pid,$match,$mismatch,$gap,$qstart,$qend,
	    $tstart,$tend,$evalue) = split(/\t/,$_);
	next if( exists $rlens->{$h} );
	my ($qstrand,$tstrand) = (1,1);
	if( $qstart > $qend ) { 
	    ($qend,$qstart) = ($qstart,$qend);
	    $qstrand = -1;
	}
	if( $tstart > $tend ) { 
	    ($tend,$tstart) = ($tstart,$tend);
	    $tstrand = -1;
	}

	if($tstart > $buffer_end_start ) { # if target alignment start is not 1 or some 
	                                   # number close to 1 (buffer_end_start)
	    $readnames{$h}++;
	} elsif( abs($qlens->{$q}-$qend) >= $buffer_end_start ) {
	    # -----|        Query
	    # ------------| Hit
	    debug("$q overhanging $qend vs length: ".$qlens->{$q}."\n");
	    my ($read_seq) = &get_read($searchdb,$h);
	    my $read_len = $read_seq->length;
	    if( $read_len > $tend) { # if the end of read align (tend) after the end of this aln
		$readnames{$h}++;
	    }
	}
    }
    debug("readnames are ",sort keys %readnames, "\n");
    if( keys %readnames ) {
	&retrieve_reads($searchdb,[sort keys %readnames],"$outputreads.add");
	my $in = Bio::SeqIO->new(-format => 'fasta',
				 -file   => "$outputreads.add");
	my $out = Bio::SeqIO->new(-format => 'fasta',
				  -file   => ">>$outputreads");
	while( my $s = $in->next_seq ) {
	    $out->write_seq($s);
	}
    }
    return scalar keys %readnames;
}

sub make_consensus_HMM {
    my ($hmmfile,$outfile) = @_;
    my $cmd = sprintf("%s -c %s > $outfile",
		      $paths->{HMMEMIT},$hmmfile);
    `$cmd`;
}

END {
    if( $cleanup ) {	
	warn("rm -rf $tmpdir\n");
    }
}
