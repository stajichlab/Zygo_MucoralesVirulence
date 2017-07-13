#!/usr/bin/perl
#
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Bio::DB::Fasta;
use Bio::SeqIO;

my $man = 0;
my $help = 0;
my ($indir,$domaindir) = qw(results domainseq);
my $out; # can be empty will print to stdout 
my $ext = '.domtbl';
my $sep = '__';
my $cutoff = 1e-4;
my $seqdb;
my $speciesorder;
my $verbose;

GetOptions('help|?'    => \$help, man => \$man,
	   'i|input:s' => \$indir,
           'o|out:s'     => \$out,
	   'v|verbose!'  => \$verbose,
	   'd|domain|domainout:s' => \$domaindir,
	   'c|cutoff|evalue:s' => \$cutoff,
	   'ext|extension:s' => \$ext,
	   'sep|separator:s' => \$sep,
	   'db|database|seqs:s' => \$seqdb,
	   'species:s' => \$speciesorder,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;

mkdir($domaindir) unless -d $domaindir;

if ( $ext !~ /^\./) {
  $ext = '.'. $ext;
}

my (%table, %specieset);
opendir(my $ind => $indir) || die "cannot opne $indir: $!";
for my $file ( readdir($ind) ) {
 next unless $file =~ /(\S+)\Q$ext\E$/;
 my $stem = $1;
 my ($domain_name,$species) = split(/\Q$sep\E/,$stem);
 warn("domain=$domain_name species=$species\n") if $verbose;
 my $filepath = File::Spec->catfile($indir,$file);
 open(my $in => $filepath ) || die "cannot open $filepath: $!";
 while(<$in>) {  # parse domtbl file
  next if /^\#/; # skip comment lines
  my ($gene_name,$qacc,$qlen,$domain,$domacc,$tlen,
      $fullevalue,$fullscore,$fullbias,$n,$ntotal,$cvalue,$ivalue,
      $score,$dombias,
      $hstart,$hend, $qstart,$qend,$envfrom,$envto,$acc,$desc) =
	     split(/\s+/,$_,23);
	my $evalue = $ivalue; 
      if( $evalue > $cutoff ) {
	# skip
	next
      }
     push @{$table{$domain}->{$species}}, [$gene_name,$hstart,$hend];
     $specieset{$species}++;
 }
}
my $db;
if ( $seqdb ) { 
 $db = Bio::DB::Fasta->new($seqdb);
}

my @species;
if ( $speciesorder ) {
 open(my $tmp => $speciesorder ) || die $!;
 while(<$tmp>) {
  my @lst = split;
  push @species, @lst;
 }
} else {
 @species = sort keys %specieset;
}
my $outfh;
if( $out ) {
 open($outfh => ">$out") || die "cannot open $out: $!";
} else {
 $outfh = \*STDOUT;
}
print $outfh join(",", qw(DOMAIN), @species), "\n";

# sort by copy number later so abundant one come first?
for my $d ( sort keys %table ) {
 print $outfh join(",",$d, map { scalar @{$table{$d}->{$_} || []} } @species), "\n";
 if( $db ) {
   my $outseq = Bio::SeqIO->new(-format => 'fasta', -file => ">$domaindir/$d.fas");
   while ( my ($sp,$domainobs) = each %{$table{$d}} ) { # process domain/species pair data
							 # which is a list of domain positions in the protein
     my %gene_count;
     for my $domain_pos ( @{$domainobs || []} ) { # there are observations of this domain
	my ($proteinid,$start,$end) = @$domain_pos;
	my $pepseq = $db->seq($proteinid,$start,$end);
         
	my $id = sprintf("%s_DOM%02d",$proteinid,$gene_count{$proteinid}++);
	$outseq->write_seq(Bio::PrimarySeq->new(-seq => $pepseq, -id => $id, -desc => sprintf("%s:%d..%d",$proteinid,$start,$end)));
     }
   }
 }
}
__END__

=head1 NAME

gather_domaincounts.pl - generate a table of domain counts

=head1 USAGE

perl gather_domaincounts.pl -i results -o table [-ext .domtbl] [-sep __]

=head1 DESCRIPTION

Generate table of domain counts for use in comparing, plotting as heatmap, etc

=head2 ARGUMENTS

-i / --input   Input folder  [One file per domain]
-o / --output  Output table 
-ext           Extension of the result files. Default is .domtbl
-sep           Separation of the species name from the domain name in the input file

=head1 AUTHOR

Jason Stajich - jasonstajich.phd[at]gmail.com
University of California-Riverside

=cut
