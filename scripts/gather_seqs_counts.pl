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
my $cutoff = 1e-10;
my $seqdb;
my $speciesorder;
my $verbose;
my $outseqs;
GetOptions('help|?'         => \$help, man => \$man,
	   'i|input:s'      => \$indir,
           'o|t|table:s'      => \$out,
	   'os|outdir:s'    => \$outseqs,
	   'v|verbose!'     => \$verbose,
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

if( $seqdb && ! -f "$seqdb.cidx" ) {
    `cdbfasta $seqdb`;
}

my (%table, %specieset);
opendir(my $ind => $indir) || die "cannot opne $indir: $!";
for my $file ( readdir($ind) ) {
 next unless $file =~ /(\S+)\Q$ext\E$/;
 my $stem = $1;
 my ($domain_name,$species) = split(/\Q$sep\E/,$stem);
 $species =~ s/\.v\d+//;
 $species =~ s/\.aa//;
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
print $outfh join("\t", qw(DOMAIN), @species), "\n";

# sort by copy number later so abundant one come first?
for my $d ( sort keys %table ) {
 print $outfh join("\t",$d, map { scalar @{$table{$d}->{$_} || []} } @species), "\n";
 if( $seqdb ) {
    
     open(my $cdbfasta => "| cdbyank $seqdb.cidx > $domaindir/$d.fas");
     while ( my ($sp,$domainobs) = each %{$table{$d}} ) 
     { # process domain/species pair data
	 # which is a list of domain positions in the protein
	 my %gene_count;
	 
	 for my $hit ( @{$domainobs || []} ) {
	     my $gene = $hit->[0];
	     print $cdbfasta join("\n", $gene), "\n";
	 }
     }
 }
}
__END__

=head1 NAME

gather_seq_counts.pl - generate a table of seq counts and seqs themselves

=head1 USAGE

perl gather_domaincounts.pl -i results -o table [-ext .domtbl] [-sep __]

=head1 DESCRIPTION

Generate table of domain counts for use in comparing, plotting as heatmap, etc

=head2 ARGUMENTS

-i / --input   Input folder  [One file per domain]
-t / --table  Output table 
-os / --output Output sequences directory
-ext           Extension of the result files. Default is .domtbl
-sep           Separation of the species name from the domain name in the input file
-db            Fasta format protein database

=head1 AUTHOR

Jason Stajich - jasonstajich.phd[at]gmail.com
University of California-Riverside

=cut
