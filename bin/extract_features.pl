#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Getopt::Long;

# Help message
sub usage {
    print <<"END_USAGE";
Usage: $0 --reference <input_file> --refdir <output_directory>

Required arguments:
  --reference   Input GenBank file (e.g., NC_002516v2.gbk)
  --refdir     Output directory for ref.fa and ref.gff

Optional arguments:
  --help       Show this help message
END_USAGE
    exit;
}

# Get command line arguments
my ($reference, $refdir, $help);
GetOptions(
    'reference=s' => \$reference,
    'refdir=s'    => \$refdir,
    'help'        => \$help
) or usage();

usage() if $help;
die "Missing required --reference argument\n" unless $reference;
die "Missing required --refdir argument\n" unless $refdir;

# Configuration
my $EXE = "feature_extractor";         # Your tool name
my $ref_fmt = "genbank";              # Input format (can be made configurable)

# Error handling subroutine
sub err { die "[ERROR] @_\n"; }

# Create output directory
mkdir $refdir or err "Cannot create output directory '$refdir': $!" unless -d $refdir;

# Main processing
print "Extracting FASTA and GFF from $reference to $refdir/\n";
my $in = Bio::SeqIO->new(-file=>$reference, -format=>$ref_fmt) 
    or err "Could not open --reference: $reference ($!)";

my $out = Bio::SeqIO->new(-file=>">$refdir/ref.fa", -format=>'fasta');
my $gff = Bio::Tools::GFF->new(-file=>">$refdir/ref.gff", -gff_version=>3);

my $nseq = 0;
my $nfeat = 0;
my %refseq;
my %tagcnt;

while (my $seq = $in->next_seq) {
    exists $refseq{$seq->id} and err("Duplicate sequence ".$seq->id." in $reference");
    
    # Process sequence
    my $dna = uc($seq->seq);
    $dna =~ s/[^AGTCN]/n/g;
    $refseq{ $seq->id } = $dna;
    $seq->seq($dna);
    $out->write_seq($seq);
    $nseq++;
    
    # Process features
    for my $f ($seq->get_SeqFeatures) {
        my $ftype = $f->primary_tag;
        next if $ftype =~ m/^(source|gene|misc_feature)$/;
        $tagcnt{ $ftype }++;
        
        $f->frame( $ftype eq 'CDS' ? '0' : '.' );
        
        if ($f->has_tag('locus_tag')) {
            my($id) = $f->get_tag_values('locus_tag');
            $f->add_tag_value('ID', $id);
        }
        else {
            $f->add_tag_value( 'ID', $ftype.'_'.$tagcnt{$ftype} );
        }
        
        if ($f->has_tag('gene')) {
            my($gene) = $f->get_tag_values('gene');
            $f->add_tag_value('Name', $gene);
        }
        
        $f->source_tag($EXE);
        $gff->write_feature($f);
        $nfeat++;
    }
}

print "Successfully processed $nseq sequences with $nfeat features\n";
print "Output files created:\n";
print "  - FASTA: $refdir/ref.fa\n";
print "  - GFF3: $refdir/ref.gff\n";