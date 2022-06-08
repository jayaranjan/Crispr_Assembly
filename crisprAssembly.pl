use strict;
use warnings;
use Test::More tests => 35;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);
use Data::Dumper;
use Text::CSV qw ( csv );
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::ExonTranscript;
use Bio::EnsEMBL::TranscriptMapper;
use HTTP::Tiny;
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',    #alternatively 'useastdb.ensembl.org'
    -user => 'anonymous',
    -port => 3337
);

my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
  my $slice = $slice_adaptor->fetch_by_region(qw/chromosome 13   32332312 32332334 1      GRCh37/);
################### To GET ALL THE GENE IN THE SLICE ########################
my @genes = @{ $slice->get_all_Genes() };
foreach my $gene (@genes) {
    printf( "In terms of slice: %d-%d (%+d)\n",
        $gene->start(), $gene->end(), $gene->strand() );
    printf(
        "In terms of seq_region: %d-%d (%+d)\n",
        $gene->seq_region_start(),
        $gene->seq_region_end(),
        $gene->seq_region_strand()
    );
}
