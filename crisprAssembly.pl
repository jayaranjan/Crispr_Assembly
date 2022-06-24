use strict;
use warnings;
use Test::More tests => 9;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);
use Data::Dumper;
use Text::CSV qw ( csv );

#use Data::Compare;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',    #alternatively 'useastdb.ensembl.org'
    -user => 'anonymous',
    -port => 3337
);

my $slice_adaptor = $registry->get_adaptor(qw/Human core slice/);
my $slice =
  $slice_adaptor->fetch_by_region(qw/chromosome 4 99304971 99352760 /);

################### To GET ALL THE GENE IN THE SLICE ########################
my @genes = @{ $slice->get_all_Genes() };
foreach my $gene (@genes) {
    printf( " In terms of slice: %d-%d (%+d)\n",
        $gene->start(), $gene->end(), $gene->strand() );
    printf(
        "In terms of seq_region: %d-%d (%+d)\n",
        $gene->seq_region_start(),
        $gene->seq_region_end(),
        $gene->seq_region_strand()
    );

}
my @projection = @{ $slice->project(qw/chromosome GRCh37/) };

my $len_projection = scalar @projection;
if ( $len_projection > 1 ) {
    die "number of projections greater than 1 ";
}
my $updated = $projection[0]->to_Slice();
printf(
    "\n%s:%d-%d %s \n",
    $updated->seq_region_name, $updated->start,
    $updated->end,             $updated->seq
);

#print $projection[0]->seq;
# prints "CTATATGAGTTGAGTTCATCTGG",
#my $slice_38 = $slice_adaptor->fetch_by_region(
#    qw/chromosome 13 31758155 317581771 1  GRCh38/);

#failing tests.
#test to get crispr site from the old assembly.

sub crispr_assembly {
    my ( $old_crispr, $new_crispr ) = @_;
    if ( !$old_crispr ) {
        return "GRCH38 output required";
    }
    if ( !$new_crispr ) {
        return "GRCH37 output required";
    }

}
ok( crispr_assembly() eq "GRCH38 output required",
    "GRCH38  crispr output required" );
ok(
    crispr_assembly() eq "GRCH37 output required",
    "GRCH37 crispr assembly output required"
);

#ok(len_projection('1') eq  "number of projections greater than 1 " , "test passed");
ok( crispr_assembly() eq " slice on old assembly",
    "slice on old assembly required" );
ok( crispr_assembly() eq " slice on new assembly",
    "slice on new assembly required" );

ok( crispr_assembly() eq "not a match",           "multiple projection" );
ok( crispr_assembly() eq "match",                 "single projection" );
ok( crispr_assembly() eq "projection into slice", " test projection slice" );
ok( crispr_assembly() eq "updated slice",         "updated slice required" );

# test it's still crispr
#my $adaptor = $registry->get_adaptor(qw/Human core slice/);
#my $slice =
#  $adaptor->fetch_by_region(qw/chromosome 13 32332292 32332314 1 GRCh37/);

#my @projection = @{$slice->project(qw/chromosome  GRCh38/)};
#print $projection[0]->seq;
# prints "CTATATGAGTTGAGTTCATCTGG",

is_deeply(
    csv( in => 'testoutput.csv', headers => 'auto' ),
    csv( in => 'miniSample.csv', headers => 'auto' ),
    'compare GRCH37 and GRCH38 (test.csv) '
);

