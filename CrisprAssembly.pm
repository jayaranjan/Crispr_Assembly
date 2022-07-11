package CrisprAssembly;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);
use Data::Dumper;
use Text::CSV qw ( csv );
use Exporter qw (import);
use Test::More tests => 5;
our @EXPORT_OK = qw ( check_crispr_assembly_input  Get_Coordinates);

#sub Get_Coordinates {

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',    #alternatively 'useastdb.ensembl.org'
    -user => 'anonymous',
    -port => 3337
);

my $slice_adaptor = $registry->get_adaptor(qw/Human core slice/);
my $slice =
  $slice_adaptor->fetch_by_region(
    qw/chromosome 13 32332370 32332992 1  GRCh38/);

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
    "\n GRCh37 Coordinates: %s:%d-%d \nGRCH37\n%s \n",
    $updated->seq_region_name, $updated->start,
    $updated->end,             $updated->seq
);

print " \nGRCH38\n", $updated->seq, "\n";
#################
my @feature = @{ $slice->get_all_Genes() };
foreach my $feature (@feature) {
    printf(
        "Feature at: %s %d-%d (%+d) projects to\n",
        $feature->seq_region_name(), $feature->start(),
        $feature->end(),             $feature->strand()
    );

    my $projection = $feature->project('clone');

    foreach my $segment ( @{$projection} ) {
        my $to_slice = $segment->to_Slice();

        printf(
            "\t%s %d-%d (%+d)\n",
            $to_slice->seq_region_name(), $to_slice->start(),
            $to_slice->end(),             $to_slice->strand()
        );
    }

}

my $coordinate_adaptor =
  $registry->get_adaptor( 'Human', 'Core', 'CoordSystem' );
my $coordinate = $coordinate_adaptor->fetch_by_name('chromosome');
printf "Coordinate system: %s %s %s %s \n", $coordinate->name(),
  $coordinate->version(), $coordinate->dbID(),
  $coordinate->is_sequence_level();

#sub Get_Coordinates {

    #print $sequence, "Test seq\n";
    $sequence = $slice->subseq( 100, 200 );
    my $coord_sys  = $slice->coord_system()->name();
    my $seq_region = $slice->seq_region_name();
    my $start      = $slice->start();
    my $end        = $slice->end();
    my $strand     = $slice->strand();
    print "GRCH38 Coordinates: $coord_sys $seq_region $start-$end ($strand)\n";

    # return $coord_sys, $seq_region, $start, $end, $strand;
#}

# if ( defined $coord_sys) {
#        print "Test 1 Slice: $coord_sys $seq_region $start-$end ($strand)\n";
#ok (1,'the slice is defined');
#   }
#   else {
#          print " ($coord_sys)Test t1 The coordinates does not exist in GRCh37\n";
#ok(0);
# }
if ($updated->seq =~ m[CTAACCCTTTCAGGTCTAAATGG])  
{ 
    print "Test Pass : CRISPR Sequence Match \n"; 
} 
else 
{ 
    print "Test Failed : CRISPR Sequence Match \n"; 
}
sub check_crispr_assembly_input {
    my ( $old_crispr, $new_crispr ) = @_;

    #my $new_crispr = $updated->seq;
    if ( !$old_crispr ) {
        return "GRCH38 output required";
    }
    if ( !$new_crispr ) {
        return "GRCH37 output required";
    }
}

1;


#ok( defined Get_Coordinates($coord_sys), 'The coordinates are defined' );
ok( defined $slice_adaptor, 'To test if GRCh37 slice is defined' );
ok( defined $slice, 'To test if GRCH38 coordinates is defined' );
ok( defined $updated->seq, 'The GRCH37 sequence projection' );
ok( defined $coord_sys, 'The coordinates are defined' );
ok (defined $seq_region, 'Test sequence region input');
ok (defined $start, 'Test slice start input');
ok (defined $end, 'Test slice end input');
ok (defined $strand, 'Test slice strand input');
ok($updated->seq =~ m[CTAACCCTTTCAGGTCTAAATGG],'crispr test');
