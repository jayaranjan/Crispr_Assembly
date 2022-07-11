use FindBin;
use lib "$FindBin::Bin";
use CrisprAssembly qw( check_crispr_assembly_input Get_Coordinates);
use strict;
use warnings;
use Test::More tests => 9;
 #failing tests.
 #test to get crispr site from the old assembly
ok( check_crispr_assembly_input() eq "GRCH38 output required",
    "GRCH38  crispr output required" );
#ok( check_crispr_assembly_input() eq "GRCH38 output required",
 #   "GRCH37 crispr assembly output required" );
#ok( check_crispr_assembly_input() eq " slice on old assembly",
#    "slice on old assembly required" );
#ok( check_crispr_assembly_input() eq " slice on new assembly",
 #   "slice on new assembly required" );
#ok( check_crispr_assembly_input() eq "not a match", "multiple projection" );
#ok( check_crispr_assembly_input() eq "match",       "single projection" );
#ok( check_crispr_assembly_input() eq "projection into slice",
#    " test projection slice" );
#ok( check_crispr_assembly_input() eq "updated slice",
 #   "updated slice required" );
# test it's still crispr
# if (defined $coord_sys and $seq_region)
# {
# print "Slice: $coord_sys $seq_region $start-$end ($strand)\n";
# ok (1,'the slice is defined');
# }
# else {
# print " The coordinates does not exist in GRCh37\n";
# ok(0);
# }
# ok (defined $coord_sys, 'The coordinates are defined');


ok (defined $updated->seq, 'projection of new sequence');
 ok (defined Get_Coordinates( $coord_sys), 'The coordinates are defined');
  ok (defined Get_Coordinates($seq_region), 'Test sequence region input');
 ok (defined Get_Coordinates($start), 'Test slice start input');
 ok (defined Get_Coordinates($end), 'Test slice end input');
ok (defined Get_Coordinates($strand), 'Test slice strand input');
    
