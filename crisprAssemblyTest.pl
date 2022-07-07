
use strict;
use warnings;
use Test::More tests => 9;
ok( crispr_assembly() eq "GRCH38 output required",
    "GRCH38  crispr output required" );
ok(
    crispr_assembly() eq "GRCH37 output required",
    "GRCH37 crispr assembly output required"
);
#test iif the coodrinates are correct  GRCh37 & 38.
#test to get the codinates of GRCh37 using GRCH38 Sequence and check if it is right using the test data ( website ) 
#ok(len_projection('1') eq  "number of projections greater than 1 " , "test passed");
ok( crispr_assembly() eq " slice on old assembly",
    "slice on old assembly required" );
ok( crispr_assembly() eq " slice on new assembly",
    "slice on new assembly required" );

ok( crispr_assembly() eq "not a match", "multiple projection" );
ok( crispr_assembly() eq "match","single projection" );
ok( crispr_assembly() eq "projection into slice", " test projection slice" );
ok( crispr_assembly() eq "updated slice", "updated slice required" );

# test it's still crispr

#is_deeply(
 #   csv( in => 'testoutput.csv', headers => 'auto' ),
  #  csv( in => 'miniSample.csv', headers => 'auto' ),
   # 'compare GRCH37 and GRCH38 (test.csv) ');



