use FindBin;
use lib "$FindBin::Bin";
use CrisprAssembly qw ( grch38_slice projection update_projection get_coordinates check_crispr_assembly_input);
use strict;
use warnings;
use Test::More tests => 18;
 
my $registry = 'Bio::EnsEMBL::Registry';
 $registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous', -port => 3337);
 my $slice_adaptor = $registry->get_adaptor(qw/Human core slice/);
my $slice = $slice_adaptor->fetch_by_region(qw/chromosome 13 32332370 32332992 1  GRCh38/);
my @projection = @{ $slice->project(qw/chromosome GRCh37/) };
 
my @grch38_genes = &grch38_slice($slice);
my $gene = $grch38_genes[0];
my $gene_start = $gene->seq_region_start();
my $gene_end = $gene->seq_region_end();
my $gene_strand = $gene->seq_region_strand();
my $gene_id =  $gene->stable_id;
my $updated = &projection($slice);

my $coordinates_name = $updated->seq_region_name;
my $coordinates_start = $updated->start;
my $coordinates_end = $updated->end;
my $coordinates_seq = $updated->seq;
my $feature = $gene;
my $feature_name = $feature->seq_region_name();
my $feature_start = $feature->start();
my $feature_end = $feature->end();
my $feature_strand = $feature->strand();
my $slice_grch38 = &get_coordinates($slice);
my $check_grch38_name = $slice->coord_system()->name();
my $check_grch38_region = $slice->seq_region_name();
my $check_grch38_start = $slice->start();
my $check_grch38_end = $slice->end();
my $check_grch38_strand = $slice->strand();

 #failing tests.
is ($gene_id, "ENSG00000139618",'Test Gene ID');
is( $gene_start , "32315474", 'Test gene slice start '); 
is( $gene_end , "32399668", 'Test gene slice end '); 
is( $gene_strand , "1", 'Test gene slice strand'); 
is( $coordinates_name , "13", 'Test slice coordinates name'); 
is( $coordinates_start , "32906507", 'Test slice coordinates start'); 
is( $coordinates_end , "32907129", 'Test slice coordinates end'); 
is( $feature_name , "13", 'Test slice feature projection name'); 
is( $feature_start , "-16895", 'Test slice feature projection start'); 
is( $feature_end , "67299", 'Test slice feature projection end'); 
is( $feature_strand , "1", 'Test slice feature projection strand'); 
is( $check_grch38_name , "chromosome", 'Test GRCh38 coordinates name'); 
is( $check_grch38_region , "13", 'Test GRCh38 coordinates region'); 
is( $check_grch38_start , "32332370", 'Test GRCh38 coordinates start'); 
is( $check_grch38_end , "32332992", 'Test GRCh38 coordinates end'); 
is( $check_grch38_strand , "1", 'Test GRCh38 coordinates strand'); 
is( check_crispr_assembly_input(), "GRCH38 output required", "DNA sequence required", );
is( $coordinates_seq ,"ACAGTTGTAGATACCTCTGAAGAAGATAGTTTTTCATTATGTTTTTCTAAATGTAGAACAAAAAATCTACAAAAAGTAAGAACTAGCAAGACTAGGAAAAAAATTTTCCATGAAGCAAACGCTGATGAATGTGAAAAATCTAAAAACCAAGTGAAAGAAAAATACTCATTTGTATCTGAAGTGGAACCAAATGATACTGATCCATTAGATTCAAATGTAGCAAATCAGAAGCCCTTTGAGAGTGGAAGTGACAAAATCTCCAAGGAAGTTGTACCGTCTTTGGCCTGTGAATGGTCTCAACTAACCCTTTCAGGTCTAAATGGAGCCCAGATGGAGAAAATACCCCTATTGCATATTTCTTCATGTGACCAAAATATTTCAGAAAAAGACCTATTAGACACAGAGAACAAAAGAAAGAAAGATTTTCTTACTTCAGAGAATTCTTTGCCACGTATTTCTAGCCTACCAAAATCAGAGAAGCCATTAAATGAGGAAACAGTGGTAAATAAGAGAGATGAAGAGCAGCATCTTGAATCTCATACAGACTGCATTCTTGCAGTAAAGCAGGCAATATCTGGAACTTCTCCAGTGGCTTCTTCATTTCAGGGTATCAAAAAGTCTAT", 'test sequence of GRCh37');
