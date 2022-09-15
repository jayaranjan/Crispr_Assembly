use FindBin;
use lib "$FindBin::Bin";
use CrisprAssembly qw ( grch38_slice pass_crispr_id fetch_crispr projection print_projection grch37_slice_feature get_coordinates check_crispr_assembly_input);
use strict;
use warnings;
use Test::More tests => 18;
use JSON;
use Data::Dumper; 
my $registry = 'Bio::EnsEMBL::Registry';
   $registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');
my $slice_adaptor = $registry->get_adaptor(qw/human core slice/);
# the slice have to get coordinates from multiple crisper ids and out put the infromation about the seq (use for loop )or (if else).
my $slice = $slice_adaptor->fetch_by_region(qw/chromosome 13 32332370 32332992 1 GRCh38 /);
my ($crispr_id_1, $crispr_id_2) = @ARGV;
    &pass_crispr_id ($crispr_id_1,$crispr_id_2);
#&fetch_crispr($crispr_id_1);
my $client = REST::Client->new();
   $client->GET('https://wge.stemcell.sanger.ac.uk/api/crispr_by_id?species=Grch38&id='.$crispr_id_1.'&id='.$crispr_id_2.'');
    &fetch_crispr($client, $crispr_id_1, $crispr_id_2);
   # &pass_crispr_id($client);
my @projection = @{ $slice->project(qw/chromosome GRCh37/) };
my @grch38_genes = &grch38_slice($slice);
my $gene = $grch38_genes[0];
my $gene_start = $gene->seq_region_start();
my $gene_end = $gene->seq_region_end();
my $gene_strand = $gene->seq_region_strand();
my $gene_id =  $gene->stable_id;
    &print_projection("GRCh38", $slice);
my $updated = &projection($slice);
    &print_projection("GRCh37", $updated);
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
is( $gene_start , "32315086", 'Test gene slice start '); 
is( $gene_end , "32400268", 'Test gene slice end '); 
is( $gene_strand , "1", 'Test gene slice strand'); 
is( $coordinates_name , "13", 'Test slice coordinates name'); 
is( $coordinates_start , "32906507", 'Test slice coordinates start'); 
is( $coordinates_end , "32907129", 'Test slice coordinates end'); 
is( $feature_name , "13", 'Test slice feature projection name'); 
is( $feature_start , "-17283", 'Test slice feature projection start'); 
is( $feature_end , "67899", 'Test slice feature projection end'); 
is( $feature_strand , "1", 'Test slice feature projection strand'); 
is( $check_grch38_name , "chromosome", 'Test GRCh38 coordinates name'); 
is( $check_grch38_region , "13", 'Test GRCh38 coordinates region'); 
is( $check_grch38_start , "32332370", 'Test GRCh38 coordinates start'); 
is( $check_grch38_end , "32332992", 'Test GRCh38 coordinates end'); 
is( $check_grch38_strand , "1", 'Test GRCh38 coordinates strand'); 
is( check_crispr_assembly_input(), "GRCH38 output required", "DNA sequence required", );
is( $coordinates_seq ,"ACAGTTGTAGATACCTCTGAAGAAGATAGTTTTTCATTATGTTTTTCTAAATGTAGAACAAAAAATCTACAAAAAGTAAGAACTAGCAAGACTAGGAAAAAAATTTTCCATGAAGCAAACGCTGATGAATGTGAAAAATCTAAAAACCAAGTGAAAGAAAAATACTCATTTGTATCTGAAGTGGAACCAAATGATACTGATCCATTAGATTCAAATGTAGCAAATCAGAAGCCCTTTGAGAGTGGAAGTGACAAAATCTCCAAGGAAGTTGTACCGTCTTTGGCCTGTGAATGGTCTCAACTAACCCTTTCAGGTCTAAATGGAGCCCAGATGGAGAAAATACCCCTATTGCATATTTCTTCATGTGACCAAAATATTTCAGAAAAAGACCTATTAGACACAGAGAACAAAAGAAAGAAAGATTTTCTTACTTCAGAGAATTCTTTGCCACGTATTTCTAGCCTACCAAAATCAGAGAAGCCATTAAATGAGGAAACAGTGGTAAATAAGAGAGATGAAGAGCAGCATCTTGAATCTCATACAGACTGCATTCTTGCAGTAAAGCAGGCAATATCTGGAACTTCTCCAGTGGCTTCTTCATTTCAGGGTATCAAAAAGTCTAT", 'test sequence of GRCh37');
