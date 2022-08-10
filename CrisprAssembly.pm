package CrisprAssembly;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);
use Data::Dumper;
use Text::CSV qw ( csv );
use Exporter qw (import);

our @EXPORT_OK = qw ( grch38_slice projection print_projection grch37_slice_feature get_coordinates check_crispr_assembly_input);

sub grch38_slice {
    my $slice = shift;  
                                            #To get all the gene from the slice#
    my @genes = @{ $slice->get_all_Genes() };
     foreach my $gene (@genes) {
        print "Gene ID: ", $gene->stable_id, "\n";
        printf( " In terms of slice: %d-%d (%+d)\n",
            $gene->start(), $gene->end(), $gene->strand() 
        );
        printf(
            "In terms of seq_region: %d-%d (%+d)\n",
            $gene->seq_region_start(),
            $gene->seq_region_end(),
            $gene->seq_region_strand()
        );
    }
  return @genes;
}

sub projection {
    my $slice = shift;
    my @projection = @{ $slice->project(qw/chromosome GRCh37/) };
    my $len_projection = scalar @projection;
       if ( $len_projection > 1 ) {
        die "number of projections greater than 1 ";
        }
    my $updated = $projection[0]->to_Slice();
 #return @projection,$updated;
}

sub print_projection {
    my $assembly = shift;
    my $updated = shift;  
        printf(
        "\n " . $assembly . " Coordinates: %s:%d-%d \n%s \n",
        $updated->seq_region_name, $updated->start,
        $updated->end,             $updated->seq);

  return $updated;   
}

sub grch37_slice_feature {
    my $slice = shift;
    my @feature = @{ $slice->get_all_Genes() };
        foreach my $feature (@feature) {
        printf(
            "Projected slice Feature at GRCh37: %s %d-%d (%+d)\n",
            $feature->seq_region_name(), $feature->start(),
            $feature->end(),             $feature->strand());
} 
    my $projection = $feature->project('clone');
        foreach my $segment ( @{$projection} ) {
    my $to_slice = $segment->to_Slice();
        printf(
            "  %s %d-%d (%+d)\n",
            $to_slice->seq_region_name(), $to_slice->start(),
            $to_slice->end(),             $to_slice->strand());
    }  
 return $slice,@feature,$to_slice;
}

sub get_coordinates {
    my ($slice,$updated) = @_;
   # $sequence = $slice->subseq( 100, 200 );
    my $slice_grch38 = $slice; 
    my $coord_sys  = $slice_grch38->coord_system()->name();
    my $seq_region = $slice_grch38->seq_region_name();
    my $start      = $slice_grch38->start();
    my $end        = $slice_grch38->end();
    my $strand     = $slice_grch38->strand();

         ## Undef Test ##
     if ( defined $coord_sys) {
            print "Test 1 Slice PASSED: $coord_sys $seq_region $start-$end ($strand)\n";
      }
     else {
          print " Test 1 Slice FAILED: The coordinates does not exist in GRCh37\n";
      }
   return $coord_sys, $seq_region, $start, $end, $strand;   
}

sub check_crispr_assembly_input {
    my ( $old_crispr, $new_crispr ) = @_;
        if ( !$old_crispr ) {
            return "GRCH38 output required";
        }
        if ( !$new_crispr ) {
            return "GRCH37 output required";
        }
}
1;

__END__
