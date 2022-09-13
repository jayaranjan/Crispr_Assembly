package CrisprAssembly;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);
use Data::Dumper;
use Text::CSV qw ( csv );
use Exporter qw (import);
use REST::Client;
use JSON;
use Data::Dumper;
our @EXPORT_OK = qw ( grch38_slice pass_crispr_id tolist fetch_crispr projection print_projection grch37_slice_feature get_coordinates check_crispr_assembly_input);

sub grch38_slice {
    my $slice = shift;  
    #To get all the gene from the slice
    my @genes = @{ $slice->get_all_Genes() };
     foreach my $gene (@genes) {
        print "\nGene ID: ", $gene->stable_id, "\n";
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

sub tolist {
    my $data = shift;
    my $key = shift;
        if (ref (data->{$key}) eq 'ARRAY') {
            $data->{$key};
  }
        elsif (ref($data->{$key}) eq 'HASH') {
            [$data->{$key}];
  }
        else {
                [];
  }
}

sub pass_crispr_id {
    #my $crispr_id = shift;
    my ($crispr_id_1,$crispr_id_2)=@_;
    my $crispr_counter = 1;
    print "Number of crispr id input ", $crispr_id, "\n";
    foreach my $crispr(@ARGV) {    
        print "Crispr ID # $crispr_counter : $crispr\n";  
        $crispr_counter++;
    }
#return crispr_id_1 , crispr_id_2;
}


sub fetch_crispr {
    my $client = shift;
#   my ($crispr_id_1,$crispr_id_2) = @_;
      # print "\nCrispr Sequence By ID\n",$client->responseContent(),"\n";
    my $response = from_json($client->responseContent());
        print Dumper ($response);
    my $chr_name = tolist ($response -> {'chr_name'});
    my $chr_start =  $response -> {'1106710999'}{'chr_start'};
   #my $chr_start =  $response -> $crispr_id_1,{'chr_start'};
    my $chr_end = $response-> {'1106710999'}{'chr_end'};
    my $pam_right = tolist ($response -> {'pam_right'});
        print "\n",$id,$char_name,"\n";
        print Dumper ($chr_start);
        print Dumper ($chr_end); 
 #return (char_name,char_start,char_end,pam_right);                   
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

# Undef Test ##
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
