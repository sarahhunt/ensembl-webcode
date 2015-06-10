=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package EnsEMBL::Web::Role::Gene::Variation;

use Role::Tiny;

sub get_fg_db {
  my $self = shift;
  my $slice = $self->get_Slice( @_ );
  my $fg_db = undef;
  my $db_type  = 'funcgen';

  unless($slice->isa("Bio::EnsEMBL::Compara::AlignSlice::Slice")) {
    $fg_db = $slice->adaptor->db->get_db_adaptor($db_type);
    if(!$fg_db) {
      warn("Cannot connect to $db_type db");
      return [];
    }
  }

  return $fg_db;
}

sub get_feature_view_link {
  my ($self, $feature) = @_;
  my $feature_id  = $feature->display_label;
  my $feature_set = $feature->feature_set->name;

  return if $feature_set =~ /cisRED|CRM/i;

  my $link = $self->hub->url({
    type   => 'Location',
    action => 'Genome',
    ftype  => 'RegulatoryFactor',
    fset   =>  $feature_set,
    id     =>  $feature_id,
  });

  return qq{<span class="small"><a href="$link">[view all]</a></span>};
}

sub get_extended_reg_region_slice {
  my $self = shift;
  ## retrieve default slice
  my $object_slice = $self->Obj->feature_Slice;
     $object_slice = $object_slice->invert if $object_slice->strand < 1; ## Put back onto correct strand!


  my $fg_db = $self->get_fg_db;
  my $fg_slice_adaptor = $fg_db->get_SliceAdaptor;
  my $fsets = $self->feature_sets;
  my $gr_slice = $fg_slice_adaptor->fetch_by_Gene_FeatureSets($self->Obj, $fsets);
  $gr_slice = $gr_slice->invert if $gr_slice->strand < 1; ## Put back onto correct strand!


  ## Now we need to extend the slice!! Default is to add 2kb to either end of slice, if gene_reg slice is
  ## extends more than this use the values returned from this
  my $start = $self->Obj->start;
  my $end   = $self->Obj->end;

  my $gr_start = $gr_slice->start;
  my $gr_end = $gr_slice->end;
  my ($new_start, $new_end);

  if ( ($start  - 2000) < $gr_start) {
     $new_start = 2000;
  } else {
     $new_start = $start - $gr_start;
  }

  if ( ($end +2000) > $gr_end) {
    $new_end = 2000;
  }else {
    $new_end = $gr_end - $end;
  }

  my $extended_slice =  $object_slice->expand($new_start, $new_end);
  return $extended_slice;
}

sub feature_sets {
  my $self = shift;

  my $available_sets = $self->species_defs->databases->{'DATABASE_FUNCGEN'}->{'FEATURE_SETS'};
  my $fg_db = $self->get_fg_db;
  my $feature_set_adaptor = $fg_db->get_FeatureSetAdaptor;
  my @fsets;

  foreach my $name ( @$available_sets){
    push @fsets, $feature_set_adaptor->fetch_by_name($name);
  }
  return \@fsets;
}

sub reg_factors {
  my $self = shift;
  my $gene = $self->gene;
  my $fsets = $self->feature_sets;
  my $fg_db= $self->get_fg_db;
  my $ext_feat_adaptor = $fg_db->get_ExternalFeatureAdaptor;
  my $fg_slice_adaptor = $fg_db->get_SliceAdaptor;
  my $slice = $self->get_extended_reg_region_slice;
  my $factors_by_gene = $ext_feat_adaptor->fetch_all_by_Gene_FeatureSets($gene, $fsets, 1);
  my $factors_by_slice = $ext_feat_adaptor->fetch_all_by_Slice_FeatureSets($slice, $fsets);

  my (%seen, @factors_to_return);

  foreach (@$factors_by_gene){
   my $label = $_->display_label .':'.  $_->start .''.$_->end;
   unless (exists $seen{$label}){
      push @factors_to_return, $_;
      $seen{$label} = 1;
   }
  }

  foreach (@$factors_by_slice){
   my $label = $_->display_label .':'. $_->start .''.$_->end;
   unless (exists $seen{$_->display_label}){
      push @factors_to_return, $_;
      $seen{$label} = 1;
   }
  }

 return \@factors_to_return;
}

sub reg_features {
  my $self = shift;
  my $gene = $self->gene;
  my $fg_db= $self->get_fg_db;
  my $slice =  $self->get_extended_reg_region_slice;
  my $reg_feat_adaptor = $fg_db->get_RegulatoryFeatureAdaptor;
  my $feats = $reg_feat_adaptor->fetch_all_by_Slice($slice);
  return $feats;

}


1;
