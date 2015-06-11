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

package EnsEMBL::Web::Object::Gene;

### Wrapper around a Bio::EnsEMBL::Gene object

use parent qw(EnsEMBL::Web::Object);

######## WRAPPERS AROUND API METHODS #############

sub gene                        { return $_[0]->api_object;             }
sub stable_id                   { return $_[0]->api_object->stable_id;  }
sub feature_type                { return $_[0]->api_object->type;       }
sub analysis                    { return $_[0]->api_object->analysis;   }
sub source                      { return $_[0]->api_object->source;     }
sub version                     { return $_[0]->api_object->version;    }
sub logic_name                  { return $_[0]->api_object->analysis->logic_name; }
sub coord_system                { return $_[0]->api_object->slice->coord_system->name; }
sub seq_region_type             { return $_[0]->coord_system;    }
sub seq_region_name             { return $_[0]->api_object->slice->seq_region_name; }
sub seq_region_start            { return $_[0]->api_object->start;      }
sub seq_region_end              { return $_[0]->api_object->end;        }
sub seq_region_strand           { return $_[0]->api_object->strand;     }
sub feature_length              { return $_[0]->api_object->feature_Slice->length; }
sub get_latest_incarnation      { return $_[0]->api_object->get_latest_incarnation; }
sub get_all_associated_archived { return $_[0]->api_object->get_all_associated_archived; }

######### METHODS USEFUL FOR TABS, ETC ################

sub type_name                   { return $_[0]->species_defs->translate('Gene'); }

sub default_action { return $_[0]->api_object->isa('Bio::EnsEMBL::ArchiveStableId') ? 'Idhistory' : $_[0]->api_object->isa('Bio::EnsEMBL::Compara::Family') ? 'Family' : 'Summary'; }

sub short_caption {
## Create short label used on local navigation (tabs and lefthand menu)
  my $self = shift;
  
  return 'Gene-based displays' unless shift eq 'global';
  
  my $dxr   = $self->api_object->can('display_xref') ? $self->api_object->display_xref : undef;
  my $label = $dxr ? $dxr->display_id : $self->api_object->stable_id;
  
  return "Gene: $label";  
}

sub caption {
## Create long label used in page headings
  my $self = shift;
  my $heading = $self->type_name.': ';
  my $subhead;

  my( $disp_id ) = $self->display_xref;
  if( $disp_id && $disp_id ne $self->stable_id ) {
    $heading .= $disp_id;
    $subhead = $self->stable_id;
  }
  else {
    $heading .= $self->stable_id;
  }

  return [$heading, $subhead];
}

sub display_xref {
  my $self = shift;
  return undef if $self->Obj->isa('Bio::EnsEMBL::Compara::Family');
  return undef if $self->Obj->isa('Bio::EnsEMBL::ArchiveStableId');
  my $trans_xref = $self->Obj->display_xref();
  return undef unless  $trans_xref;
  (my $db_display_name = $trans_xref->db_display_name) =~ s/(.*HGNC).*/$1 Symbol/; #hack for HGNC name labelling, remove in e58
  return ($trans_xref->display_id, $trans_xref->dbname, $trans_xref->primary_id, $db_display_name, $trans_xref->info_text );
}


1
