=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package EnsEMBL::Web::Configuration::Phenotype;

use strict;

use base qw(EnsEMBL::Web::Configuration);

sub caption { return 'Phenotype'; }

sub modify_page_elements { $_[0]->page->remove_body_element('summary'); }

sub set_default_action {
  my $self = shift;
  $self->{'_data'}->{'default'} = 'Locations'; 
}

sub tree_cache_key {
  my $self = shift;
  my $desc = $self->object ? $self->object->get_phenotype_desc : 'All phenotypes';
  return join '::', $self->SUPER::tree_cache_key(@_), $desc;
}

sub populate_tree {
  my $self = shift;
  my $hub  = $self->hub;

  $self->create_node('All', 'List of Phenotypes', [qw(all_phenotypes EnsEMBL::Web::Component::Phenotype::All )] );

  my $avail = ($self->object && $self->object->phenotype_id) ? 1 : 0;

  my $title = $self->object ? $self->object->long_caption : '';

  ## This can only be drawn for limited numbers of features
  $self->create_node('Karyotype', "Karyotype View",
    [qw(sum EnsEMBL::Web::Component::Phenotype::Sum  locations EnsEMBL::Web::Component::Phenotype::Locations )],
    { 'availability' => $avail, 'concise' => $title },
  );

  $self->create_node('Features', "Feature View",
    [qw(sum EnsEMBL::Web::Component::Phenotype::Sum  locations EnsEMBL::Web::Component::Phenotype::Features )],
    { 'availability' => $avail, 'concise' => $title },
  );



  my $ot_avail = ($self->object->pheno && $self->object->pheno->ontology_accessions) ? 1 : 0;
  $self->create_node('CoMapped', 'Similar Phenotypes',
    [qw( sum EnsEMBL::Web::Component::Phenotype::Sum  
         comapped EnsEMBL::Web::Component::Phenotype::MatchingPhenotypes  )],
    { 'availability' => $ot_avail, 'concise' => 'Matching Phenotypes' }
  );


}

1;
