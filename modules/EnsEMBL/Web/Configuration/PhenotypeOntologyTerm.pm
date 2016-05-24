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


package EnsEMBL::Web::Configuration::PhenotypeOntologyTerm;

use strict;

use base qw(EnsEMBL::Web::Configuration);
use Data::Dumper;
sub caption { 
  my $self = shift;
  my $type = ($self->hub->param('oa') =~/HP/           ? 'Phenotype ': 
              $self->hub->param('oa') =~ /Orphanet|DO/ ? 'Disease '  : 
              '');
  return $type . 'Ontology Accession ' . $self->hub->param('oa'); 

}

sub modify_page_elements { $_[0]->page->remove_body_element('summary'); }

sub set_default_action {
  my $self = shift;
  $self->{'_data'}->{'default'} = 'Summary'; 
}

sub tree_cache_key {
  my $self = shift;

  my $desc = $self->object ? $self->object->Obj->name : '';
  return join '::', $self->SUPER::tree_cache_key(@_), $desc;
}

sub populate_tree {
  my $self = shift;
  my $hub  = $self->hub;

  $self->create_node('Summary', 'Summary',
    [qw( phenotypes EnsEMBL::Web::Component::PhenotypeOntologyTerm::Phenotypes  )],
  );
#  print "<p> In pop tree<p>";
  #my $avail = $self->object->availability();

  ## availability!
  $self->create_node('Genes', 'Annotated Genes',
    [qw( genes EnsEMBL::Web::Component::PhenotypeOntologyTerm::Genes )],
#     { 'availability' => 'phenotypeontologyterm has_genes', 'concise' => 'Genes' }
  );
  ## availability!
  $self->create_node('Variants', 'Annotated Variants',
    [qw( variants EnsEMBL::Web::Component::PhenotypeOntologyTerm::Variants )],
#     { 'availability' => 'phenotypeontologyterm has_variants', 'concise' => 'Variant' }
  );

  ## availability!
  $self->create_node('StructuralVariants', 'Annotated Structural Variants',
    [qw( structuralvariants EnsEMBL::Web::Component::PhenotypeOntologyTerm::StructuralVariants )],
#     { 'availability' => 'phenotypeontologyterm has_structuralvariants', 'concise' => 'StructuralVariant' }
  );

  ## availability! No data for human
  $self->create_node('QTL', 'Annotated QTL',
    [qw( QTL EnsEMBL::Web::Component::PhenotypeOntologyTerm::QTL )],
     { 'availability' => 0}
  );

  ## availability!
#  $self->create_node('Citations', 'Citations',
#    [qw(citations EnsEMBL::Web::Component::PhenotypeOntologyTerm::Citations )],
#     { 'availability' => 0}
#  );






}

1;

