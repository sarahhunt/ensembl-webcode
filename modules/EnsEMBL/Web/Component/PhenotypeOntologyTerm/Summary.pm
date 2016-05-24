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

package EnsEMBL::Web::Component::PhenotypeOntologyTerm::Summary;



use strict;

use HTML::Entities qw(encode_entities);

use EnsEMBL::Web::Controller::SSI;
use EnsEMBL::Web::Exceptions;

use base qw(EnsEMBL::Web::Component::PhenotypeOntologyTerm);
use Data::Dumper;
sub _init {
  my $self = shift;
  $self->cacheable(0);
  $self->ajaxable(0);
}

sub content {
  my $self      = shift;

  #return unless  $self->hub->param('oa');
  my $ontology_accession = $self->hub->param('oa');

  ## look up term info
  my $adaptor = $self->hub->get_databases('go')->{'go'}->get_OntologyTermAdaptor;
  my $ontologyterm = $adaptor->fetch_by_accession($ontology_accession);

  my $external_link = $self->ext_link();
  my $html;

  if(defined $ontologyterm ){
    my $summary_table = $self->new_twocol( [ 'Term',       '<b>'. $ontologyterm->name().'</b>'],
                                           [ 'Accession',  $external_link],
                                           [ 'Definition', (split/\"/,$ontologyterm->definition())[1]]
                                         );
    $html .= $summary_table->render;
  }
  else{
    $html .= '<b>View '. $external_link . 'at source</b><p>';
  }  


  return $html;
}

# put somewhere central
sub ext_link{
  my $self = shift;

  my $iri_form = $self->hub->param('oa');
  $iri_form =~ s/\:/\_/ unless $iri_form  =~ /^GO/;;

  my $ontology_link;
  $ontology_link = $self->hub->get_ExtURL_link($self->hub->param('oa'), 'EFO',  $iri_form) if $iri_form =~ /^EFO/;
  $ontology_link = $self->hub->get_ExtURL_link($self->hub->param('oa'), 'ORDO', $iri_form) if $iri_form =~ /^Orphanet/;
  $ontology_link = $self->hub->get_ExtURL_link($self->hub->param('oa'), 'DOID', $iri_form) if $iri_form =~ /^DO/;
  $ontology_link = $self->hub->get_ExtURL_link($self->hub->param('oa'), 'HPO',  $iri_form) if $iri_form =~ /^HP/;
  $ontology_link = $self->hub->get_ExtURL_link($self->hub->param('oa'), 'GO',   $iri_form) if $iri_form =~ /^GO/;


  return $ontology_link;
}

1;
