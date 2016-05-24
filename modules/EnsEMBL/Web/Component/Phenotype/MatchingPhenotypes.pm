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

package EnsEMBL::Web::Component::Phenotype::MatchingPhenotypes;



use strict;

use HTML::Entities qw(encode_entities);

use EnsEMBL::Web::Controller::SSI;
use EnsEMBL::Web::Exceptions;

use base qw(EnsEMBL::Web::Component::Phenotype);
use Data::Dumper;
sub _init {
  my $self = shift;
  $self->cacheable(0);
  $self->ajaxable(0);
}

sub content {
  my $self      = shift;



#this doesn't work:
#  my $pheno = $self->object->Obj;
  my $pheno = $self->object->pheno;

  my $html;

  foreach my $acc (@{$pheno->ontology_accessions()}){ 
    my $link = $self->SUPER::pheno_ont_url($acc);
    $html .= '<d><b>View features for ' . $link . '</b>';
  }
  $html = '<p>';

  ## get summary table with links to other phenos
  $html .= $self->SUPER::get_phenotypes_mapping_to_same_terms( $pheno->ontology_accessions() );  

  ## get summary table with links to child phenos
  $html .= $self->SUPER::get_phenotypes_mapping_to_child_terms( $pheno->ontology_accessions() );


  return $html; 
 
}
 
1;
