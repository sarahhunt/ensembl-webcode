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

package EnsEMBL::Web::Component::PhenotypeOntologyTerm;

use strict;

use base qw(EnsEMBL::Web::Component::Shared);


##cross reference to phenotype entries
sub phenotype_url{
  my $self  = shift;
  my $pheno = shift;
  my $pid   = shift;
 
  my $params = {
      'type'      => 'Phenotype',
      'action'    => 'Features',
      'ph'        => $pid,
      __clear     => 1
    };
  
  return sprintf('<a href="%s">%s</a>', $self->hub->url($params), $pheno );

}

##move to object
sub get_child_features{
  my $self = shift;
  my $type = shift;

  my $adaptor = $self->hub->get_databases('go')->{'go'}->get_OntologyTermAdaptor;

  my $ontologyterm = $adaptor->fetch_by_accession($self->hub->param('oa'));
  return undef unless defined $ontologyterm ;

  my $terms = $adaptor->fetch_all_by_parent_term( $ontologyterm );
  return undef unless $terms;

  my @pheno_feats;
  foreach my $term (@{$terms}){
    push @pheno_feats, @{$self->get_all_PhenotypeFeatures_by_type($term->accession(), $type)};
  }
  return @pheno_feats ?  \@pheno_feats : undef;
}

## move to object
sub get_all_PhenotypeFeatures_by_type{
  my $self = shift;
  my $acc  = shift;
  my $type = shift;

  my $vardb   = $self->hub->database('variation');
  ## put convenience method in adaptor for speed?
  my $phen_ad  = $vardb->get_adaptor('Phenotype');
  my $ps       = $phen_ad->fetch_all_by_ontology_accession($acc);


  my $phenf_ad = $vardb->get_adaptor('PhenotypeFeature');
  my @pheno_feats;
  foreach my $p(@{$ps}){
    my $pfs = $phenf_ad->fetch_all_by_phenotype_id_feature_type($p->dbID, $type );
    push @pheno_feats, @{$pfs} if @{$pfs} ;
  }
  return \@pheno_feats ;
}



1;
