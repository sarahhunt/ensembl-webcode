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

package EnsEMBL::Web::Object::PhenotypeOntologyTerm;
### NAME: EnsEMBL::Web::Object::PhenotypeOntologyTerm
### Wrapper around a Bio::EnsEMBL::OntologyTerm object

use strict;
use warnings;
no warnings "uninitialized";

use EnsEMBL::Web::Cache;
use HTML::Entities qw(encode_entities);

use base qw(EnsEMBL::Web::Object);
use Data::Dumper;

sub default_action { return 'Summary'; }

sub short_caption {
  my $self = shift;
  my $name = $self->name;
  my $caption = 'Phenotype-based displays for ' . $self->name . ' ('. $self->hub->param('oa') .')';

  return $caption;
}

sub caption {
  my $self = shift;
  return $self->name;
}

sub long_caption {
  my $self = shift;
  return 'Annotations associated with '.$self->name();
}


sub Obj{
  my $self = shift;

  return unless $self->hub->param('oa');

  my $ontolterm  = $self->hub->database('ontology')->get_adaptor('OntologyTerm')->fetch_by_stable_id( $self->hub->param('oa') );

  return $ontolterm && $ontolterm->name() || '';
}



sub get_all_Phenotypes {
  my $self = shift;

  my $vardb = $self->hub->database('variation');
  my $pa    = $vardb->get_adaptor('Phenotype');
  return $pa->fetch_by_ontology_accession( $self->hub->param('oa') );
};




sub name{
  my $self = @_;
warn "In Object::PhenotypeOntologyTerm name\n";
#print Dumper $self->Obj;
  return $self->Obj && $self->Obj->name() || '';
}



sub get_all_PhenotypeFeatures{
  my $self = shift;
  my $acc  = shift;
  my $type = shift;

  my $vardb   = $self->hub->database('variation');
  ## put convenience method in adaptor for speed?
  my $phen_ad  = $vardb->get_adaptor('Phenotype');
  my $ps       = $phen_ad->fetch_by_ontology_accession($acc);


  my $phenf_ad = $vardb->get_adaptor('PhenotypeFeature');
  my @pheno_feats;
  foreach my $p(@{$ps}){
    my $pfs     = $pfa->fetch_all_by_phenotype_id_feature_type($p->dbID, $type );
    push @pheno_feats, $pfs ;
  }
  return \@pheno_feats;
}

sub availability {
  my $self = shift;
  
  if (!$self->{'_availability'}) {
    my $availability;

    $availability->{'phenotypeontologyterm'} = 1;

    my $counts = $self->counts;      
    $availability->{"has_$_"} = $counts->{$_} for qw(genes variation structuralvariation qtl); 
 #   print Dumper $availability;
    $self->{'_availability'} = $availability;
  }
  
  return $self->{'_availability'};
}

sub counts {
  my $self = shift;

  my $counts = $self->{'_counts'};

  unless ($counts) {

    my $vardb = $self->hub->database('variation');
    my $pa    = $vardb->get_adaptor('Phenotype');
    my $count_by_feat = $pfa->count_features_by_ontology_accession( $self->hub->param('oa') );

    $counts = {};
    $counts->{'genes'}              = $count_by_feat->{'Genes'}     || 0;
    $counts->{'variants'}           = $count_by_feat->{'Variation'} || 0;
    $counts->{'structuralvariants'} = $count_by_feat->{'StructuralVariation'}  +  $count_by_feat->{'SupportingStructuralVariation'} || 0;
    $counts->{'qtl'}                = $count_by_feat->{'QTL'}       || 0;

    $self->{'_counts'} = $counts;
  }

  return $counts;
}

1;
