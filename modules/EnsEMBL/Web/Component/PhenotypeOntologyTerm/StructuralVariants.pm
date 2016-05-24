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

package EnsEMBL::Web::Component::PhenotypeOntologyTerm::StructuralVariants;



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

  my $phenotype_features = $self->SUPER::get_all_PhenotypeFeatures_by_type($self->hub->param('oa'),  'StructuralVariation');


  my $html;

  ## look for PhenotypeFeatures attached to this term
  if( @$phenotype_features){
    $html .= $self->format_table($phenotype_features);
  }
  else{
    $html .= '<p>No structural variant annotations found for ' .$self->hub->param('oa');
  }


  ## look for PhenotypeFeatures attached to child terms
  my $child_features = $self->SUPER::get_child_features('StructuralVariation');

  if( $child_features ){
    $html .= '<p><h3>Vtructural variant annotations for child terms of ' .$self->hub->param('oa') .'</h3><p>';
    $html .= $self->format_table($child_features);
  }

  return $html;
}


sub format_table{

  my $self = shift;
  my $phenotype_features = shift;

  my $data;
  my $pid;

  my $total = scalar(@{$phenotype_features});
  return "<p><b>Found $total features - too many to display<p></b>" if scalar(@{$phenotype_features}) > 1000;

  ## Group gene <-> phenotype combinations yucky structure
  foreach my $pf (@{$phenotype_features}){
    $data->{$pf->phenotype()->description()}{$pf->object_id()}{loc} = $pf->seq_region_name() .":".$pf->seq_region_start() ."-". $pf->seq_region_end() . "(" .$pf->seq_region_strand() . ")";
    push @{$data->{$pf->phenotype()->description()}{$pf->object_id()}{source}},  $pf->source_name();
    push @{$data->{$pf->phenotype()->description()}{$pf->object_id()}{studies}},  $pf->study_name();
 
    $pid->{$pf->phenotype()->description()} = $pf->phenotype()->dbID();
  }

  my $table = $self->new_table([], [], { data_table => 1, sorting => [ 'names' ] });
  $table->add_columns(
    { key => 'desc',       'label' => 'Phenotype Description',    'title' => 'The phenotype description supplied with the annotation', 'sort' => 'html'},
    { key => 'name',       'label' => 'Variant Name',             'title' => 'The variant ID supplied in the annotations', 'sort' => 'html'},
    { key => 'loc',        'label' => 'Genomic location (strand)','title' => 'Position of the feature (e.g. chromosome number, start and end coordinates, forward or reverse strand)', 'sort' => 'position_html
'},
    { key => 'source',     'label' => 'Annotation source',        'title' => 'Project or database reporting the association', 'sort' => 'html'},

  );

  foreach my $pheno( keys %{$data}){
    my $pheno_link = $self->SUPER::phenotype_url( $pheno, $pid->{$pheno}); ## find better way..
    foreach my $feat( keys %{$data->{$pheno}}){


      my $feat_link = $self->variation_url( $feat);
      $table->add_rows( {  desc         => $pheno_link,
                           name         => $feat_link,
                           loc          => $data->{$pheno}->{$feat}->{loc},
                           source       => join( ", ", @{$data->{$pheno}->{$feat}->{source}}), 
                         });
    }
  } 


  return $table->render;
}

sub variation_url{
  my $self = shift;
  my $var  = shift;

  my $params = {
      'type'   => 'Variation',
      'action' => 'Phenotype',
      'v'      => $var,
      __clear     => 1
    };

  return sprintf('<a href="%s">%s</a>', $self->hub->url($params), $var);
}


1;

