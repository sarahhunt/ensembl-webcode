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

package EnsEMBL::Web::Component::PhenotypeOntologyTerm::Variants;



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

  my $phenotype_features = $self->SUPER::get_all_PhenotypeFeatures_by_type($self->hub->param('oa'),  'Variation');

  my $html;

  ## look for PhenotypeFeatures attached to this term
  if( @$phenotype_features){
    $html .= $self->format_table($phenotype_features);
  }
  else{
    $html .= '<p>No variants annotations found for ' .$self->hub->param('oa');
  }


  ## look for PhenotypeFeatures attached to child terms
  my $child_features = $self->SUPER::get_child_features('Variation');

  if( $child_features ){
    $html .= '<p><h3>Variant annotations for child terms of ' .$self->hub->param('oa') .'</h3><p>';
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
  #return "<p><b>Found $total features - too many to display<p></b>" if scalar(@{$phenotype_features}) > 5000;

  ## Group gene <-> phenotype combinations yucky structure
  foreach my $pf (@{$phenotype_features}){
    $data->{$pf->phenotype()->description()}{$pf->object_id()}{loc} = $pf->seq_region_name() .":".$pf->seq_region_start() ."-". $pf->seq_region_end() . "(" .$pf->seq_region_strand() . ")";
    push @{$data->{$pf->phenotype()->description()}{$pf->object_id()}{source}},  $pf->source_name();
    push @{$data->{$pf->phenotype()->description()}{$pf->object_id()}{studies}},  $pf->study_name();
#$data->{$pf->phenotype()->description()}{$pf->object_id()}{p_value} = [];
    push @{$data->{$pf->phenotype()->description()}{$pf->object_id()}{p_value}},  $pf->p_value();
    push @{$data->{$pf->phenotype()->description()}{$pf->object_id()}{risk_allele}},  $pf->risk_allele() if $pf->risk_allele() =~ /\D+/;
 
    $pid->{$pf->phenotype()->description()} = $pf->phenotype()->dbID();
  }

  my $table = $self->new_table([], [], { data_table => 1, sorting => [ 'names' ] });
  $table->add_columns(
    { key => 'desc',       'label' => 'Phenotype Description',    'title' => 'The phenotype description supplied with the annotation', 'sort' => 'html'},
    { key => 'name',       'label' => 'Variant Name',             'title' => 'The variant ID supplied in the annotations', 'sort' => 'html'},
    { key => 'loc',        'label' => 'Genomic location (strand)','title' => 'Position of the feature (e.g. chromosome number, start and end coordinates, forward or reverse strand)', 'sort' => 'position_html
'},
    { key => 'source',     'label' => 'Annotation source',        'title' => 'Project or database reporting the association', 'sort' => 'html'},
    { key => 'risk_allele', 'label' => 'Risk Allele',             'title' => 'Allele reported to be associated with the trait', 'sort' => 'html'},
    { key => 'studies',    'label' => 'Study',                    'title' => 'Link to the pubmed article or other source showing the association', 'sort' => 'html'},
    { key => 'p_value',    'label' => 'P value (negative log)',   'title' => 'The probability that the association is significant (a higher number indicates a higher probability)'},

  );

  foreach my $pheno( keys %{$data}){
    my $pheno_link = $self->SUPER::phenotype_url( $pheno, $pid->{$pheno}); ## find better way..
    foreach my $feat( keys %{$data->{$pheno}}){


      my $feat_link = $self->variation_url( $feat);
      $table->add_rows( {  desc         => $pheno_link,
                           name         => $feat_link,
                           loc          => $data->{$pheno}->{$feat}->{loc},
                           source       => unique( $data->{$pheno}->{$feat}->{source} ),
                           risk_allele  => unique( $data->{$pheno}->{$feat}->{risk_allele}), 
                           p_value      => unique( $data->{$pheno}->{$feat}->{p_value}),
                           studies      => unique( $data->{$pheno}->{$feat}->{studies}), 
                         });
     }
  } 
  return $table->render;
}

## only mention same source etc once
sub unique{

  my $data = shift;
  return "-" unless ref($data) eq 'ARRAY' && defined $data->[0];

  my %a;
  map { $a{$_} = 1; } @{$data};
  return join(", ", sort keys %a);
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

