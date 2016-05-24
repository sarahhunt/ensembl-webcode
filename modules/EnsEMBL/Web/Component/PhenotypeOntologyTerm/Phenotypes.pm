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

package EnsEMBL::Web::Component::PhenotypeOntologyTerm::Phenotypes;



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


  my $html .= '<h3>Records linked to ontology term '. $ontology_accession .'</h3>';
  $html .= '<p>';
  $html .= $self->count_features([$ontology_accession]);

  $html .= '<p>';
  $html .= '<h3>Records linked to child terms of '. $ontology_accession .'</h3>';
  $html .= $self->count_child_features();


  ## Add phenotype Features
#  $html .= $self->get_features();

  return $html;
}

sub count_child_features{
  my $self      = shift;

  my @data;

  my $adaptor = $self->hub->get_databases('go')->{'go'}->get_OntologyTermAdaptor;

  my $ontologyterm = $adaptor->fetch_by_accession($self->hub->param('oa'));
  return undef unless $ontologyterm;

  my $terms = $adaptor->fetch_all_by_parent_term( $ontologyterm );
  return undef unless $terms;

  my @accessions;
  foreach my $term (@{$terms}){
    push @accessions, $term->accession();
  }

  return $self->count_features( \@accessions);

}

sub count_features{

  my $self       = shift;
  my $accessions = shift;

  my @rows;

  my $vardb   = $self->hub->database('variation');
  my $phen_ad = $vardb->get_adaptor('Phenotype');
  my $pf_ad   = $vardb->get_adaptor('PhenotypeFeature');

  my $total_annotations = 0; ## keep an overall count

  foreach my $accession (@{$accessions}){
    my $phenotypes = $phen_ad->fetch_all_by_ontology_accession( $accession );
    foreach my $pheno(@{$phenotypes}){
      my $number_of_features = $pf_ad->count_all_by_phenotype_id($pheno->dbID());
      $total_annotations += $number_of_features;
      push @rows, [ $accession, $self->SUPER::phenotype_url($pheno->description,$pheno->dbID()), $number_of_features];
    }
  }

  return '<p>No annotations available'  unless @rows;

  ## format
  my $table = $self->new_table();
  $table->add_columns( 
    { key => 'accession',   title => 'Accession',   align => 'left', sort => 'html'    },
    { key => 'description', title => 'Description', align => 'left', sort => 'html'    },  
    { key => 'features',    title => 'Annotations', align => 'left', sort => 'html'    },
  );

  foreach my $row (@rows){  $table->add_rows($row);}

  return $table->render;
}





sub get_features{

  my $self = shift;
  
  my $phenos  = $self->hub->get_adaptor('get_PhenotypeAdaptor', 'variation')->fetch_all_by_ontology_accession($self->hub->param('oa'));


  my $table = $self->new_table([], [], { data_table => 1, sorting => [ 'desc' ] });
  $table->add_columns(
   { key => 'desc',        'label' => 'Description',          'title' => 'Phenotype description supplied by data souce', 'sort' => 'html'},
   { key => 'accessions',  'label' => 'Ontology Accessions',  'title' => 'Ontology term the description has been mapped to', 'sort' => 'position_html'},
   { key => 'feat_count',  'label' => 'Feature Count',        'title' => 'Number of reported Variant, gene, or QTL', 'sort' => 'html'},
   { key => 'phe_sources', 'label' => 'Annotation source(s)', 'title' => 'Project or database reporting the associations', 'sort' => 'html'}
  );

  my %full_data;
  foreach my $pheno (@{$phenos}) {
    my %source; ## store number of sources using each description
    my $feat       = $self->hub->get_adaptor('get_PhenotypeFeatureAdaptor', 'variation')->fetch_all_by_phenotype_id_source_name($pheno->dbID());
    my $feat_count = scalar(@{$feat});

    foreach my $pf(@{$feat}){
      $source{$pf->source_name()} = 1;
      next if $pf->seq_region_name =~/HSC/;  ## HACK to remove ALTs
      ## save full data set

     # my $variation_url   = $hub->url({ action => "Variation", action => "Phenotype", vf => $vf });
#      my $phenotype_url   = $hub->url({ action => "Phenotype", action => "Summary", pf => $pheno->dbID() });
#warn "Got phenot URL :  $phenotype_url\n";
      push @{$full_data{$pf->source_name()}}, { 'desc'      => {'value' =>  $pheno->description() },
                                                'name'      => {'value' =>  $pf->object_id},
                                                'loc'       => {'value' =>  $pf->seq_region_name.':'.$pf->seq_region_start .'-'.$pf->seq_region_end }, 
                                                'feat_type' => {'value' =>  $pf->type}, 
                                                'source'    => {'value' =>  $pf->source_name()} }; 
    }
    my $sources = join(", ", keys %source);
    my $row = {
              'desc'        => {'value' => '<a href="Locations?ph='.$pheno->dbID() . '">' . $pheno->description() . '</a>' },
              'accessions'  => {'value' => join(", ", @{$pheno->ontology_accessions()} ) },
              'feat_count'  => {'value' => $feat_count},
              'phe_sources' => {'value' => $sources}
              };
     $table->add_rows($row);
  }

  return '<b> No annotations available for this phenotype<\b>' unless scalar(keys %full_data) > 0;

  my $html .= $table->render;
  $html .= $self->full_data(\%full_data);

  return $html;
}

1;

