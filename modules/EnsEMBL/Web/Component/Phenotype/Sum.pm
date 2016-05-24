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

## Don't want this to display with list of Phenotypes
package EnsEMBL::Web::Component::Phenotype::Sum;



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
  my $self = shift;

  return unless $self->hub->param('ph');

  ## create table of ontology terms & definitions
  my $html = $self->ontology_mappings();

  ## summary of available features & sources
  $html .= $self->feature_counts();

  return $html;
}

sub feature_counts{

  my $self = shift;

  my $ph_id     = $self->hub->param('ph');

  ## Count features annotated with this phenotype
  my $feats = $self->hub->get_adaptor('get_PhenotypeFeatureAdaptor', 'variation')->fetch_all_by_phenotype_id_source_name($ph_id);
  my %count_source;
  my $tot;
  foreach my $pf(@{$feats}){
    $tot++;
    $count_source{$pf->source_name()}{$pf->type()}++;
  }
  my $html = "<h3>This phenotype is linked to $tot features</h3>";

  my @features;
  foreach my $source( keys %count_source){
    foreach my $type (keys %{$count_source{$source}}){ 
      push @features, ["$source $type", $count_source{$source}{$type}];
    }
  }
  ## links needed
  my $features_table = $self->new_twocol( @features );
  $html .= $features_table->render; 
  $html .= "<p>";

  return $html; 
}


sub ontology_mappings{

  my $self = shift;

  ## Find matching ontology terms
  my $ontol_data =  $self->get_all_ontology_data();

  return '<b> No mappings to ontologies available in current Ensembl database</b>' unless $ontol_data;

  my $html = "<h3> Ontology mappings: </h3><p>";

  my @rows;
  foreach my $acc (keys %{$ontol_data}){
    push @rows, [$ontol_data->{$acc}->{name},  $ontol_data->{$acc}->{link}, $ontol_data->{$acc}->{definition}, $ontol_data->{$acc}->{external} ];
 }

  my $table = $self->new_table();
  $table->add_columns(
    { key => 'term',       title => 'Term',       align => 'left', sort => 'html'    },
    { key => 'accession',  title => 'Accession',  align => 'left', sort => 'html'    },
    { key => 'definition', title => 'Definition', align => 'left', sort => 'html'    },
    { key => 'external',   title => 'Ontology',   align => 'left', sort => 'html'    },
  );

  foreach my $row (@rows){  $table->add_rows($row);}
  $html .=  $table->render;
  $html .= "<p>";

  return $html;
}


## get ontology information from object & format
sub get_all_ontology_data{

  my $self = shift;

  my %data;

  my $ontologyterms = $self->object->get_OntologyTerms();

  if ($ontologyterms){

    foreach my $ot (@{$ontologyterms}){

      ## build link out to Ontology source
      my $acc = $ot->accession();
      next unless $acc; 

      $data{$acc}{link}       = $self->SUPER::pheno_ont_url($acc);
      $data{$acc}{name}       = $ot->name();
      $data{$acc}{definition} = (split/\"/, $ot->definition())[1]; 
      $data{$acc}{external}   = $self->external_ontology_link($acc);
    }
  }

  ## temp -remove before release -  DOID and HPO not in current ontology DB
  my $accessions = $self->object->pheno->ontology_accessions();
  foreach my $acc (@{$accessions}){
    unless ($data{$acc}{link}){ 
    $data{$acc}{link}     = $self->SUPER::pheno_ont_url($acc);
    $data{$acc}{external} = $self->external_ontology_link($acc);
    }
  }
  ## temp -remove before release 

  return \%data;
}

## build link out to Ontology source
sub external_ontology_link{

  my $self = shift;
  my $acc  = shift;

  my $iri_form = $acc;
  $iri_form =~ s/\:/\_/ unless $acc =~ /^GO/;

  my $ontology_link;
  $ontology_link = $self->hub->get_ExtURL_link('EFO',   'EFO',  $iri_form)    if $iri_form =~ /^EFO/;
  $ontology_link = $self->hub->get_ExtURL_link('Orphanet', 'ORDO', $iri_form) if $iri_form =~ /^Orphanet/;
  $ontology_link = $self->hub->get_ExtURL_link('DO',   'DOID', $iri_form)     if $iri_form =~ /^DO/;
  $ontology_link = $self->hub->get_ExtURL_link('HP',   'HPO',  $iri_form)     if $iri_form =~ /^HP/;
  $ontology_link = $self->hub->get_ExtURL_link('GO',   'GO',  $iri_form)      if $iri_form =~ /^GO/;

  return $ontology_link;

}
1;
