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

## Creates table of ontology terms which map to this phenotype description
## Display on all but list of phenotypes tab
package EnsEMBL::Web::Component::Phenotype::OntologyTerm;



use strict;

use HTML::Entities qw(encode_entities);

use EnsEMBL::Web::Controller::SSI;
use EnsEMBL::Web::Exceptions;

use base qw(EnsEMBL::Web::Component::Phenotype);

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

  return $html;
}


sub ontology_mappings{

  my $self = shift;

  ## Find matching ontology terms
  my $ontol_data =  $self->get_all_ontology_data();

  return undef unless $ontol_data;

  my $html;

  my @rows;
  foreach my $acc (keys %{$ontol_data}){
    push @rows, [$ontol_data->{$acc}->{name},  $ontol_data->{$acc}->{link}, $ontol_data->{$acc}->{definition} ];
  }

  my $table = $self->new_table([], [], { data_table => 1 });
  $table->add_columns(
    { key => 'term',       title => 'Term',       align => 'left', sort => 'html'    },
    { key => 'accession',  title => 'Accession',  align => 'left', sort => 'html'    },
    { key => 'definition', title => 'Definition', align => 'left', sort => 'html'    },
  );

  foreach my $row (@rows){  $table->add_rows($row);}



  if(scalar(@rows) < 6){
    $html = "<h3> Ontology mappings: </h3><p>";
    $html .=  $table->render;
    $html .= "<p>";

  }
  else{
    my $title = 'Ontology mappings:';
    my $id = 'ontology_mappings';
    $html .= $self->toggleable_table($title, $id, $table, 0);
  }
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

      $data{$acc}{link}       = $self->external_ontology_link($acc);
      $data{$acc}{name}       = $ot->name();
      $data{$acc}{definition} = (split/\"/, $ot->definition())[1]; 
    }
  }

  ## temp - for terms not in the ontology db
  my $accessions = $self->object->pheno->ontology_accessions();
  foreach my $acc (@{$accessions}){
    $data{$acc}{link}  = $self->external_ontology_link($acc)
      unless $data{$acc}{link};
  }
  return \%data;
}
## build link out to Ontology source
sub external_ontology_link{
  my $self = shift;
  my $acc  = shift;

  my $iri_form = $acc;
  $iri_form =~ s/\:/\_/ unless $acc =~ /^GO/;

  my $ontology_link;
  $ontology_link = $self->hub->get_ExtURL_link( $acc, 'EFO',  $iri_form) if $iri_form =~ /^EFO/;
  $ontology_link = $self->hub->get_ExtURL_link( $acc, 'ORDO', $iri_form) if $iri_form =~ /^Orphanet/;
  $ontology_link = $self->hub->get_ExtURL_link( $acc, 'HPO',  $iri_form) if $iri_form =~ /^HP/;
  $ontology_link = $self->hub->get_ExtURL_link( $acc, 'GO',   $iri_form) if $iri_form =~ /^GO/;

  return $ontology_link;
}

1;

