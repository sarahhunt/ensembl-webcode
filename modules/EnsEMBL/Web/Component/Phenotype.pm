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

package EnsEMBL::Web::Component::Phenotype;

use strict;

use base qw(EnsEMBL::Web::Component::Shared);


##  Return a table of phenotype association records 
##  for an array of ontology accessions
##  Returns data at same or child level
sub get_phenotypes_mapping_to_same_terms{

  my $self       = shift;
  my $accessions = shift;
  my $is_child   = shift;

  my @phenos;
  my %names;
  my $adaptor = $self->hub->get_databases('go')->{'go'}->get_OntologyTermAdaptor;

  foreach my $acc (@{$accessions}){ 
    next if $acc =~/HP:0000006|HP:0000007/;# Autosomal dominant & Autosomal recessive
    my $ontologyterm = $adaptor->fetch_by_accession($acc);
    $names{$acc} = $ontologyterm->name() if defined $ontologyterm;
 
    my $other_phenos = $self->object->get_phenotype_by_ontology_accession($acc);
    push @phenos, @{$other_phenos}

  }
  return undef unless @phenos && scalar @phenos > 1;


  my $table = $self->new_table([], [], { data_table => 1, sorting => [ 'desc' ] });
   
  $table->add_columns(
   { key => 'desc',        'label' => 'Description',          'title' => 'Phenotype description supplied by data souce', 'sort' => 'html'}, 
   { key => 'accession',   'label' => 'Ontology Accession',   'title' => 'Accession of the ontology term the description has been mapped to', 'sort' => 'position_html'},
   { key => 'term',        'label' => 'Ontology Term',        'title' => 'Name of the ontology term the description has been mapped to', 'sort' => 'position_html'},
   { key => 'feat_count',  'label' => 'Feature Count',        'title' => 'Number of reported Variant, gene, or QTL', 'sort' => 'html'},
   { key => 'phe_sources', 'label' => 'Annotation source(s)', 'title' => 'Project or database reporting the associations', 'sort' => 'html'}
  );


  my $count_results = 0;
  my %done; ## if looking for match to phenotype description, may be linked by different accessions to the same phenotype

  foreach my $pheno (@phenos) {

    ## skip if the same
    next if defined $self->hub->param('ph') &&  $pheno->dbID() eq  $self->hub->param('ph');
    ## skip if already seen
    next if $done{$pheno->dbID()};
    $done{$pheno->dbID()} = 1;

    $count_results++;

    ## extract feature counts & sources
#    my $feat       = $self->hub->get_adaptor('get_PhenotypeFeatureAdaptor', 'variation')->fetch_all_by_phenotype_id_source_name($pheno->dbID());
#    my $feat_count = scalar(@{$feat});
#    my %source;
#    foreach my $pf(@{$feat}){
#      $source{$pf->source_name()} = 1;
#    }
 
    my $feat_count;
    my $counts  = $self->hub->get_adaptor('get_PhenotypeFeatureAdaptor', 'variation')->count_all_with_source_by_Phenotype($pheno);
    foreach my $s(keys %{$counts}){
      $feat_count +=  $counts->{$s};
    }

    my $sources    = join(", ", keys %{$counts});
    my @accessions = @{$pheno->ontology_accessions()};

    my $terms;
    foreach my $accn (@accessions){
      $terms .= $names{$accn} . ","
    }
    chop($terms);

    my $row = {
              'desc'        => {'value' => '<a href="Locations?ph='.$pheno->dbID() . '">' . $pheno->description() . '</a>' },
              'accession'   => {'value' => join(", ", @accessions) },
              'term'        => {'value' => $terms},
              'feat_count'  => {'value' => $feat_count},
              'phe_sources' => {'value' => $sources}
              };
     $table->add_rows($row);
  }

  ## Is there anything to report?
  return undef unless $count_results > 0;

   my $html =  ($is_child == 1 ?
    ('<h3>' . $count_results. ' phenotypes map to child ontology terms </h3>') :
    ('<h3>' . $count_results. ' phenotypes map to the same ontology terms </h3>'));

  $html .= $table->render;

  return $html;
}

sub get_phenotypes_mapping_to_child_terms{

  my $self       = shift;
  my $accessions = shift;

  my $adaptor = $self->hub->get_databases('go')->{'go'}->get_OntologyTermAdaptor;

  my @child_accessions;
  foreach my $acc (@{$accessions}){

    my $ontologyterm = $adaptor->fetch_by_accession($acc);
    next unless defined $ontologyterm ;

    my $terms = $adaptor->fetch_all_by_parent_term( $ontologyterm );
    next unless $terms;

    foreach my $term (@{$terms}){
      push @child_accessions, $term->accession();
    }
  }
  return $self->get_phenotypes_mapping_to_same_terms(\@child_accessions, 1);

}

1;
