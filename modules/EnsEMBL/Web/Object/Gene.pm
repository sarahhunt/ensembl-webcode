=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package EnsEMBL::Web::Object::Gene;

### Wrapper around a Bio::EnsEMBL::Gene object

use strict;

use Time::HiRes qw(time);

use EnsEMBL::Web::Constants; 
use EnsEMBL::Web::Cache;
use Bio::EnsEMBL::Compara::Homology;

use parent qw(EnsEMBL::Web::Object);

our $MEMD = EnsEMBL::Web::Cache->new;

######## WRAPPERS AROUND API METHODS #############

sub gene                        { return $_[0]->Obj;             }
sub stable_id                   { return $_[0]->Obj->stable_id;  }
sub feature_type                { return $_[0]->Obj->type;       }
sub analysis                    { return $_[0]->Obj->analysis;   }
sub source                      { return $_[0]->Obj->source;     }
sub version                     { return $_[0]->Obj->version;    }
sub logic_name                  { return $_[0]->Obj->analysis->logic_name; }
sub coord_system                { return $_[0]->Obj->slice->coord_system->name; }
sub seq_region_type             { return $_[0]->coord_system;    }
sub seq_region_name             { return $_[0]->Obj->slice->seq_region_name; }
sub seq_region_start            { return $_[0]->Obj->start;      }
sub seq_region_end              { return $_[0]->Obj->end;        }
sub seq_region_strand           { return $_[0]->Obj->strand;     }
sub feature_length              { return $_[0]->Obj->feature_Slice->length; }
sub get_latest_incarnation      { return $_[0]->Obj->get_latest_incarnation; }
sub get_all_associated_archived { return $_[0]->Obj->get_all_associated_archived; }
sub gxa_check                   { return; } #implemented in widget plugin, to check for gene expression atlas availability

######### WEB-SPECIFIC METHODS ################

sub type_name                   { return $_[0]->species_defs->translate('Gene'); }

sub default_action { return $_[0]->Obj->isa('Bio::EnsEMBL::ArchiveStableId') ? 'Idhistory' : $_[0]->Obj->isa('Bio::EnsEMBL::Compara::Family') ? 'Family' : 'Summary'; }

sub short_caption {
  my $self = shift;
  
  return 'Gene-based displays' unless shift eq 'global';
  
  my $dxr   = $self->Obj->can('display_xref') ? $self->Obj->display_xref : undef;
  my $label = $dxr ? $dxr->display_id : $self->Obj->stable_id;
  
  return "Gene: $label";  
}

sub caption {
  my $self = shift;
  my $heading = $self->type_name.': ';
  my $subhead;

  my( $disp_id ) = $self->display_xref;
  if( $disp_id && $disp_id ne $self->stable_id ) {
    $heading .= $disp_id;
    $subhead = $self->stable_id;
  }
  else {
    $heading .= $self->stable_id;
  }

  return [$heading, $subhead];
}

sub availability {
  my $self = shift;
  my ($database_synonym) = @_;
  
  if (!$self->{'_availability'}) {
    my $availability = $self->_availability;
    my $obj = $self->Obj;
    
    if ($obj->isa('Bio::EnsEMBL::ArchiveStableId')) {
      $availability->{'history'} = 1;
    } elsif ($obj->isa('Bio::EnsEMBL::Gene')) {
      my $member      = $self->database('compara') ? $self->database('compara')->get_GeneMemberAdaptor->fetch_by_stable_id($obj->stable_id) : undef;
      my $pan_member  = $self->database('compara_pan_ensembl') ? $self->database('compara_pan_ensembl')->get_GeneMemberAdaptor->fetch_by_stable_id($obj->stable_id) : undef;
      my $counts      = $self->counts($member, $pan_member);
      my $rows        = $self->table_info($self->get_db, 'stable_id_event')->{'rows'};
      my $funcgen_res = $self->database('funcgen') ? $self->table_info('funcgen', 'feature_set')->{'rows'} ? 1 : 0 : 0;

      $availability->{'history'}              = !!$rows;
      $availability->{'gene'}                 = 1;
      $availability->{'core'}                 = $self->get_db eq 'core';
      $availability->{'has_gene_tree'}        = $member ? $member->has_GeneTree : 0;
      $availability->{'can_r2r'}              = $self->hub->species_defs->R2R_BIN;
      if ($availability->{'can_r2r'}) {
        my $tree = $availability->{'has_gene_tree'} ? $self->database('compara')->get_GeneTreeAdaptor->fetch_default_for_Member($member) : undef;
        $availability->{'has_2ndary_cons'}    = $tree && $tree->get_tagvalue('ss_cons') ? 1 : 0;
        $availability->{'has_2ndary'}         = ($availability->{'has_2ndary_cons'} || ($obj->canonical_transcript && scalar(@{$obj->canonical_transcript->get_all_Attributes('ncRNA')}))) ? 1 : 0;
      }
      $availability->{'has_gxa'}              = $self->gxa_check;

      $availability->{'alt_allele'}           = $self->table_info($self->get_db, 'alt_allele')->{'rows'};
      $availability->{'regulation'}           = !!$funcgen_res; 
      $availability->{'has_species_tree'}     = $member ? $member->has_GeneGainLossTree : 0;
      $availability->{'family'}               = !!$counts->{families};
      $availability->{'family_count'}         = $counts->{families};
      $availability->{'not_rnaseq'}           = $self->get_db eq 'rnaseq' ? 0 : 1;
      $availability->{"has_$_"}               = $counts->{$_} for qw(transcripts alignments paralogs orthologs similarity_matches operons structural_variation pairwise_alignments);
      $availability->{'multiple_transcripts'} = $counts->{'transcripts'} > 1;
      $availability->{'not_patch'}            = $obj->stable_id =~ /^ASMPATCH/ ? 0 : 1; ## TODO - hack - may need rewriting for subsequent releases
      $availability->{'has_alt_alleles'} =  scalar @{$self->get_alt_alleles};
      
      if ($self->database('variation')) {
        $availability->{'has_phenotypes'} = $self->get_phenotype;
      }

      if ($self->database('compara_pan_ensembl')) {
        $availability->{'family_pan_ensembl'} = !!$counts->{families_pan};
        $availability->{'has_gene_tree_pan'}  = !!($pan_member && $pan_member->has_GeneTree);
        $availability->{"has_$_"}             = $counts->{$_} for qw(alignments_pan paralogs_pan orthologs_pan);
      }
    } elsif ($obj->isa('Bio::EnsEMBL::Compara::Family')) {
      $availability->{'family'} = 1;
    }
    $self->{'_availability'} = $availability;
  }

  return $self->{'_availability'};
}

sub counts {
  my ($self, $member, $pan_member) = @_;
  my $obj = $self->Obj;

  return {} unless $obj->isa('Bio::EnsEMBL::Gene');
  
  my $key = sprintf '::COUNTS::GENE::%s::%s::%s::', $self->species, $self->hub->core_param('db'), $self->hub->core_param('g');
  my $counts = $self->{'_counts'};
  $counts ||= $MEMD->get($key) if $MEMD;
  
  if (!$counts) {
    $counts = {
      transcripts        => scalar @{$obj->get_all_Transcripts},
      exons              => scalar @{$obj->get_all_Exons},
#      similarity_matches => $self->count_xrefs
      similarity_matches => $self->get_xref_available,
      operons => 0,
      alternative_alleles =>  scalar @{$self->get_alt_alleles},
    };
    if ($obj->feature_Slice->can('get_all_Operons')){
      $counts->{'operons'} = scalar @{$obj->feature_Slice->get_all_Operons};
    }
    $counts->{structural_variation} = 0;

    if ($self->database('variation')){ 
      my $vdb = $self->species_defs->get_config($self->species,'databases')->{'DATABASE_VARIATION'};
      $counts->{structural_variation} = $vdb->{'tables'}{'structural_variation'}{'rows'};
      $counts->{phenotypes} = $self->get_phenotype;
    }
    if ($member) {
      $counts->{'orthologs'}  = $member->number_of_orthologues;
      $counts->{'paralogs'}   = $member->number_of_paralogues;
      $counts->{'families'}   = $member->number_of_families;
    }
    my $alignments = $self->count_alignments;
    $counts->{'alignments'} = $alignments->{'all'} if $self->get_db eq 'core';
    $counts->{'pairwise_alignments'} = $alignments->{'pairwise'} + $alignments->{'patch'};

    ## Add pan-compara if available 
    if ($pan_member) {
      my $compara_dbh = $self->database('compara_pan_ensembl')->dbc->db_handle;

      $counts->{'orthologs_pan'}  = $pan_member->number_of_orthologues;
      $counts->{'paralogs_pan'}   = $pan_member->number_of_paralogues;
      $counts->{'families_pan'}   = $pan_member->number_of_families;

      $counts->{'alignments_pan'} = $self->count_alignments('DATABASE_COMPARA_PAN_ENSEMBL')->{'all'} if $self->get_db eq 'core';
    }    

    ## Add counts from plugins
    $counts = {%$counts, %{$self->_counts($member, $pan_member)}};

    $MEMD->set($key, $counts, undef, 'COUNTS') if $MEMD;
    $self->{'_counts'} = $counts;
  }
  
  return $counts;
}
sub get_phenotype {
  my $self = shift;
  
  my $phen_count  = 0;
  my $pfa         = Bio::EnsEMBL::Registry->get_adaptor($self->species, 'variation', 'PhenotypeFeature');
  $phen_count     = $pfa->count_all_by_Gene($self->Obj);

  if (!$phen_count) {
    my $hgncs = $self->obj->get_all_DBEntries('hgnc') || [];

    if(scalar @$hgncs && $hgncs->[0]) {
      my $hgnc_name = $hgncs->[0]->display_id;
      $phen_count   = $pfa->_check_gene_by_HGNC($hgnc_name) if $hgnc_name; # this method is super-fast as it uses some direct SQL on a nicely indexed table
    }
  }
  
  return $phen_count;
}
sub get_xref_available{
  my $self=shift;
  my $available = ($self->count_xrefs > 0);
  if(!$available){
    my @my_transcripts= @{$self->Obj->get_all_Transcripts};
    my @db_links;
    for (my $i=0; !$available && ($i< scalar @my_transcripts); $i++) {
      eval { 
        @db_links = @{$my_transcripts[$i]->get_all_DBLinks};
      };
            
      for (my $j=0; !$available && ($j< scalar @db_links); $j++) {
        $available = $available || ($db_links[$j]->type eq 'MISC') || ($db_links[$j]->type eq 'LIT');
      }      
    }
  }
  return $available;
}

sub count_xrefs {
  my $self = shift;
  my $type = $self->get_db;
  my $dbc = $self->database($type)->dbc;

  # xrefs on the gene
  my $xrefs_c = 0;
  my $sql = '
    SELECT x.display_label, edb.db_name, edb.status
      FROM gene g, object_xref ox, xref x, external_db edb
     WHERE g.gene_id = ox.ensembl_id
       AND ox.xref_id = x.xref_id
       AND x.external_db_id = edb.external_db_id
       AND ox.ensembl_object_type = "Gene"
       AND g.gene_id = ?';

  my $sth = $dbc->prepare($sql);
  $sth->execute($self->Obj->dbID);
  while (my ($label,$db_name,$status) = $sth->fetchrow_array) {
    #these filters are taken directly from Component::_sort_similarity_links
    #code duplication needs removing, and some of these may well not be needed any more
    next if ($status eq 'ORTH');                        # remove all orthologs
    next if (lc($db_name) eq 'medline');                # ditch medline entries - redundant as we also have pubmed
    next if ($db_name =~ /^flybase/i && $type =~ /^CG/ ); # Ditch celera genes from FlyBase
    next if ($db_name eq 'Vega_gene');                  # remove internal links to self and transcripts
    next if ($db_name eq 'Vega_transcript');
    next if ($db_name eq 'Vega_translation');
    next if ($db_name eq 'GO');
    next if ($db_name eq 'OTTP') && $label =~ /^\d+$/; #ignore xrefs to vega translation_ids
    next if ($db_name =~ /ENSG|OTTG/);
    $xrefs_c++;
  }
  return $xrefs_c;
}

sub count_gene_supporting_evidence {
  #count all supporting_features and transcript_supporting_features for the gene
  #- not used in the tree but keep the code just in case we change our minds again!
  my $self = shift;
  my $obj = $self->Obj;
  my $o_type = $self->get_db;
  my $evi_count = 0;
  my %c;
  foreach my $trans (@{$obj->get_all_Transcripts()}) {
    foreach my $evi (@{$trans->get_all_supporting_features}) {
      my $hit_name = $evi->hseqname;
      $c{$hit_name}++;
    }
    foreach my $exon (@{$trans->get_all_Exons()}) {
      foreach my $evi (@{$exon->get_all_supporting_features}) {
        my $hit_name = $evi->hseqname;
        $c{$hit_name}++;
      }
    }
  }
  return scalar(keys(%c));
}


sub can_export {
  my $self = shift;
  
  return $self->action =~ /^(Export|Sequence|TranscriptComparison|Compara_Alignments|Compara_Tree|SpeciesTree|Compara_Ortholog|Compara_Paralog|Family)$/ ? 0 : $self->availability->{'gene'};
}

1;
