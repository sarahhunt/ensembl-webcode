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

package EnsEMBL::Web::Role::Gene;

## Role to provide munging of data from a Bio::EnsEMBL::Gene object

use Role::Tiny;

######### PAGE-RELATED METHODS ################

sub can_export {
  my $self = shift;

  return $self->action =~ /^(Export|Sequence|TranscriptComparison|Compara_Alignments|Compara_Tree|SpeciesTree|Compara_Ortholog|Compara_Paralog|Family)$/ ? 0 : $self->availability->{'gene'};
}

sub availability {
  my $self = shift;
  my ($database_synonym) = @_;

  if (!$self->{'_availability'}) {
    my $availability = $self->_availability;
    my $obj = $self->api_object;

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
  my $obj = $self->api_object;

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
  $phen_count     = $pfa->count_all_by_Gene($self->api_object);

  if (!$phen_count) {
    my $hgncs = $self->api_object->get_all_DBEntries('hgnc') || [];

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
    my @my_transcripts= @{$self->api_object->get_all_Transcripts};
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
  $sth->execute($self->api_object->dbID);
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
  my $obj = $self->api_object;
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

########## GENERAL GENE-MUNGING METHODS ##############

sub date_format {
## FIXME - a generic formatting method shouldn't be in an Object! 
  my( $self, $time, $format ) = @_;
  my( $d,$m,$y) = (localtime($time))[3,4,5];
  my %S = ('d'=>sprintf('%02d',$d),'m'=>sprintf('%02d',$m+1),'y'=>$y+1900);
  (my $res = $format ) =~s/%(\w)/$S{$1}/ge;
  return $res;
}

sub get_Slice {
  my ($self, $context, $ori) = @_;

  my $slice = $self->Obj->feature_Slice;
  $context  = $slice->length * $1 / 100 if $context =~ /(\d+)%/;
  $slice    = $slice->invert if $ori && $slice->strand != $ori;

  return $slice->expand($context, $context);
}

sub get_all_transcripts {
  my $self = shift;
  unless ($self->{'data'}{'_transcripts'}){
    foreach my $transcript (@{$self->gene()->get_all_Transcripts}){
      my $transcriptObj = $self->new_object(
        'Transcript', $transcript, $self->__data
      );
      $transcriptObj->gene($self->gene);
      push @{$self->{'data'}{'_transcripts'}} , $transcriptObj;
    }
  }
  return $self->{'data'}{'_transcripts'};
}

sub display_xref {
  my $self = shift;
  return undef if $self->Obj->isa('Bio::EnsEMBL::Compara::Family');
  return undef if $self->Obj->isa('Bio::EnsEMBL::ArchiveStableId');
  my $trans_xref = $self->Obj->display_xref();
  return undef unless  $trans_xref;
  (my $db_display_name = $trans_xref->db_display_name) =~ s/(.*HGNC).*/$1 Symbol/; #hack for HGNC name labelling, remove in e58
  return ($trans_xref->display_id, $trans_xref->dbname, $trans_xref->primary_id, $db_display_name, $trans_xref->info_text );
}

sub gene_type {
  my $self = shift;
  my $db = $self->get_db;
  my $type = '';
  if( $db eq 'core' ){
    $type = ucfirst(lc($self->Obj->status))." ".$self->Obj->biotype;
    $type =~ s/_/ /;
    $type ||= $self->db_type;
  } elsif ($db =~ /vega/) {
    my $biotype = ($self->Obj->biotype eq 'tec') ? uc($self->Obj->biotype) : ucfirst(lc($self->Obj->biotype));
    $type = ucfirst(lc($self->Obj->status))." $biotype";
    $type =~ s/_/ /g;
    $type =~ s/unknown //i;
    return $type;
  } else {
    $type = $self->logic_name;
    if ($type =~/^(proj|assembly_patch)/ ){
      $type = ucfirst(lc($self->Obj->status))." ".$self->Obj->biotype;
    }
    $type =~ s/_/ /g;
    $type =~ s/^ccds/CCDS/;
  }
  $type ||= $db;
  if( $type !~ /[A-Z]/ ){ $type = ucfirst($type) } #All lc, so format
  return $type;
}

sub mod_date {
  my $self = shift;
  my $time = $self->gene()->modified_date;
  return $self->date_format( $time,'%d/%m/%y' );
}

sub created_date {
  my $self = shift;
  my $time = $self->gene()->created_date;
  return $self->date_format( $time,'%d/%m/%y' );
}

sub get_author_name {
  my $self = shift;
  my $attribs = $self->Obj->get_all_Attributes('author');
  if (@$attribs) {
    return $attribs->[0]->value;
  } else {
    return undef;
  }
}

sub retrieve_remarks {
  my $self = shift;
  my @remarks = grep {$_ ne 'EnsEMBL merge exception'} map { $_->value } @{ $self->Obj->get_all_Attributes('remark') };
  return \@remarks;
}

sub vega_projection {
  my $self = shift;
  my $alt_assembly = shift;
  my $alt_projection = $self->Obj->feature_Slice->project('chromosome', $alt_assembly);
  my @alt_slices = ();
  foreach my $seg (@{ $alt_projection }) {
    my $alt_slice = $seg->to_Slice;
    push @alt_slices, $alt_slice;
  }
  return \@alt_slices;
}

sub get_similarity_hash {
  my ($self, $recurse, $obj) = @_;
  $obj ||= $self->Obj;
  my $DBLINKS;
  eval { $DBLINKS = $obj->get_all_DBEntries; };
  warn ("SIMILARITY_MATCHES Error on retrieving gene DB links $@") if ($@);
  return $DBLINKS  || [];
}

sub get_database_matches {
  my $self = shift;
  my $dbpat = shift;
  my @DBLINKS;
  eval { @DBLINKS = @{$self->Obj->get_all_DBLinks($dbpat)};};
  return \@DBLINKS  || [];
}

sub get_alternative_locations {
  my $self = shift;
  my @alt_locs = map { [ $_->slice->seq_region_name, $_->start, $_->end, $_->slice->coord_system->name ] }
     @{$self->Obj->get_all_alt_locations};
  return \@alt_locs;
}

sub get_rnaseq_tracks {
  my $self = shift;
  my $tracks = [];
  my $rnaseq_db = $self->hub->database('rnaseq');
  if ($rnaseq_db) {
    my $aa = $self->hub->get_adaptor('get_AnalysisAdaptor', 'rnaseq');
    $tracks = [ grep { $_->displayable } @{$aa->fetch_all} ];
  }
  return $tracks;
}


sub insdc_accession {
  my $self = shift;

  my $csv = $self->Obj->slice->coord_system->version;
  my $csa = Bio::EnsEMBL::Registry->get_adaptor($self->species,'core',
                                                'CoordSystem');
  # 0 = look on chromosome
  # 1 = look on supercontig/scaffold
  # maybe in future 2 = ... ?
  for(my $method = 0;$method < 2;$method++) {
    my $slice;
    if($method == 0) {
      $slice = $self->Obj->slice->sub_Slice($self->Obj->start,
                                            $self->Obj->end);
    } elsif($method == 1) {
      # Try to project to supercontig (aka scaffold)
      foreach my $level (qw(supercontig scaffold)) {
        next unless $csa->fetch_by_name($level,$csv);
        my $gsa = $self->Obj->project($level,$csv);
        if(@$gsa==1) {
          $slice = $gsa->[0]->to_Slice;
          last;
        }
      }
    }
    if($slice) {
      my $name = $self->_insdc_synonym($slice,'INSDC');
      if($name) {
        return join(':',$slice->coord_system->name,$csv,$name,
                      $slice->start,$slice->end,$slice->strand);
      }
    }
  }
  return undef;
}

sub _insdc_synonym {
  my ($self,$slice,$name) = @_;

  my $dbc = $self->database($self->get_db)->dbc;
  my $sql = qq(
    SELECT external_db_id FROM external_db WHERE db_name = ?
  );
  my $sth = $dbc->prepare($sql);
  $sth->execute($name);
  my ($dbid) = $sth->fetchrow_array;
  foreach my $s (@{$slice->get_all_synonyms()}) {
    return $s->name if $s->external_db_id == $dbid;
  }
  return undef;
}

sub get_gene_supporting_evidence {
  #get supporting evidence for the gene: transcript_supporting_features support the
  #whole transcript or the translation, supporting_features provide depth the the evidence
  my $self    = shift;
  my $obj     = $self->Obj;
  my $species = $self->species;
  my $ln      = $self->logic_name;
  my $dbentry_adap = Bio::EnsEMBL::Registry->get_adaptor($species, "core", "DBEntry");
  my $o_type  = $self->get_db;
  my $e;
  foreach my $trans (@{$obj->get_all_Transcripts()}) {
    my $tsi = $trans->stable_id;
    my %t_hits;
    my %vega_evi;
  EVI:
    foreach my $evi (@{$trans->get_all_supporting_features}) {
      my $name = $evi->hseqname;
      my $db_name = $dbentry_adap->get_db_name_from_external_db_id($evi->external_db_id);
      #save details of evidence for vega genes for later since we need to combine them 
      #before we can tell if they match the CDS / UTR
      if ($ln =~ /otter/) {
        push @{$vega_evi{$name}{'data'}}, $evi;
        $vega_evi{$name}->{'db_name'} = $db_name;
        $vega_evi{$name}->{'evi_type'} = ref($evi);
        next EVI;
      }

      #for e! genes...
      #use coordinates to check if the transcript evidence supports the CDS, UTR, or just the transcript
      #for protein features give some leeway in matching to transcript - +- 3 bases
      if ($evi->isa('Bio::EnsEMBL::DnaPepAlignFeature')) {
        if ((abs($trans->coding_region_start-$evi->seq_region_start) < 4)
                 || (abs($trans->coding_region_end-$evi->seq_region_end) < 4)) {
          $e->{$tsi}{'evidence'}{'CDS'}{$name} = $db_name;
          $t_hits{$name}++;
        }
        else {
          $e->{$tsi}{'evidence'}{'UNKNOWN'}{$name} = $db_name;
          $t_hits{$name}++;
        }
      }
      elsif ( $trans->coding_region_start == $evi->seq_region_start
                || $trans->coding_region_end == $evi->seq_region_end ) {
        $e->{$tsi}{'evidence'}{'CDS'}{$name} = $db_name;
        $t_hits{$name}++;
      }

      elsif ( $trans->seq_region_start  == $evi->seq_region_start
                || $trans->seq_region_end == $evi->seq_region_end ) {
        $e->{$tsi}{'evidence'}{'UTR'}{$name} = $db_name;
        $t_hits{$name}++;
      }
      else {
        $e->{$tsi}{'evidence'}{'UNKNOWN'}{$name} = $db_name;
        $t_hits{$name}++;
      }
    }
    $e->{$tsi}{'logic_name'} = $trans->analysis->logic_name;

    foreach my $isf (@{$trans->get_all_IntronSupportingEvidence}) {
      push @{$e->{$tsi}{'intron_supporting_evidence'}},$isf->hit_name;
    }

 #make a note of the hit_names of the supporting_features (but don't bother for vega db genes)
    if ($ln !~ /otter/) {
      foreach my $exon (@{$trans->get_all_Exons()}) {
        foreach my $evi (@{$exon->get_all_supporting_features}) {
          my $hit_name = $evi->hseqname;
          if (! exists($t_hits{$hit_name})) {
            $e->{$tsi}{'extra_evidence'}{$hit_name}++;
          }
        }
      }
    }

    #look at vega evidence to see if it can be assigned to 'CDS' 'UTR' etc
    while ( my ($hit_name,$rec) = each %vega_evi ) {
      my ($min_start,$max_end) = (1e8,1);
      my $db_name  = $rec->{'db_name'};
      my $evi_type = $rec->{'evi_type'};
      foreach my $hit (@{$rec->{'data'}}) {
        $min_start = $hit->seq_region_start <= $min_start ? $hit->seq_region_start : $min_start;
        $max_end   = $hit->seq_region_end   >= $max_end   ? $hit->seq_region_end   : $max_end;
      }
      if ($evi_type eq 'Bio::EnsEMBL::DnaPepAlignFeature') {
        #protein evidence supports CDS
        $e->{$tsi}{'evidence'}{'CDS'}{$hit_name} = $db_name;
      }
      else {
        if ($min_start < $trans->coding_region_start && $max_end > $trans->coding_region_end) {
          #full length DNA evidence supports CDS
          $e->{$tsi}{'evidence'}{'CDS'}{$hit_name} = $db_name;
        }
        if (  $max_end   < $trans->coding_region_start
           || $min_start > $trans->coding_region_end
           || $trans->seq_region_start  == $min_start
           || $trans->seq_region_end    == $max_end ) {
          #full length DNA evidence or that exclusively in the UTR supports the UTR
          $e->{$tsi}{'evidence'}{'UTR'}{$hit_name} = $db_name;
        }
        elsif (! $e->{$tsi}{'evidence'}{'CDS'}{$hit_name}) {
          $e->{$tsi}{'evidence'}{'UNKNOWN'}{$hit_name} = $db_name;
        }
      }
    }  
  }
  return $e;
}

# generate URLs for evidence links
sub add_evidence_links {
  my $self = shift;
  my $ids  = shift;
  my $links = [];
  foreach my $hit_name (sort keys %$ids) {
    my $db_name = $ids->{$hit_name};
    my $display = $self->hub->get_ExtURL_link( $hit_name, $db_name, $hit_name );
    push @{$links}, [$display,$hit_name];
  }
  return $links;
}

sub get_gene_slices {
  my ($self, $master_config, @slice_configs) = @_;
  foreach my $array (@slice_configs) {
    if ($array->[1] eq 'normal') {
      my $slice = $self->get_Slice($array->[2], 1);
      $self->__data->{'slices'}{$array->[0]} = [ 'normal', $slice, [], $slice->length ];
    } else {
      $self->__data->{'slices'}{$array->[0]} = $self->get_munged_slice($master_config, $array->[2], 
1);
    }
  }
}

sub get_munged_slice {
  my $self          = shift;
  my $master_config = undef;
  if (ref($_[0]) =~ /ImageConfig/) {
    $master_config = shift;
  } else {
    shift;
  }

  my $slice         = $self->get_Slice(@_);
  my $stable_id     = $self->stable_id;
  my $length        = $slice->length;
  my $munged        = '0' x $length;
  my $context       = $self->param('context') || 100;
  my $extent        = $context eq 'FULL' ? 5000 : $context;
  my $features      = $slice->get_all_Genes(undef, $self->param('opt_db'));
  my @lengths;
  my $offset;

  # Allow return of data for a single transcript
  my $page_type     = $self->hub->type;
  my $focus_transcript = $page_type eq 'Transcript' ? $self->param('t') : undef;

  if ($context eq 'FULL') {
    @lengths = ($length);
  } else {
    foreach my $gene (grep { $_->stable_id eq $stable_id } @$features) {
      foreach my $transcript (@{$gene->get_all_Transcripts}) {
        next if $focus_transcript && $transcript->stable_id ne $focus_transcript;
        if (defined($offset)) {
          if ($offset > $transcript->start-$extent) {
            $offset = $transcript->start-$extent;
          }
        } else {
          $offset = $transcript->start-$extent;
        }
        foreach my $exon (@{$transcript->get_all_Exons}) {
          my $start       = $exon->start - $extent;
          my $exon_length = $exon->end   - $exon->start + 1 + 2 * $extent;
          if ($start-1 >= 0) {
            substr($munged, $start - 1, $exon_length) = '1' x $exon_length;
          } else {
            warn "Got negative substr when munging slices - don't think this should happen\n";
            substr($munged, 0, $exon_length - $start) = '1' x $exon_length;
          }
        }
      }
    }
    $munged =~ s/^0+//;
    $munged =~ s/0+$//;
    @lengths = map length($_), split /(0+)/, $munged;
  }

  # @lengths contains the sizes of gaps and exons(+/- context)

  $munged = undef;

  my $collapsed_length = 0;
  my $flag             = 0;
  my $subslices        = [];
  my $pos              = $offset;

  foreach (@lengths, 0) {
    if ($flag = 1 - $flag) {
      push @$subslices, [ $pos + 1, 0, 0 ];
      $collapsed_length += $_;
    } else {
      $subslices->[-1][1] = $pos;
    }

    $pos += $_;
  }
  # compute the width of the slice image within the display
  my $pixel_width =
    ($master_config ? $master_config->get_parameter('image_width') : 800) -
    ($master_config ? $master_config->get_parameter('label_width') : 100) -
    ($master_config ? $master_config->get_parameter('margin')      :   5) * 3;

  # Work out the best size for the gaps between the "exons"
  my $fake_intron_gap_size = 11;
  my $intron_gaps          = $#lengths / 2;

  if ($intron_gaps * $fake_intron_gap_size > $pixel_width * 0.75) {
    $fake_intron_gap_size = int($pixel_width * 0.75 / $intron_gaps);
  }

  # Compute how big this is in base-pairs
  my $exon_pixels    = $pixel_width - $intron_gaps * $fake_intron_gap_size;
  my $scale_factor   = $collapsed_length / $exon_pixels;
  my $padding        = int($scale_factor * $fake_intron_gap_size) + 1;
  $collapsed_length += $padding * $intron_gaps;

  # Compute offset for each subslice
  my $start = 0;
  foreach (@$subslices) {
    $_->[2] = $start  - $_->[0];
    $start += $_->[1] - $_->[0] - 1 + $padding;
  }

  return [ 'munged', $slice, $subslices, $collapsed_length ];
}


# Calls for HistoryView

sub get_archive_object {
  my $self = shift;
  my $id = $self->stable_id;
  my $archive_adaptor = $self->database('core')->get_ArchiveStableIdAdaptor;
  my $archive_object = $archive_adaptor->fetch_by_stable_id($id, 'Gene');
  return $archive_object;
}

=head2 history

 Arg1        : data object
 Description : gets the deduplicated archive id history tree based around this ID
 Return type : listref of Bio::EnsEMBL::ArchiveStableId
               As every ArchiveStableId knows about it's successors, this is
                a linked tree.

=cut

sub history {
  my $self = shift;

  my $archive_adaptor = $self->database('core')->get_ArchiveStableIdAdaptor;
  return unless $archive_adaptor;

  my $history = $archive_adaptor->fetch_history_tree_by_stable_id($self->stable_id);

  return $history;
}

=head2 get_predecessors

 Arg1        : data object
 Description : gets the complete archive id history tree based around this ID
 Return type : listref of Bio::EnsEMBL::ArchiveStableId

=cut

sub get_predecessors {
  my $self = shift;
  my $archive_adaptor = $self->database('core')->get_ArchiveStableIdAdaptor;
  my $archive = $archive_adaptor->fetch_by_stable_id($self->stable_id, 'Gene');
  return [] unless $archive;
  my $predecessors = $archive_adaptor->fetch_predecessor_history($archive);
  return $predecessors;
}

1;
