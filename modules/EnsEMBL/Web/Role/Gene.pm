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

## TODO - factor out compara, variation and regulation methods

use Role::Tiny;

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
