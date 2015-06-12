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

package EnsEMBL::Web::Role::Bio::Gene::Variation;

### Variation-specific data-munging for gene pages

use Role::Tiny;

# Calls for GeneSNPView

# Valid user selections
sub valids {
  my $self = shift;
  my %valids = (); # Now we have to create the snp filter

  foreach ($self->param) {
    $valids{$_} = 1 if $_ =~ /opt_/ && $self->param($_) eq 'on';
  }

  return \%valids;
}

sub getVariationsOnSlice {
  my( $self, $slice, $subslices, $gene, $so_terms ) = @_;
  my $sliceObj = $self->new_object('Slice', $slice, $self->__data);

  my ($count_snps, $filtered_snps, $context_count) = $sliceObj->getFakeMungedVariationFeatures($subslices,$gene,$so_terms);
  $self->__data->{'sample'}{"snp_counts"} = [$count_snps, scalar @$filtered_snps];
  $self->__data->{'SNPS'} = $filtered_snps;
  return ($count_snps, $filtered_snps, $context_count);
}

sub store_TransformedTranscripts {
  my( $self ) = @_;
  my $focus_transcript = $self->hub->type eq 'Transcript' ? $self->param('t') : undef;

  my $offset = $self->__data->{'slices'}{'transcripts'}->[1]->start -1;
  foreach my $trans_obj ( @{$self->get_all_transcripts} ) {
    next if $focus_transcript && $trans_obj->stable_id ne $focus_transcript;
    my $transcript = $trans_obj->api_object;
    my ($raw_coding_start,$coding_start);
    if (defined( $transcript->coding_region_start )) {
      $raw_coding_start = $transcript->coding_region_start;
      $raw_coding_start -= $offset;
      $coding_start = $raw_coding_start + $self->munge_gaps( 'transcripts', $raw_coding_start );
    }
    else {
      $coding_start  = undef;
      }

    my ($raw_coding_end,$coding_end);
    if (defined( $transcript->coding_region_end )) {
      $raw_coding_end = $transcript->coding_region_end;
      $raw_coding_end -= $offset;
      $coding_end = $raw_coding_end   + $self->munge_gaps( 'transcripts', $raw_coding_end );
    }
    else {
      $coding_end = undef;
    }
    my $raw_start = $transcript->start;
    my $raw_end   = $transcript->end  ;
    my @exons = ();
    foreach my $exon (@{$transcript->get_all_Exons()}) {
      my $es = $exon->start - $offset;
      my $ee = $exon->end   - $offset;
      my $O = $self->munge_gaps( 'transcripts', $es );
      push @exons, [ $es + $O, $ee + $O, $exon ];
    }
    $trans_obj->__data->{'transformed'}{'exons'}        = \@exons;
    $trans_obj->__data->{'transformed'}{'coding_start'} = $coding_start;
    $trans_obj->__data->{'transformed'}{'coding_end'}   = $coding_end;
    $trans_obj->__data->{'transformed'}{'start'}        = $raw_start;
    $trans_obj->__data->{'transformed'}{'end'}          = $raw_end;
  }
}

sub get_included_so_terms {
  my $self     = shift;

  # map the selected consequence type to SO terms
  my %cons = %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;
  my %options =  EnsEMBL::Web::Constants::VARIATION_OPTIONS;

  my $hub  = $self->hub;

  my %selected_so;

  foreach ($hub->param) {
    if ($hub->param($_) eq 'on' && $_ =~ /opt_/ && exists($options{'type'}{$_})) {
      foreach my $con (keys %cons) {
        my $consequence = "opt_" . lc $cons{$con}->SO_term;
        $selected_so{$con} = 1 if $_ eq $consequence;
      }
    }
  }

  my @so_terms = keys %selected_so;
  return \@so_terms;
}

sub store_TransformedSNPS {
  my $self     = shift;
  my $so_terms = shift;
  my $vfs      = shift;

  my $valids   = $self->valids;

  my $tva = $self->get_adaptor('get_TranscriptVariationAdaptor', 'variation');

  my @transcripts = @{$self->get_all_transcripts};
  if ($self->hub->type eq 'Transcript'){
    @transcripts = ($self->hub->core_object('transcript'));
  }

  my $included_so;
  if ($self->need_consequence_check) {
    $included_so = $self->get_included_so_terms;
  }

  # get all TVs and arrange them by transcript stable ID and VF ID, ignore non-valids
  my $tvs_by_tr;

  my $method        = 'fetch_all_by_VariationFeatures';
  my $have_so_terms = (defined $so_terms && scalar @$so_terms);
  my $filtered_vfs  = $vfs;

  if ($have_so_terms) {
    # tva needs an ontology term adaptor to fetch by SO term
    $tva->{_ontology_adaptor} ||= $self->hub->get_databases('go')->{'go'}->get_OntologyTermAdaptor;

    $method .= '_SO_terms';

    my %term_hash;
    foreach my $so_term (@$so_terms) {
      $term_hash{$so_term} = 1;
    }

    my @vfs_with_term = grep { scalar map { $term_hash{$_} ? 1 : () } @{$_->consequence_type('SO')} } @$vfs;
    $filtered_vfs = \@vfs_with_term;
  }

  my $tvs;
  if (!$have_so_terms && $included_so ) {
    $tva->{_ontology_adaptor} ||= $self->hub->get_databases('go')->{'go'}->get_OntologyTermAdaptor;
    $tvs = $tva->fetch_all_by_VariationFeatures_SO_terms($filtered_vfs,[map {$_->transcript} @transcripts],$included_so,1) ;
  } else {
    $tvs = $tva->$method($filtered_vfs,[map {$_->transcript} @transcripts],$so_terms,0, $included_so) ;
  }

  if (!$self->need_consequence_check) {
    foreach my $tv (@$tvs) {
      $tvs_by_tr->{$tv->transcript->stable_id}->{$tv->variation_feature->dbID} = $tv;
    }
  } else {
    foreach my $tv (@$tvs) {
      my $found = 0;
      foreach my $type(@{$tv->consequence_type || []}) {
        if (exists($valids->{'opt_'.lc($type)})) {
          $tvs_by_tr->{$tv->transcript->stable_id}->{$tv->variation_feature->dbID} = $tv;
          $found=1;
          last;
        }
      }
    }
  }

  # then store them in the transcript's data hash
  my $total_tv_count = 0;
  foreach my $trans_obj (@{$self->get_all_transcripts}) {
    $trans_obj->__data->{'transformed'}{'snps'} = $tvs_by_tr->{$trans_obj->stable_id};
    $total_tv_count += scalar(keys %{$tvs_by_tr->{$trans_obj->stable_id}});
  }
}

sub store_ConsequenceCounts {
  my $self     = shift;
  my $so_term_sets = shift;
  my $vfs      = shift;

  my $tva = $self->get_adaptor('get_TranscriptVariationAdaptor', 'variation');

  my @transcripts = @{$self->get_all_transcripts};
  if ($self->hub->type eq 'Transcript'){
    @transcripts = ($self->hub->core_object('transcript'));
  }

  my $included_so;

  if ($self->need_consequence_check) {
    # Can't use counts with consequence check so clear any existing stored conscounts and return - no longer true
    # if (exists($self->__data->{'conscounts'})) { delete $self->__data->{'conscounts'}; }
    #return;
    $included_so = $self->get_included_so_terms;
  }

  $tva->{_ontology_adaptor} ||= $self->hub->get_databases('go')->{'go'}->get_OntologyTermAdaptor;

  my %conscounts;

  foreach my $cons (keys %$so_term_sets) {
    my $filtered_vfs = $vfs;

    my $so_terms = $so_term_sets->{$cons};

    my %term_hash = map {$_ => 1} @$so_terms;

    my @vfs_with_term = grep { scalar map { $term_hash{$_} ? 1 : () } @{$_->consequence_type()} } @$vfs;
    $filtered_vfs = \@vfs_with_term;

    $conscounts{$cons} = $tva->count_all_by_VariationFeatures_SO_terms($filtered_vfs,[map {$_->transcript} @transcripts],$so_terms,$included_so) ;;
  }

  if (!$included_so) {
    $conscounts{'ALL'} = $tva->count_all_by_VariationFeatures($vfs,[map {$_->transcript} @transcripts]) ;
  } else {
    $conscounts{'ALL'} = $tva->count_all_by_VariationFeatures_SO_terms($vfs,[map {$_->transcript} @transcripts], $included_so) ;
  }

  # then store them in the gene's data hash
  $self->__data->{'conscounts'} = \%conscounts;
}

sub need_consequence_check {
  my( $self ) = @_;

  my %options =  EnsEMBL::Web::Constants::VARIATION_OPTIONS;

  my $hub  = $self->hub;

  foreach ($hub->param) {
    if ($hub->param($_) eq 'off' && $_ =~ /opt_/ && exists($options{'type'}{$_})) {
      return 1;
    }
  }

  return 0;
}

sub store_TransformedDomains {
  my( $self, $key ) = @_;
  my %domains;
  my $focus_transcript = $self->hub->type eq 'Transcript' ? $self->param('t') : undef;

  my $offset = $self->__data->{'slices'}{'transcripts'}->[1]->start -1;
  foreach my $trans_obj ( @{$self->get_all_transcripts} ) {
    next if $focus_transcript && $trans_obj->stable_id ne $focus_transcript;
    my %seen;
    my $transcript = $trans_obj->api_object;
    next unless $transcript->translation;
    foreach my $pf ( @{$transcript->translation->get_all_ProteinFeatures( lc($key) )} ) {
## rach entry is an arry containing the actual pfam hit, and mapped start and end co-ordinates
      if (exists $seen{$pf->display_id}{$pf->start}){
        next;
      } else {
        $seen{$pf->display_id}->{$pf->start} =1;
        my @A = ($pf);
        foreach( $transcript->pep2genomic( $pf->start, $pf->end ) ) {
          my $O = $self->munge_gaps( 'transcripts', $_->start - $offset, $_->end - $offset) - $offset;
          push @A, $_->start + $O, $_->end + $O;
        }
        push @{$trans_obj->__data->{'transformed'}{lc($key).'_hits'}}, \@A;
      }
    }
  }
}

sub munge_gaps {
  my( $self, $slice_code, $bp, $bp2  ) = @_;
  my $subslices = $self->__data->{'slices'}{ $slice_code }[2];
  foreach( @$subslices ) {

    if( $bp >= $_->[0] && $bp <= $_->[1] ) {
      return defined($bp2) && ($bp2 < $_->[0] || $bp2 > $_->[1] ) ? undef : $_->[2] ;
    }
  }
  return undef;
}

1;

