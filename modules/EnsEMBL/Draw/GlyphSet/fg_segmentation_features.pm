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

package EnsEMBL::Draw::GlyphSet::fg_segmentation_features;

### Draw regulatory segmentation features track (semi-continuous
### track of colour blocks)

use strict;

use base qw(EnsEMBL::Draw::GlyphSet::bigbed);

sub get_data {
  my $self    = shift;
  my $slice   = $self->{'container'};
  my $db_type = $self->my_config('db_type') || 'funcgen';
  my $fg_db;

  if (!$slice->isa('Bio::EnsEMBL::Compara::AlignSlice::Slice')) {
    $fg_db = $slice->adaptor->db->get_db_adaptor($db_type);

    if(!$fg_db) {
      warn "Cannot connect to $db_type db";
      return [];
    }
  }

  my $f = $self->fetch_features_from_file($fg_db);
  return $f if defined $f;
  return $self->fetch_features_from_db($fg_db);
}

sub _feature_set {
  my ($self) = @_;

  my $slice   = $self->{'container'};
  my $db_type = $self->my_config('db_type') || 'funcgen';
  return undef if $slice->isa('Bio::EnsEMBL::Compara::AlignSlice::Slice');
  my $fgh = $slice->adaptor->db->get_db_adaptor($db_type);
  my $fsa       = $fgh->get_FeatureSetAdaptor();
  my $cta       = $fgh->get_CellTypeAdaptor;
  return undef unless $fsa and $cta;
  my $cell_line = $self->my_config('cell_line');
  return undef unless $cell_line;
  my $ctype = $cta->fetch_by_name($cell_line);
  my $fsets = $fsa->fetch_all_displayable_by_type('segmentation', $ctype);
  return undef unless $fsets and @$fsets;
  return $fsets->[0];
}

sub fetch_features_from_db {
  my ($self, $db) = @_;

  my $fset = $self->_feature_set();
  return [] unless $fset;


  $self->{'legend'}{'fg_regulatory_features_legend'} ||= { priority => 1020, legend => [] };

  my $slice     = $self->{'container'};
  my $features = $fset->get_Features_by_Slice($slice);
  my @dff;
  foreach my $f (@$features) {
    my $colour_key = $self->colour_key($f);
    my $colour = $self->my_colour($colour_key) || '#e1e1e1';
    my $text = $self->my_colour($colour_key,'text');
    push @dff,{
      colour => $colour,
      label_colour => $colour,
      join_colour => $colour,
      start => $f->start,
      end => $f->end,
      score => 1000,
      label => $text,
    };
  }

  return [{
    features => { '-1' => \@dff },
    metadata => {
      force_strand => '-1',
      default_strand => 1,
      omit_feature_links => 1,
      display => 'normal'
    }
  }];
}

sub _result_set {
  my ($self) = @_;

  my $slice   = $self->{'container'};
  my $db_type = $self->my_config('db_type') || 'funcgen';
  return undef if $slice->isa('Bio::EnsEMBL::Compara::AlignSlice::Slice');
  my $fgh = $slice->adaptor->db->get_db_adaptor($db_type);
  return undef unless $fgh;
  my $cell_line = $self->my_config('celltype');
  return undef unless $cell_line;
  my $cta = $fgh->get_CellTypeAdaptor;
  my $ct = $cta->fetch_by_dbID($cell_line);
  my $rsa = $fgh->get_ResultSetAdaptor;
  my $rsets = $rsa->fetch_all_by_CellType($ct);
  my @segs = grep { $_->feature_class eq 'segmentation' } @$rsets;
  return undef unless @segs;
  return $segs[0];
}

sub fetch_features_from_file {
  my ($self,$fgh) = @_;

  my $rs = $self->_result_set();
  return undef unless $rs;
  my $bigbed_file = $rs->dbfile_path;
  my $file_path = join('/',$self->species_defs->DATAFILE_BASE_PATH,
                           lc $self->species,
                           $self->species_defs->ASSEMBLY_VERSION);
  $bigbed_file = "$file_path/$bigbed_file" unless $bigbed_file =~ /^$file_path/;
  $bigbed_file =~ s/\s//g;
  my $out = $self->SUPER::get_data($bigbed_file);
  return $out;
}

sub href { return undef; }
sub bg_link {
  my ($self, $strand) = @_;

  my $rs = $self->_result_set();
  my $fs = $self->_feature_set();
  if($rs) {
    return $self->_url({
      action   => 'SegFeature',
      ftype    => 'Regulation',
      dbid     => $rs->dbID,
      species  => $self->species,
      fdb      => 'funcgen',
      scalex   => $self->scalex,
      width    => $self->{'container'}->length,
      celldbid => $self->my_config('celltype'),
    });
  } elsif($fs) {
    return $self->_url({
      action   => 'SegFeature',
      ftype    => 'Regulation',
      dbid     => $fs->dbID,
      species  => $self->species,
      fdb      => 'funcgen',
      scalex   => $self->scalex,
      width    => $self->{'container'}->length,
      cl       => $self->my_config('cell_line'),
    });
  } else {
    return undef;
  }
}

sub colour_key {
  my ($self, $f) = @_;
  my $type = $f->feature_type->name;

  if ($type =~ /Repressed/ or $type =~ /low activity/) {
    $type = 'repressed';
  } elsif ($type =~ /CTCF/) {
    $type = 'ctcf';
  } elsif ($type =~ /Enhancer/) {
    $type = 'enhancer';
  } elsif ($type =~ /Flank/) {
    $type = 'promoter_flanking';
  } elsif ($type =~ /TSS/) {
    $type = 'promoter';
  } elsif ($type =~ /Transcribed/) {
    $type = 'region';
  } elsif ($type =~ /Weak/) {
    $type = 'weak';
  } elsif ($type =~ /Heterochr?omatin/i) { # ? = typo in e76
    $type = 'heterochromatin';
  } else {
    $type = 'default';
  }
  return lc $type;
}

sub render {
  my ($self) = @_;

  $self->{'my_config'}->set('link_on_bgd', 1);
  $self->render_compact;
}

1;
