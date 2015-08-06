=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package EnsEMBL::Web::ImageConfig::sequence_logo;

use strict;

use base qw(EnsEMBL::Web::ImageConfig);

sub init {
  my $self = shift;

  $self->set_parameters({
    sortable_tracks => 1, # allow the user to reorder tracks
    opt_lines       => 1, # draw registry lines
    spritelib       => { default => $self->species_defs->ENSEMBL_WEBROOT . '/htdocs/img/sprites' },
  });

  $self->create_menus(qw(
    sequence
    variation
    other
    information
  ));
  $self->image_resize = 1;

  $self->add_tracks('other',
    [ 'scalebar',  '', 'scalebar',  { display => 'normal', strand => 'b', name => 'Scale bar', description => 'Shows the scalebar' }],
    [ 'ruler',     '', 'ruler',     { display => 'normal', strand => 'b', name => 'Ruler',     description => 'Shows the length of the region being displayed' }],
    [ 'draggable', '', 'draggable', { display => 'normal', strand => 'b', menu => 'no' }],
  );
  
  $self->add_tracks('information',
    [ 'missing', '', 'text', { display => 'normal', strand => 'r', name => 'Disabled track summary', description => 'Show counts of number of tracks turned off by the user' }],
    [ 'info',    '', 'text', { display => 'normal', strand => 'r', name => 'Information',            description => 'Details of the region shown in the image' }]
  );
 
  $self->add_tracks('sequence',
    [ 'seq',       'Sequence',            'sequence', { display => 'normal', strand => 'b', description => 'Track showing reference sequence',       colourset => 'seq',      threshold => 1,   depth => 1      }],
    [ 'logo',       'Logo',            'seq_logo', { display => 'normal', strand => 'b', description => 'Track showing sequence logo', }],
  );
 
  $self->load_tracks;

  $self->modify_configs(
    [ 'fg_regulatory_features_legend', 'fg_segmentation_features_legend', 'fg_methylation_legend' ],
    { display => 'off' }
  );

  $self->modify_configs(
    [ 'variation_feature_variation' ],
    { display => 'normal' }
  );

}

1;
