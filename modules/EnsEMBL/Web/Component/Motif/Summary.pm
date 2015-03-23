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

package EnsEMBL::Web::Component::Motif::Summary;

use strict;
use warnings;
no warnings "uninitialized";
use base qw(EnsEMBL::Web::Component);


sub _init {
  my $self = shift;
  $self->cacheable( 0 );
  $self->ajaxable(  0 );
}

sub content {
  my $self    = shift;
  my $object  = $self->object;
  my $html;

  ## Basic information about the motif
  my $twocol  = $self->new_twocol;
  my %summary = %{$object->summary||{}};

  $twocol->add_row('Name', $object->name);

  my $description = $object->binding_matrix->description || '[waiting for API]';
  $twocol->add_row('Description', $description);

  my $r = sprintf('%s:%s-%s', $summary{'seq_region_name'}, $summary{'start'}, $summary{'end'});
  my $url = $self->hub->url({'type' => 'Location', 'action' => 'View', 'r' => $r});
  my $location = sprintf('<a href="%s">%s</a> (%s)', $url, $r, $summary{'strand'});
  $twocol->add_row('Location', $location);

  my $bm_link = $self->hub->get_ExtURL_link($self->object->binding_matrix_name, 'JASPAR', $self->object->binding_matrix_name); 
  $twocol->add_row('Binding matrix', $bm_link);

  $twocol->add_row('Score', $summary{'score'});

  my $cell_lines = '[waiting for API]';
  $twocol->add_row('Evidence', $cell_lines);

  $html .= $twocol->render; 

  ## Regulatory context
  my $rf            = $self->hub->core_object('regulation');
  my $feature_slice = $rf->Obj->feature_Slice;

  my $image_config = $self->hub->get_imageconfig('motif_feature');

  $image_config->set_parameters({
    container_width => $feature_slice->length,
    image_width     => $self->image_width || 800,
    slice_number    => '1|1',
  });

  ## Show the ruler only on the same strand as the motif
  $image_config->modify_configs(
    [ 'ruler', 'scalebar' ],
    { 'strand', $rf->Obj->strand > 0 ? 'f' : 'r' }
  );

  my $image = $self->new_image($feature_slice, $image_config, [ $rf->stable_id ]);
  return if $self->_export_image($image);

  $image->imagemap         = 'yes';
  $image->{'panel_number'} = 'top';
  $image->set_button('drag', 'title' => 'Drag to select region');

  $html .= $image->render;

  return $html;
}

1;
