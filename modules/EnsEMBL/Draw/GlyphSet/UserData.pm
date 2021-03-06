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

package EnsEMBL::Draw::GlyphSet::UserData;

use strict;

use parent qw(EnsEMBL::Draw::GlyphSet::Generic);

sub bg_link {
  my $self = shift;
  return $self->_url({ action => 'UserData' }); 
}

### Renderers

sub render_gradient {
### Features coloured on a gradient by score, e.g. pvalues
  my $self = shift;
  $self->{'my_config'}->set('drawing_style', ['Graph::Heatmap']);
  $self->{'my_config'}->set('height', 8);
  $self->{'my_config'}->set('no_axis', 1);
  $self->{'my_config'}->set('use_pvalue', 1);
  $self->_render_aggregate;
}

1;
