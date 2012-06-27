# $Id$

package EnsEMBL::Web::Component::Location::ChangeSpecies;

### Module to replace part of the former SyntenyView, in this case 
### the lefthand menu dropdown of syntenous species

use strict;

use EnsEMBL::Web::Form;

use base qw(EnsEMBL::Web::Component::Location);

sub _init {
  my $self = shift;
  $self->cacheable(0);
  $self->ajaxable(1);
}

sub content {
  my $self             = shift;
  my $hub              = $self->hub;
  my $species_defs     = $hub->species_defs;
  my $url              = $hub->url({ otherspecies => undef }, 1);
  my $image_config     = $hub->get_imageconfig('Vsynteny');
  my $vwidth           = $image_config->image_height;
  my $form             = $self->new_form({ id => 'change_sp', action => $url->[0], method => 'get', class => 'nonstd autocenter labels_right check', style => $vwidth ? "width:${vwidth}px" : undef });
  my %synteny_hash     = $species_defs->multi('DATABASE_COMPARA', 'SYNTENY');
  my %synteny          = %{$synteny_hash{$hub->species} || {}};
  my @sorted_by_common = sort { $a->{'common'} cmp $b->{'common'} } map {{ name => $_, common => $species_defs->get_config($_, 'SPECIES_COMMON_NAME') }} keys %synteny;
  my @values;

  foreach my $next (@sorted_by_common) {
    next if $next->{'name'} eq $hub->species;
    push @values, { name => $next->{'common'}, value => $next->{'name'} };
  }

  $form->add_hidden({ name => $_, value => $url->[1]->{$_} }) for keys %{$url->[1]};
  $form->add_element(
    type         => 'DropDownAndSubmit',
    select       => 'select',
    style        => 'narrow',
    name         => 'otherspecies',
    label        => 'Change Species',
    values       => \@values,
    value        => $hub->param('otherspecies') || $hub->param('species') || $self->default_otherspecies,
    button_value => 'Go'
  );

  return $form->render;
}


1;
