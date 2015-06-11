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

package EnsEMBL::Web::Role::Bio::Gene::Family;

### Compara-specific data-munging for gene pages

use Role::Tiny;

sub get_all_families {
  my $self = shift;
  my $compara_db = shift || 'compara';

  my $families;
  if (ref($self->gene) =~ /Family/) { ## No gene in URL, so CoreObjects fetches a family instead
    ## Explicitly set db connection, as registry is buggy!
    my $family = $self->gene;
    my $dba = $self->database('core', $self->species);
    my $genome_db = $self->database($compara_db)->get_GenomeDBAdaptor->fetch_by_name_assembly($self->species);
    my $members = $family->get_Member_by_source_GenomeDB('ENSEMBLPEP', $genome_db);
    my $info = {'description' => $family->description};
    my $genes = [];
    my $prots = {};
    foreach my $member (@$members) {
      my $gene = $member->gene_member->get_Gene;
      push @$genes, $gene;
      my $protein = $member->get_Translation;
      if ($prots->{$gene->stable_id}) {
        push @{$prots->{$gene->stable_id}}, $protein;
      }
      else {
        $prots->{$gene->stable_id} = [$protein];
      }
    }
    $info->{'genes'}    = $genes;
    $info->{'proteins'} = $prots;
    $info->{'count'}    = @$genes;
    $families->{$self->param('family')} = {'info' => $info};
  }
  else {
    foreach my $transcript (@{$self->get_all_transcripts}) {
      my $trans_families = $transcript->get_families($compara_db);
      while (my ($id, $info) = each (%$trans_families)) {
        if (exists $families->{$id}) {
          push @{$families->{$id}{'transcripts'}}, $transcript;
        }
        else {
          $families->{$id} = {'info' => $info, 'transcripts' => [$transcript]};
        }
      }
    }
  }
  return $families;
}

sub create_family {
  my ($self, $id, $cmpdb) = @_;
  $cmpdb ||= 'compara';
  my $databases = $self->database($cmpdb) ;
  my $family_adaptor;
  eval{ $family_adaptor = $databases->get_FamilyAdaptor };
  if ($@){ warn($@); return {} }
  return $family_adaptor->fetch_by_stable_id($id);
}

1;
