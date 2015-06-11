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

package EnsEMBL::Web::Role::Bio::Gene::Tree;

### Compara-specific data-munging for gene pages

use Role::Tiny;

sub get_GeneTree {
  my $self       = shift;
  my $compara_db = shift || 'compara';
  my $whole_tree = shift;
  my $clusterset_id = $self->hub->param('clusterset_id') || 'default';
  my $cache_key  = sprintf('_protein_tree_%s_%s', $compara_db, $clusterset_id);

  if (!$self->{$cache_key}) {
    my $member  = $self->get_compara_Member($compara_db)           || return;
    my $adaptor = $member->adaptor->db->get_adaptor('GeneTree')    || return;
    my $tree    = $adaptor->fetch_all_by_Member($member, -clusterset_id => $clusterset_id)->[0];
    unless ($tree) {
        $tree = $adaptor->fetch_default_for_Member($member);
    }
    return unless $tree;
    return $tree if $whole_tree;

    $tree->preload;
    $self->{$cache_key} = $tree->root;
    $self->{"_member_$compara_db"} = $member;

    my $parent      = $adaptor->fetch_parent_tree($tree);
    if ($parent->tree_type ne 'clusterset') {
      my %subtrees;
      my $total_leaves = 0;
      foreach my $subtree (@{$adaptor->fetch_subtrees($parent)}) {
        $subtrees{$subtree->{_parent_id}} = ($tree->root_id eq $subtree->root_id ? $tree : $subtree);
      }
      $parent->preload;
      foreach my $leaf (@{$parent->root->get_all_leaves}) {
        my $subtree = $subtrees{$leaf->node_id};
        $leaf->{'_subtree'} = $subtree;
        $leaf->{'_subtree_size'} = $subtree->get_tagvalue('gene_count');
        $total_leaves += $leaf->{'_subtree_size'};
      }
      $parent->{'_total_num_leaves'} = $total_leaves;
      $tree->{'_supertree'} = $parent;
    }
  }
  return $self->{$cache_key};
}

sub get_SpeciesTree {
  my $self       = shift;
  my $compara_db = shift || 'compara';

  my $hub            = $self->hub;
  my $collapsability = $hub->param('collapsability');
  my $cache_key      = "_species_tree_".$collapsability."_".$compara_db;
  my $database       = $self->database($compara_db);

  if (!$self->{$cache_key}) {
    my $cafeTree_Adaptor = $database->get_CAFEGeneFamilyAdaptor();
    my $geneTree_Adaptor = $database->get_GeneTreeAdaptor();

    my $member   = $self->get_compara_Member($compara_db)           || return;
    my $geneTree = $geneTree_Adaptor->fetch_default_for_Member($member) || return;
    my $cafeTree = $cafeTree_Adaptor->fetch_by_GeneTree($geneTree) || return;

    $cafeTree->multifurcate_tree();
    $cafeTree    = $cafeTree->root($cafeTree->root->lca_reroot($cafeTree->lca_id)) if($collapsability eq 'part');

    $self->{$cache_key} = $cafeTree;
  }

  return $self->{$cache_key};
}

1;
