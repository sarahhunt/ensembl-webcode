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

package EnsEMBL::Web::Role::Bio;

### Generic role to hold functionality shared by _all_ genomic objects

use Role::Tiny;

sub counts            { return {};        }
sub _counts           { return {};        } # Implemented in plugins
sub availability      { return {};        }
sub can_export        { return 0;         }
sub default_action    { return 'Summary'; }
sub __data            { return $_[0]{'data'};                  }
sub __objecttype      { return $_[0]{'data'}{'_objecttype'};   }

## Old accessor - deprecate in due course
sub Obj               { return $_[0]{'data'}{'_object'};       } # Gets the underlying Ensembl object wrapped by the web object
## New accessor - to keep the webcode cleaner and more readable
sub api_object        { return $_[0]{'data'}{'_object'};       } # Gets the underlying API object
