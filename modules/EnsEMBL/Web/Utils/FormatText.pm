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

package EnsEMBL::Web::Utils::FormatText;

## Handy methods for formatting strings, dates, etc

use base qw(Exporter);

our @EXPORT = our @EXPORT_OK = qw(date_format pretty_date);

sub date_format {
### Generic method for formatting a unix timestamp, according to a simple format
### E.g. given a format '%d/%m/%y', return '31/12/99'
### @param time - Unix timestamp
### @param format (optional) - String containing POSIX datetime variables and
###                            optional punctuation characters
### @return String
  my ($time, $format) = @_;
  $format ||= '%d/%m/%y';
  my ($d,$m,$y) = (localtime($time))[3,4,5];
  my %S = ('d' => sprintf('%02d',$d), 'm' => sprintf('%02d',$m+1), 'y' => $y+1900);
  (my $res = $format) =~ s/%(\w)/$S{$1}/ge;
  return $res;
}

sub pretty_date {
### Format a unix timestamp as e.g. 'Mon 1 Jan, 2015'
### @param timestamp - Unix timestamp
### @return String
  my $timestamp = shift;
  my @date = localtime($timestamp);
  my @days = ('Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat');
  my @months = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec');
  return $days[$date[6]].' '.$date[3].' '.$months[$date[4]].', '.($date[5] + 1900);
}

1;
