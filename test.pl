use strict;
use SAPserver;

my $sapServer = SAPserver->new(singleton => 1);
my $genomeIDs = $sapServer->all_genomes(-complete => 1);
$genomeIDs = [ keys %$genomeIDs ];
$genomeIDs = ["360095.3"];
for my $genomeID (@$genomeIDs) {
    my $repID = $sapServer->representative(-ids => $genomeID);
    my $taxonomy = $sapServer->taxonomy_of(-ids => $genomeID,
					   -format => 'names');
    # Format the taxonomy string.
    my $taxonomyString = join(" ", @$taxonomy);
    # Print the result.
    print join("\t", $genomeID, $repID, $taxonomyString) . "\n";
}
