use strict;
use SAPserver;



if (0)
{
my $sapServer = SAPserver->new();
my $rxns = $sapServer->all_reactions();
for my $rxn (@$rxns)
{
    #print "$rxn\n";
}
my $rxnsSmall = [$rxns->[0]];
#$rxnsSmall = ["rxn00001", "rxn00002"];
#print length(@$rxnsSmall) . "\n";
#print $rxnsSmall->[0] . "\n";
}

if (0)
{
my $rxnHash = $sapServer->reaction_strings({-ids => $rxnsSmall});
for my $rxn (sort keys %$rxnHash)
{
    print "$rxn\n";
    my $rolesList = $rxnHash->{$rxn};
    print "$rxn\t" . join(",",@$rolesList) . "\n";
}
}

if (0)
{
    my $genomeIDs = ["360095.3"];
    my $repID = $sapServer->representative(-ids => $genomeIDs);
    my $couples = $sapServer->coupled_reactions({-ids => $rxnsSmall});
}

if (0)
{
my $roles = $sapServer->all_roles_used_in_models();
#my @array = @{$roles};
#print "$array[0]\n";
my $rolesHash = $sapServer->role_reactions({-ids => $roles});
for my $role (@$roles)
{
    print "$role\n";
}
}

if (0)
{
my $roles = $sapServer->all_roles_used_in_models();
my $hash = $sapServer->occ_of_role({-roles => $roles});
for my $id (sort keys %$hash)
{
    my $list = $hash->{$id};
    print "$id\n";
    print join(",",@$list) . "\t";
}
}
