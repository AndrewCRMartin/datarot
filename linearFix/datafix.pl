#!/usr/bin/perl -s

use strict;

$::m = 0.285;
$::c = -32.477;

UsageDie() if(!defined($::m) || !defined($::c) || defined($::h));


# Read the data
<>;
while(<>)
{
    my @fields = split(/,/);
    my $x = $fields[1];
    my $y = $fields[2];

    $y = ($y-$::c)/$::m;
    print "$x $y\n";
}




sub UsageDie
{
    print <<__EOF;
    Usage;;;;;;;
__EOF
}
