#!/usr/bin/perl

while(<STDIN>) {
    chomp;
    my @arr = split(/ /);
    $vocab{$arr[2]}++ if(@arr > 3);
}

my @vs = sort keys %vocab;
unshift @vs, '<phi>';
unshift @vs, '<eps>';
my $num = 0;
for(@vs) {
    print "$_\t".$num++."\n";
}
