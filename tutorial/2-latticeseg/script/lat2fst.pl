#!/usr/bin/perl

if(not @ARGV) {
    print STDERR "You must specify a scaling factor for the acoustic model (example lat2fst.pl 0.1)\n";
    exit;
}

my @state;
while(<STDIN>) {
    chomp;
    if(/I=(\S*)\s*t=.*W=(\S*)/) {
        $state[$1] = $2;
    } elsif(/S=(\S*)\s*E=(\S*)\s*a=(\S*)/) {
        print "$1 $2 $state[$2] $state[$2] ". ($3*-1*$ARGV[0]) ."\n";
    }
}
print "$#state 0\n";
