#!/usr/bin/perl

binmode STDIN, ":utf8";
binmode STDOUT, ":utf8";
while(<STDIN>) {
    chomp;
    s/ //g;
    print join(' ',split(//))."\n";
}
