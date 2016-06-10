#!/usr/bin/perl

use strict;
use List::Util qw(max min);

# return a minimum edit distance path with the following notation, plus cost
#  d=delete i=insert s=substitute e=equal
sub levenshtein {
    my ($s, $t) = @_;
    my (@sc, $m, @tc, $n, %d, @str, $i, $j, $id, $cost, $type, $aid, $bid, $cid, $c);
    my (@d, %str);
    @sc = split(/ /, $s);
    $m = @sc;
    @tc = split(/ /, $t);
    $n = @tc;
    # initialize
    @d = (); @str = ();
    foreach $i (0 .. $m) { $id = pack('S2', $i, 0); $d{$id} = $i; $str{$id} = 'd'x$i; }
    foreach $j (1 .. $n) { $id = pack('S2', 0, $j); $d{$id} = $j; $str{$id} = 'i'x$j; }
    
    foreach $i (1 .. $m) {
        foreach $j (1 .. $n) {
            if($sc[$i-1] eq $tc[$j-1]) {
                $cost = 0; $type = 'e'; # equal
            } else {
                $cost = 1.1; $type = 's'; # substitution
            }

            $aid = pack('S2', $i-1, $j); $a = $d{$aid} + 1; # deletion
            $bid = pack('S2', $i, $j-1); $b = $d{$bid} + 1; # insertion
            $cid = pack('S2', $i-1, $j-1); $c = $d{$cid} + $cost; # insertion
            
            $id = pack('S2', $i, $j);
            
            # we want matches to come at the end, so do deletions/insertions first
            if($a <= $b and $a <= $c) {
                $d{$id} = $a;
                $type = 'd';
                $str{$id} = $str{$aid}.'d';
            }
            elsif($b <= $c) {
                $d{$id} = $b;
                $type = 'i';
                $str{$id} = $str{$bid}.'i';
            }
            else {
                $d{$id} = $c;
                $str{$id} = $str{$cid}.$type;
            }

            delete $d{$cid};
            delete $str{$cid};
            # print "".$sc[$i-1]." ".$tc[$j-1]." $i $j $a $b $c $d[$id] $type\n"
        }
    }

    $id = pack('S2', $m, $n);
    return ($str{$id}, $d{$id}); 
}

sub pad {
    my ($s, $l, $m) = @_;
    return $s . (' ' x ($l-length($s)*$m));
}

use strict;
binmode STDOUT, ":utf8";
open REF, "<:utf8", $ARGV[0];
open TEST, "<:utf8", $ARGV[1];

my ($reflen, $testlen);
my %scores = ();
my($ref, $test);
while($ref = <REF> and $test = <TEST>) {
    chomp $ref;
    chomp $test;
    my ($hist, $score) = levenshtein($ref, $test);

    my @ra = split(/ /, $ref);
    $reflen += @ra;
    my @ta = split(/ /, $test);
    $testlen += @ta;
    my @ha = split(//, $hist);
    my ($rd, $td, $hd, $h, $r, $t, $l);
    while(@ha) {
        $h = shift(@ha);
        $scores{$h}++;
        if($h eq 'e' or $h eq 's') {
            $r = shift(@ra);
            $t = shift(@ta);
        } elsif ($h eq 'i') {
            $r = '';
            $t = shift(@ta);
        } elsif ($h eq 'd') {
            $r = shift(@ra);
            $t = '';
        } else { die "bad history value $h"; }
        # find the length
        $l = max(length($r), length($t))*2 + 1;
        $rd .= pad($r, $l, 2);
        $td .= pad($t, $l, 2);
        $hd .= pad($h, $l, 1);
    }
    print "$rd\n$td\n$hd\n\n";
}

my $total = 0;
for (values %scores) { $total += $_; }
foreach my $k (keys %scores) {
    print "$k: $scores{$k} (".$scores{$k}/$total*100 . "%)\n";
}
my $wer = ($scores{'s'}+$scores{'i'}+$scores{'d'})/$reflen*100;
my $prec = $scores{'e'}/$testlen*100;
my $rec = $scores{'e'}/$reflen*100;
my $fmeas = (2*$prec*$rec)/($prec+$rec);
printf ("WER: %.4f%%\nPrec: %.4f%%\nRec: %.4f%%\nF-meas: %.4f%%\n", $wer, $prec, $rec, $fmeas);

