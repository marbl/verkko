#!/usr/bin/env perl

my $c = "0000";

while (<STDIN>) {
    open(F, "> part$c.txt");   $c++;
    print F $_;
    close(F);
}
