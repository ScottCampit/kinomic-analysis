#!/usr/bin/perl -w
use strict;
use warnings;
use SOAP::Lite;
use Digest::SHA qw(sha256_hex);

my $password = sha256_hex("SctCmpt93");
my $parameters = "scampit\@umich.edu,".$password.",organism*Homo sapiens#";

my $resultString_kcatkm = SOAP::Lite
-> uri('http://www.brenda-enzymes.org/soap')
-> proxy('http://www.brenda-enzymes.org/soap/brenda_server.php')
-> getKcatKmValue($parameters)
-> result;

print $resultString_kcatkm, "\n";
