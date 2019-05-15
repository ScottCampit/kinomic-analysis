#!/usr/bin/perl -w
use strict;
use warnings;
use SOAP::Lite;
use Digest::SHA qw(sha256_hex);

my $password = sha256_hex("SctCmpt93");
my $parameters = "scampit\@umich.edu,".$password.",organism*Mus musculus#";

my $resultString = SOAP::Lite
-> uri('http://brenda-enzymes.org/soap')
-> proxy('http://brenda-enzymes.org/soap/brenda_server.php')
-> getTurnoverNumber($parameters)
-> result;

print $resultString, "\n";
