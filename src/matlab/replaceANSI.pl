#!/usr/bin/perl

use strict;
use warnings;

$^I = '';       # write into the same file

while (<>) {
   s/-ansi//g; # do the replacement
   print;       # print to the modified file
}