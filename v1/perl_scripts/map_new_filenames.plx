#!/usr/bin/perl

# Teresita M. Porter, June 18, 2019

# Script to use a mapping file to rename raw read files
# Usage perl map_new_filenames.plx map.txt


use strict;
use warnings;

# declare var

# declare array
my @map;
my @files;

# declare hash

open (MAP, "<", $ARGV[0]) || die "Error can't open mapping file: $!\n";
@map = <MAP>;
close MAP;

print "Enter path to dir containing raw files (.gz) including final '/':\n";
@files = <STDIN>;


