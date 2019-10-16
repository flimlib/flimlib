#!/usr/bin/perl

###
# #%L
# FLIMLib package for exponential curve fitting of fluorescence lifetime data.
# %%
# Copyright (C) 2010 - 2015 University of Oxford and Board of Regents of the
# University of Wisconsin-Madison.
# %%
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public
# License along with this program.  If not, see
# <http://www.gnu.org/licenses/gpl-3.0.html>.
# #L%
###

use strict;
use warnings;

$^I = '';       # write into the same file

while (<>) {
   s/-ansi//g; # do the replacement
   print;       # print to the modified file
}
