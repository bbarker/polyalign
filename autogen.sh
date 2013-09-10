#!/bin/sh
libtoolize
aclocal -I m4
autoconf
automake -a
./configure 
make
