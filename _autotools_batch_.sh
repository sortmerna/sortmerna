#!/bin/bash
# needed : 
# -> Makefile.am
# -> configure.ac

# 0) cleaning
echo "Cleaning ..."
rm -rf *.o .deps configure config.cache config.sub config.guess config.log config.status aclocal.m4 Makefile Makefile.in depcomp autom4te.cache sortmedna-*.tar.gz  sortmedna-*.zip sortmedna  indexdb include/config.h include/config.h.in install-sh missing COPYING INSTALL

find . -name "*~" -exec rm \{\} \;
find . -name "stamp-h1" -exec rm \{\} \;

# 1) autotools

aclocal
autoconf -I m4
autoheader
automake --add-missing

# 2) create distrib
./configure
make clean
make 
make dist
