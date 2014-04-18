#!/bin/bash
# needed : 
# -> Makefile.am
# -> configure.ac

# 0) cleaning
echo "Nettoyage ..."
rm -rf *.o .deps configure config.cache config.sub config.guess config.log config.status aclocal.m4 Makefile Makefile.in depcomp autom4te.cache sortmedna-*.tar.gz  sortmedna-*.zip sortmedna  indexdb include/config.h include/config.h.in install-sh missing COPYING INSTALL

find . -name "*~" -exec rm \{\} \;
find . -name "stamp-h1" -exec rm \{\} \;

#echo "Appuyer la touche <Entrée> pour Autotools"
#read touche
#case $touche in
#*)	echo "Reprise du script..."
#	;;
#esac

# 1) autotools

aclocal
autoconf
autoheader
automake --add-missing

# 2) create distrib
sh _create_missing.sh
./configure
make clean
make 
make dist
