AC_PREREQ([2.50])
AC_INIT([sortmerna],[2.2],[jenya DT kopylov AT gmail DT fr],[sortmerna],[http://bioinfo.lifl.fr/sortmerna])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([include/config.h])
m4_include([m4/ax_check_zlib.m4])
AM_INIT_AUTOMAKE([subdir-objects])
# AM_SILENT_RULES([yes])
if test -z $CXXFLAGS; then
    CXXFLAGS='-O3'
fi
if test -z $CFLAGS; then
    CFLAGS='-O3'
fi
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O
AC_OPENMP
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_FUNC_MMAP

CHECK_ZLIB

AC_SUBST(SSE2_EXTRA_CFLAGS)
CXXFLAGS="$CXXFLAGS -pipe -funroll-loops -Wall -msse2 -D_FILE_OFFSET_BITS=64 $OPENMP_CFLAGS"
CFLAGS="$CFLAGS -pipe -funroll-loops -msse2 -Wall $OPENMP_CFLAGS"
LIBS="$LIBS -lm"
AC_SUBST(CXXFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(LIBS)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
