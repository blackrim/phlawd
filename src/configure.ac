#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([PHLAWD], [3.0], [eebsmith@umich.edu])
AC_CONFIG_SRCDIR([fasta_util.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [main])
AC_CHECK_LIB([sqlite3], [sqlite3_trace])
#THESE NEED TO BE WRAPPED IN WITH THE MAKE
#AC_CHECK_LIB([sqlitewrapped],[_ZN5Query12GetLastQueryEv])
#AC_SEARCH_LIBS([_ZN5Query12GetLastQueryEv], [sqlitewrapped], [have_adolc=yes])

AC_CHECK_PROG([muscle_var],[muscle],[yes],[no]) 
#change to WARN when wrapped into make
if test "x${muscle_var}" = xno; then
   AC_MSG_ERROR([
--------------------------------------------------------
The muscle program is required to be in the path
--------------------------------------------------------])
fi
AC_CHECK_PROG([mafft_var],[mafft],[yes],[no])
#change to WARN when wrapped into make
if test "x${mafft_var}" = xno; then
   AC_MSG_ERROR([
--------------------------------------------------------
The mafft program is required to be in the path
--------------------------------------------------------])
fi

#############################
# Platform specific setup
#############################

ARCH=`uname`
HOST="linux"
if [[ $ARCH == "Darwin" ]]; then
   HOST="mac"
fi
AC_SUBST(HOST)


# Checks for header files.
AC_CHECK_HEADERS([float.h malloc.h stdint.h stdlib.h string.h syslog.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memset mkdir sqrt strerror strtol])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
