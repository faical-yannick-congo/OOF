#######################################################################
# $Header: /users/langer/FE/CVSoof/OOF2/OOFSWIG/SWIG/Makefile.in,v 1.2 2014/06/26 16:05:05 lck Exp $
# Simplified Wrapper and Interface Generator (SWIG)
# 
# Makefile for version 1.0 Final
# Dave Beazley
# August 1, 1996
#
# This makefile is now mostly constructed by ./configure.
#
# $Log: Makefile.in,v $
# Revision 1.2  2014/06/26 16:05:05  lck
# *** empty log message ***
#
# Revision 1.1.1.1  1999/02/28 02:00:51  beazley
# Swig1.1
#
# Revision 1.1  1996/08/12 01:55:02  dmb
# Initial revision
#
#######################################################################

#.KEEP_STATE:


srcdir = .


# Set your C++ compiler here.   g++ works on most machines,
# but you might have to change it depending on your installation.
#
CC      = g++

#
# Set the prefix below to indicate where you want SWIG to install its
# files.   Normally this is /usr/local
# 

prefix     = /home/fyc/Documents/Projects/OOFMigrate/_CVS/OOF2/OOFSWIG/BUILD

# Location of the SWIG library.  Is normally put in /usr/local/lib/swig_lib
# The SWIG library contains configuration files and library modules
# so you should install it someplace where it can be easily accessed.

SWIG_LIB = $(prefix)/lib/swig_lib

# YACC parser.  Use bison by default.  if this doesn't work, switch
# it over to yacc.   If that still doesn't work, let me know...

YACC    = bison -y 

# Comment out the following line if you're on an SGI or don't have ranlib!
RANLIB  = ranlib
AR      = ar

########################################################################
# Normally, you shouldn't have to change anything below this point     #
########################################################################

LIBOBJS = main.o scanner.o symbol.o include.o types.o parms.o emit.o newdoc.o ascii.o \
          html.o latex.o cplus.o lang.o hash.o sstring.o wrapfunc.o getopt.o comment.o \
          typemap.o naming.o

LIBSRCS = main.cxx scanner.cxx symbol.cxx include.cxx types.cxx parms.cxx emit.cxx \
	  newdoc.cxx ascii.cxx html.cxx latex.cxx cplus.cxx lang.cxx hash.cxx \
          sstring.cxx wrapfunc.cxx getopt.cxx comment.cxx typemap.cxx naming.cxx

LIBHEADERS = internal.h ../Include/swig.h latex.h ascii.h html.h nodoc.h
LIB      = ../libswig.a
PARSER   = parser.y
INCLUDE  = -I../Include
CFLAGS   = -g -O2 -DSWIG_LIB='"$(SWIG_LIB)"' -DSWIG_CC='"$(CC)"' -DSTDC_HEADERS=1 -DHAVE_ALLOCA_H=1 -DHAVE_ALLOCA=1 -DSWIG_DOC=ASCII -DHAVE_LIBDL=1 
SHELL  = /bin/sh

#
#
#
# Rules for creation of a .o file from .cxx
.SUFFIXES: .cxx
.cxx.o: 
	$(CC) $(INCLUDE) $(CFLAGS) -c -o $*.o $<

all:    $(LIB) 

$(LIB): parser.o $(LIBOBJS) 
	@echo "Building library"
	$(AR) cr $(LIB) $(LIBOBJS) parser.o
	$(RANLIB) $(LIB)

parser.o: parser.cxx $(LIBHEADERS)
	$(CC) $(INCLUDE) $(CFLAGS) parser.cxx -c -o parser.o

parser.cxx: $(PARSER)
	$(YACC) -d parser.y
	@cp y.tab.h parser.h
	@cp y.tab.c parser.cxx

parser::  
	@cp y.tab.c.bison parser.cxx
	@cp y.tab.h.bison parser.h
	@cp y.tab.h.bison y.tab.h
	$(CC) $(CFLAGS) parser.cxx -c -o parser.o

Makefile:	$(srcdir)/Makefile.in ../config.status
		(cd ..; CONFIG_FILES=SWIG/Makefile $(SHELL) config.status)

.PRECIOUS:	Makefile

clean::
	rm -f *.o libswig.a y.tab.c y.tab.h

nuke::
	rm -f Makefile *~ #* core a.out

wc::	
	wc $(LIBSRCS) *.h parser.y

# I use this to build final releases from the RCS directory.   Probably
# not very useful to anyone else...

co::
	co $(LIBSRCS) $(LIBHEADERS) $(PARSER)



