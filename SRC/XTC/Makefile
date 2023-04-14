# Define ARCH to be one the of architectures found in the conf directory.
# The names are the same as the ones used for compiling pvm, so you might
# want to use the value of PVM_ARCH.
ARCH	= SUN3

# Set FORTRAN compiler and flags for ARCH
F77     = ifort
FFLAGS  = -O3

# Set C compiler and flags for ARCH
CC      = cc
CFLAGS 	= 

# Set HASRANLIB to t if your system needs ranlib
HASRANLIB = t	

# Set LIBS to libraries to link your c program with
# for SUN Solaris 2:
# LIBS = -lnsl

SHELL	= /bin/sh
# set RM to rm -f if you want to remove the intermediate C source file
# generated for the FORTRAN interface library
RMCMD 	= rm -f
#RMCMD	= /bin/true
M4 	= m4
M4FILE	= conf/$(ARCH).m4

LFOBS 	= libxdrf.o xtciof.o
LOBS	= xdr3dfcoord.o

default: libxdrf.a 


libxdrf.a:  $(LFOBS) ftocstr.o
	ar cr libxdrf.a $?
	case x$(HASRANLIB) in xt ) echo ranlib; ranlib libxdrf.a ;; esac

clean:
	rm -f $(LFOBS) $(LOBS) ftocstr.o libxdrf.a

tidy:
	rm -f $(LOBS) $(LFOBS)

ftocstr.o: ftocstr.c
	$(CC) $(CFLAGS) -c ftocstr.c

libxdrf.o:	libxdrf.m4 $(M4FILE)
	$(M4) $(M4FILE) libxdrf.m4 > libxdrf.c
	$(CC) $(CFLAGS) -c libxdrf.c
	$(RMCMD) libxdrf.c

xtciof.o:	xtciof.f
	$(F77) $(FFLAGS) -c xtciof.f

conf/.m4:
	@echo "ERROR: you didn't set ARCH in the Makefile"
