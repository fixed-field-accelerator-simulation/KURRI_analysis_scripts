# Makefile for Tomography compilation

# MacOSX

# Intel settings
# FC=ifort
# CC=icc
# FFLAGS=-O -qopenmp -w -xHOST -fast

# gcc settings
 FC=gfortran
 CC=gcc
 FFLAGS=-O -w

# For debugger, add -g

#--------------------------------------------------------------

 SRCDIR=./tomo_v2
 MODS=./modules


.POSIX:

all :
	cd $(MODS) ; $(MAKE) FC="$(FC)" FFLAGS="$(FFLAGS)" CFLAGS="$(CFLAGS)"
	cd $(SRCDIR) ; $(MAKE) FC="$(FC)" FFLAGS="$(FFLAGS)" CFLAGS="$(CFLAGS)"
	@mv $(SRCDIR)/tomo ./

clean :
	@rm tomo
	cd $(SRCDIR) ; $(MAKE) clean
	cd $(MODS) ; $(MAKE) clean

#docs :
#	cd manual ; $(MAKE) ; $(MAKE) clean
