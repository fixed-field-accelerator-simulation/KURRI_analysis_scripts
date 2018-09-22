

.POSIX:

tomosubs_v2.o: tomosubs_v2.f90 long_v2.o nrstuff.o
	$(FC) $(FFLAGS) $(CFLAGS) -c tomosubs_v2.f90

long_v2.o: long_v2.f90 nrstuff.o
	$(FC) $(FFLAGS) $(CFLAGS) -c long_v2.f90

nrstuff.o: nrstuff.f90
	$(FC) $(FFLAGS) $(CFLAGS) -c nrstuff.f90

clean:
	@rm -f  *.o *.mod
