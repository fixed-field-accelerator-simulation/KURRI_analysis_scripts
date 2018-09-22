#
MOD = ../modules
OBJS = $(MOD)/tomosubs_v2.o  $(MOD)/long_v2.o $(MOD)/nrstuff.o


.POSIX:

all: tomo

tomo: tomo_v2.f90 $(OBJS)
	$(FC) $(FFLAGS) $(CFLAGS) tomo_v2.f90 -I$(MOD) $(OBJS) -o $@

clean:
	@rm -f tomo *.o
