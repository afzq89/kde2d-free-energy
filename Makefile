OBJECTS = kernel.o
FCOMPLR = gfortran -O2
FCFLAGS =  -c   -w -Wall -ffree-form

kernel : $(OBJECTS) Makefile
	$(FCOMPLR) $(LDFLAGS) $(OBJECTS) -o $@


#  specific dependencies

kernel.o: kernel.f
	$(FCOMPLR) $(FCFLAGS) kernel.f

clean :
	rm -f *.o *~ *.mod
