# Makefile to compile Quantum Monte Carlo simulation for
# Computational Physics

# COMPHY directory path:
#COMPHY	=	$(HOME)/comphy
#LIBS	=	-L$(COMPHY)/lib -lran
FFLAGS	=	-O -c
LDFLAGS	=	-O
FC	=	gfortran

gmm01s: gmm01s.o
	$(FC) -o gmm01s $(LDFLAGS) gmm01s.o $(LIBS)

gmm01s.o: gmm01s.f gmm01f.par
	$(FC) $(FFLAGS) gmm01s.f
clean:
	rm -f *.o gmm01s *.out
