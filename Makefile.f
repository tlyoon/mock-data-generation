# Makefile to compile Quantum Monte Carlo simulation for
# Computational Physics

# COMPHY directory path:
#COMPHY	=	$(HOME)/comphy
#LIBS	=	-L$(COMPHY)/lib -lran
FFLAGS	=	-O -c
LDFLAGS	=	-O
FC	=	gfortran

gmm01f: gmm01f.o
	$(FC) -o gmm01f $(LDFLAGS) gmm01f.o $(LIBS)

gmm01f.o: gmm01f.f gmm01f.par
	$(FC) $(FFLAGS) gmm01f.f
clean:
	rm -f *.o gmm01f *.out
