FC = gfortran
FFLAGS = -g -Wall -O3 --free-form
modFFLAGS = -c  --free-form


CORD : trascor.f90 ./trascomod.f90
	$(FC) $(FFLAGS) -o $@ $^
tramod.mod : ./trascomod.f90
	$(FC) $(modFFLAGS) @ $^
#secdmo.mod : secdmo.f
#	$(FC) $(modFFLAGS) @ $<

.phony: clean

clean:
	rm -rf CORD tramod.mod *.plo 
