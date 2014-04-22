FC = gfortran
FFLAGS = -g -Wall -O3 --free-form
modFFLAGS = -c  --free-form


secdbl : ./secdbls.f90 ./secdmod.f90
	$(FC) $(FFLAGS) -o $@ $^
secd.mod : ./secdmod.f90
	$(FC) $(modFFLAGS) @ $^
#secdmo.mod : secdmo.f
#	$(FC) $(modFFLAGS) @ $<

.phony: clean

clean:
	rm -rf secdbl *.mod *.plo  fort*
