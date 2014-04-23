#
FC = gfortran
FFLAGS = -g -Wall -O3
secd : ./secdml.f ./param.f
	$(FC) $(FFLAGS) $< -o $@

.PHONY : clean

clean: 
	rm -rf ./errdis* ./fort* ./accuracy.plo ./accuracy.plo ./funct.plo secd
