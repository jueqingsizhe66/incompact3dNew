FC=gfortran
FFLAGS=-O3
split : ./f90split.f90 ./output/cmd.sh
	$(FC) $(FFLAGS) $< -o $@
	mv $@ ./output/
	cd ./output/ &&sh cmd.sh

.PHONY: clean

