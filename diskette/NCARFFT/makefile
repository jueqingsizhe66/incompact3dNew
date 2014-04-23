FC = gfortran
FFLAGS = -W -O2 -c
All : fishpack.o  ffts.o

fishpack.o : fishpack.f
	$(FC) $(FFLAGS) $^
ffts.o : ffts.f
	$(FC) $(FFLAGS) $^
