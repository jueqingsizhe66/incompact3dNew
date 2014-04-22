FC = gfortran
modFFLAGS= -c -O3 
FFLAGS = -g -Wall -O3
capro: capro.o pfishck.o psncs.o pcssn.o pftpe.o p2fft.o  ../../NCARFFT/fishpack.o ../../NCARFFT/ffts.o
	$(FC) -o capro capro.o pfishck.o psncs.o pcssn.o pftpe.o p2fft.o  ../../NCARFFT/fishpack.o ../../NCARFFT/ffts.o
capro.o: capro.f param.f
	$(FC) $(modFFLAGS) capro.f
psncs.o: psncs.f param.f
	$(FC) $(modFFLAGS) psncs.f
pcssn.o: pcssn.f param.f
	$(FC) $(modFFLAGS) pcssn.f
p2fft.o: p2fft.f param.f
	$(FC) $(modFFLAGS) p2fft.f
pftpe.o: pftpe.f param.f
	$(FC) $(modFFLAGS) pftpe.f
pfishck.o: pfishck.f param.f
	$(FC) $(modFFLAGS) pfishck.f
../../NCARFFT/fishpack.o: ../../NCARFFT/fishpack.f param.f
	$(FC) $(modFFLAGS) ../../NCARFFT/fishpack.f
../../NCARFFT/ffts.o: ../../NCARFFT/ffts.f param.f
	$(FC) $(modFFLAGS) ../../NCARFFT/ffts.f

.PHONY: clean

clean:
	rm *.o capro
