SWI= -c -O3    -qarch=pwr2 
capro: capro.o pfishck.o psncs.o pcssn.o pftpe.o p2fft.o p2fftess.o ../../NCARFFT/fishpack.o ../../NCARFFT/ffts.o
	xlf -o capro capro.o pfishck.o psncs.o pcssn.o pftpe.o p2fft.o p2fftess.o ../../NCARFFT/sfishpack.o ../../NCARFFT/ffts.o -l esslp2
capro.o: capro.f param.f
	xlf $(SWI) capro.f
psncs.o: psncs.f param.f
	xlf $(SWI) psncs.f
pcssn.o: pcssn.f param.f
	xlf $(SWI) pcssn.f
p2fft.o: p2fft.f param.f
	xlf $(SWI) p2fft.f
pftpe.o: pftpe.f param.f
	xlf $(SWI) pftpe.f
pfishck.o: pfishck.f param.f
	xlf $(SWI) pfishck.f
p2fftess.o: p2fftess.f param.f
	xlf $(SWI) p2fftess.f
../../NCARFFT/fishpack.o: ../../NCARFFT/fishpack.f param.f
	xlf $(SWI) ../../NCARFFT/fishpack.f
../../NCARFFT/ffts.o: ../../NCARFFT/ffts.f param.f
	xlf $(SWI) ../../NCARFFT/ffts.f
