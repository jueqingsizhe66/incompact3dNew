SWI= -c -O3   -qarch=pwr2 
visiso:  visiso.o ../../NCARFFT/ffts.o  
	xlf -o visiso visiso.o ../../NCARFFT/ffts.o 
visiso.o: visiso.f param.f 
	xlf $(SWI)  visiso.f
