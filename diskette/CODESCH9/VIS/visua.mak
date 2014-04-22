SWI= -c -O3   -qarch=pwr2 
visua:  visua.o ../../NCARFFT/ffts.o  
	xlf -o visua visua.o ../../NCARFFT/ffts.o 
visua.o: visua.f param.f 
	xlf $(SWI)  visua.f
