SWI= -c -O3     -qarch=pwr2
pistat: pistnn.o pistts.o pistio.o pisttu.o  pistco.o ../../NCARFFT/ffts.o
	xlf  -o pistat pistnn.o pistts.o pistio.o  pisttu.o  pistco.o  ../../NCARFFT/ffts.o 
pistts.o: pistts.f param.f
	xlf $(SWI) pistts.f
pistnn.o: pistnn.f param.f
	xlf $(SWI) pistnn.f
pistio.o: pistio.f param.f
	xlf $(SWI) pistio.f
pistco.o :pistco.f param.f
	xlf $(SWI) pistco.f 
pisttu.o :pisttu.f param.f
	xlf $(SWI) pisttu.f 
 
