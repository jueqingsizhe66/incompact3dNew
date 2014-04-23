SWI= -C -c -O3 -g  -qarch=pwr2
statiso:  staisoma.o staisonn.o staisosp.o staisopd.o staisoou.o ../../NCARFFT/ffts.o  
	xlf -o statiso staisoma.o staisonn.o staisosp.o staisopd.o staisoou.o ../../NCARFFT/ffts.o 
staisoma.o: staisoma.f param.f 
	xlf $(SWI)  staisoma.f
staisonn.o: staisonn.f param.f 
	xlf $(SWI)  staisonn.f
staisosp.o: staisosp.f param.f 
	xlf $(SWI)  staisosp.f
staisopd.o: staisopd.f param.f 
	xlf $(SWI)  staisopd.f
staisoou.o: staisoou.f param.f 
	xlf $(SWI)  staisoou.f
