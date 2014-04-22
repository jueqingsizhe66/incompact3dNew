SWI= -c -O3   -qarch=pwr2
turb2dcu: turb2dma.o turb2dnn.o turb2dtn.o turb2din.o turb2dph.o  ../../NCARFFT/ffts.o 
	xlf -o turb2dcu turb2dma.o turb2dnn.o turb2dtn.o turb2din.o turb2dph.o ../../NCARFFT/ffts.o 
turb2dma.o: turb2dma.f
	xlf $(SWI) turb2dma.f
turb2dnn.o: turb2dnn.f
	xlf $(SWI) turb2dnn.f
turb2dtn.o: turb2dtn.f
	xlf $(SWI) turb2dtn.f
turb2din.o: turb2din.f
	xlf $(SWI) turb2din.f
turb2dph.o: turb2dph.f
	xlf $(SWI) turb2dph.f
