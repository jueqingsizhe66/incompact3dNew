SWI= -C -c -O3   -qarch=pwr
psomgbcn: psomgbcma.o psomgbctn.o psomgbcbc.o psomgbcps.o psomgbcnn.o psomgbcin.o matfun.o ../../NCARFFT/ffts.o 
	xlf  -o psomgbcn psomgbcma.o psomgbctn.o psomgbcbc.o  psomgbcps.o psomgbcnn.o psomgbcin.o matfun.o ../../NCARFFT/ffts.o  
psomgbcma.o: psomgbcma.f param.f
	xlf $(SWI) psomgbcma.f
psomgbctn.o: psomgbctn.f param.f
	xlf $(SWI) psomgbctn.f
psomgbcbc.o: psomgbcbc.f param.f
	xlf $(SWI) psomgbcbc.f
psomgbcnn.o: psomgbcnn.f param.f
	xlf $(SWI) psomgbcnn.f
psomgbcin.o: psomgbcin.f param.f
	xlf $(SWI) psomgbcin.f
psomgbcps.o: psomgbcps.f param.f
	xlf $(SWI) psomgbcps.f
matfun.o: matfun.f param.f
	xlf $(SWI) matfun.f
../../NCARFFT/ffts.o: ../../NCARFFT/ffts.f 
	xlf $(SWIf) ../../NCARFFT/ffts.f
