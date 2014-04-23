SWI= -C -c -O3 -g  -qarch=pwr2
chatra: chatrann.o chatrats.o chatraio.o chatratu.o  chatrasp.o ../../NCARFFT/ffts.o
	xlf  -o chatra chatrann.o chatrats.o chatraio.o  chatratu.o  chatrasp.o  ../../NCARFFT/ffts.o
chatrats.o: chatrats.f param.f
	xlf $(SWI) chatrats.f
chatrann.o: chatrann.f param.f
	xlf $(SWI) chatrann.f
chatraio.o: chatraio.f param.f
	xlf $(SWI) chatraio.f
chatrasp.o :chatrasp.f param.f
	xlf $(SWI) chatrasp.f 
chatratu.o :chatratu.f param.f
	xlf $(SWI) chatratu.f 
 
