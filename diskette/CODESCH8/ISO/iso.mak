SWI= -C -c -O3 -g   -qarch=pwr2 -qposition=append
iso:  isoma.o isonn.o isohn.o isoph.o  isotr.o ../../NCARFFT/ffts.o  isosp.o isotn.o isoindip.o isoio.o isointu.o matfun.o isoinpefftps.o
	xlf -o iso isoma.o isonn.o isohn.o isoph.o isotr.o ../../NCARFFT/ffts.o  isosp.o isotn.o isoindip.o isoio.o isointu.o matfun.o isoinpefftps.o
isoma.o: isoma.f param.f 
	xlf $(SWI)  isoma.f
isonn.o: isonn.f param.f   
	xlf $(SWI)  isonn.f
isohn.o: isohn.f param.f
	xlf $(SWI)  isohn.f
isoph.o: isoph.f param.f
	xlf $(SWI)  isoph.f
isotr.o: isotr.f param.f
	xlf $(SWI)  isotr.f
isosp.o: isosp.f param.f
	xlf $(SWI)  isosp.f
isotn.o: isotn.f param.f
	xlf $(SWI)  isotn.f
isoindip.o: isoindip.f param.f
	xlf $(SWI)  isoindip.f
isoinpefftps.o: isoinpefftps.f param.f
	xlf $(SWI)  isoinpefftps.f
isointu.o: isointu.f param.f
	xlf $(SWI)  isointu.f
isoio.o: isoio.f param.f
	xlf $(SWI)  isoio.f
matfun.o: matfun.f param.f
	xlf $(SWI)  matfun.f
