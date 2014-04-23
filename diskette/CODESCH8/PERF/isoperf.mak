SWI= -c -p -O3   -qarch=pwr2
isoper:  isoma.o isonn.o isohn.o isoph.o  isotr.o matfun.o ffts.o  isosp.o isotn.o isointu.o isoindipn.o isoio.o 
	xlf -p -o isoper isoma.o isonn.o isohn.o isoph.o isotr.o matfun.o ffts.o isosp.o isotn.o isointu.o isoindipn.o isoio.o
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
isointu.o: isointu.f param.f
	xlf $(SWI)  isointu.f
isoindipn.o: isoindipn.f param.f
	xlf $(SWI)  isoindipn.f
isoio.o: isoio.f param.f
	xlf $(SWI)  isoio.f
ffts.o: ffts.f 
	xlf $(SWI)  ffts.f
matfun.o: matfun.f 
	xlf $(SWI)  matfun.f
