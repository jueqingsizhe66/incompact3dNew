SWI= -c -O3   -qarch=pwr2
chadses:  chadsma.o chadsnn.o chadshn.o chadstn.o chadstinv.o  chadsphnc.o chadsou.o chadsin.o chadstr.o chadsru.o  ../../NCARFFT/ffts.o
	xlf  -o chadses chadsma.o chadshn.o chadsnn.o chadstn.o chadstinv.o chadsphnc.o chadsou.o chadsin.o chadstr.o chadsru.o ../../NCARFFT/ffts.o
chadsma.o: chadsma.f param.f 
	xlf $(SWI)  chadsma.f
chadsnn.o: chadsnn.f param.f   
	xlf $(SWI)  chadsnn.f
chadshn.o: chadshn.f param.f   
	xlf $(SWI)  chadshn.f
chadstn.o: chadstn.f param.f
	xlf $(SWI)  chadstn.f
chadstinv.o: chadstinv.f param.f
	xlf $(SWI)  chadstinv.f
chadsphnc.o: chadsphnc.f param.f
	xlf $(SWI)  chadsphnc.f
chadsou.o: chadsou.f param.f
	xlf $(SWI)  chadsou.f
chadsin.o: chadsin.f param.f
	xlf $(SWI)  chadsin.f
chadstr.o: chadstr.f param.f
	xlf $(SWI)  chadstr.f
chadsru.o: chadsru.f param.f
	xlf $(SWI)  chadsru.f
