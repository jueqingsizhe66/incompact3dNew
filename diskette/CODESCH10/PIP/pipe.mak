SWI= -c -O3 -qarch=pwr2
piprot:  pipco.o pipge.o piphn.o piptn.o pipphnc.o pipio.o pipintur.o pipindip.o piptr.o piptu.o pipma.o pipti.o ../../NCARFFT/ffts.o
	xlf  -o piprot pipco.o pipge.o piphn.o piptn.o pipphnc.o pipio.o pipintur.o pipindip.o piptr.o piptu.o pipma.o pipti.o ../../NCARFFT/ffts.o
pipco.o: pipco.f param.f   
	xlf $(SWI)  pipco.f
piphn.o: piphn.f param.f   
	xlf $(SWI)  piphn.f
pipti.o: pipti.f param.f
	xlf $(SWI)  pipti.f
piptn.o: piptn.f param.f
	xlf $(SWI)  piptn.f
pipphnc.o: pipphnc.f param.f
	xlf $(SWI)  pipphnc.f
pipio.o: pipio.f param.f
	xlf $(SWI)  pipio.f
pipintur.o: pipintur.f param.f
	xlf $(SWI)  pipintur.f
pipindip.o: pipindip.f param.f
	xlf $(SWI)  pipindip.f
piptr.o: piptr.f param.f
	xlf $(SWI)  piptr.f
piptu.o: piptu.f param.f
	xlf $(SWI)  piptu.f
pipma.o: pipma.f param.f
	xlf $(SWI)  pipma.f
pipge.o: pipge.f param.f
	xlf $(SWI)  pipge.f
