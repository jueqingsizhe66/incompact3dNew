SWI= -c -O3 -qarch=pwr2
tr3dnunu: tr3dnuco.o tr3dnutn.o  tr3dnuma.o tr3dnuphnu.o  tr3dnuinri.o tr3dnuintr.o   tr3dnuhn.o tr3dnuio.o tr3dnutu.o tr3dnutr.o  tr3dnuge.o /u/orlandi/ncarfft/sffts.o ../../NCARFFT/fishpack.o
	xlf  -o tr3dnunu tr3dnuco.o tr3dnutn.o tr3dnuma.o tr3dnuphnu.o  tr3dnuinri.o tr3dnuintr.o   tr3dnuhn.o tr3dnuio.o  tr3dnutu.o tr3dnutr.o tr3dnuge.o /u/orlandi/ncarfft/sffts.o ../../NCARFFT/fishpack.o
tr3dnuhn.o: tr3dnuhn.f param.f
	xlf $(SWI) tr3dnuhn.f
tr3dnutn.o: tr3dnutn.f param.f
	xlf $(SWI) tr3dnutn.f
tr3dnuma.o: tr3dnuma.f param.f
	xlf $(SWI) tr3dnuma.f
tr3dnuco.o: tr3dnuco.f param.f
	xlf $(SWI) tr3dnuco.f
tr3dnuio.o: tr3dnuio.f param.f
	xlf $(SWI) tr3dnuio.f
tr3dnuintr.o: tr3dnuintr.f param.f
	xlf $(SWI) tr3dnuintr.f
tr3dnuinri.o: tr3dnuinri.f param.f
	xlf $(SWI) tr3dnuinri.f
tr3dnuphnu.o :tr3dnuphnu.f param.f
	xlf $(SWI) tr3dnuphnu.f 
tr3dnutu.o :tr3dnutu.f param.f
	xlf $(SWI) tr3dnutu.f 
tr3dnutr.o :tr3dnutr.f paramdi.f
	xlf $(SWI) tr3dnutr.f 
tr3dnuge.o :tr3dnuge.f param.f
	xlf $(SWI) tr3dnuge.f 
 
