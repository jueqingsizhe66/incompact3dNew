SWI= -c -O3 
chalesnc:  chalesma.o chalesnn.o chaleshn.o chalestn.o chalesphnc.o chalesou.o chalesin.o chalestr.o chalesru.o chalessma.o chalesdyn.o  ../../NCARFFT/ffts.o
	gfortran  -o chalesnc chalesma.o chaleshn.o chalesnn.o chalestn.o chalesphnc.o chalesou.o chalesin.o chalestr.o chalesru.o chalessma.o chalesdyn.o -l esslp2 ../../NCARFFT/ffts.o
chalesma.o: chalesma.f param.f 
	gfortran $(SWI)  chalesma.f
chalesnn.o: chalesnn.f param.f   
	gfortran $(SWI)  chalesnn.f
chaleshn.o: chaleshn.f param.f   
	gfortran $(SWI)  chaleshn.f
chalestn.o: chalestn.f param.f
	gfortran $(SWI)  chalestn.f
chalesphnc.o: chalesphnc.f param.f
	gfortran $(SWI)  chalesphnc.f
chalesou.o: chalesou.f param.f
	gfortran $(SWI)  chalesou.f
chalesin.o: chalesin.f param.f
	gfortran $(SWI)  chalesin.f
chalestr.o: chalestr.f param.f
	gfortran $(SWI)  chalestr.f
chalesru.o: chalesru.f param.f
	gfortran $(SWI)  chalesru.f
chalessma.o: chalessma.f param.f
	gfortran $(SWI)  chalessma.f
chalesdyn.o: chalesdyn.f param.f
	gfortran $(SWI)  chalesdyn.f
