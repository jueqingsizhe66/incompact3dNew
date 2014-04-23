SWI= -c -p -g -O5
trasco:  trasco.o 
	gfortran -p -g  -o trasco trasco.o 
trasco.o: trasco.f param.f 
	gfortran $(SWI)  trasco.f
