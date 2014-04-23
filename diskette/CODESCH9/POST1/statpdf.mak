SWI= -C -c -O3 -g  -qarch=pwr2
stapdfbook: statpdfl.o 
	xlf  -o stapdfbook statpdfl.o 
statpdfl.o: statpdfl.f 
	xlf $(SWI) statpdfl.f
 
