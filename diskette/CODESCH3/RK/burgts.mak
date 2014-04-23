SWI= -c -O3    -qarch=pwr2 
burgts: burgts.o  
	xlf -o burgts burgts.o  
burgts.o: burgts.f param.f
	xlf $(SWI) burgts.f
