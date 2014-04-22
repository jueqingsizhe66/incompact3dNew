SWI= -c -O3    -qarch=pwr2 
burgexim: burgexim.o  
	xlf -o burgexim burgexim.o  
burgexim.o: burgexim.f param.f
	xlf $(SWI) burgexim.f
