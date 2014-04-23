SWI= -c -O3   -qarch=pwr2
visua: vispirot.o 
	xlf  -o visua vispirot.o 
vispirot.o: vispirot.f param.f
	xlf $(SWI) vispirot.f
 
