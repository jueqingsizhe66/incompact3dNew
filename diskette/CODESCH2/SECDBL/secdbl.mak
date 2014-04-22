SWI= -c -O3    -qarch=pwr2 
secdbls: secdbls.o  
	xlf -o secdbls secdbls.o  
secdbls.o: secdbls.f param.f
	xlf $(SWI) secdbls.f
