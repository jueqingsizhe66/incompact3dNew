SWI= -C -c -O3 -g  -qarch=pwr2
filterbook:  filterbook.o isoph.o /u/orlandi/ncarfft/sffts.o  
	xlf -o filterbook filterbook.o isoph.o /u/orlandi/ncarfft/sffts.o  
filterbook.o: filterbook.f param.f 
	xlf $(SWI)  filterbook.f
isoph.o: isoph.f param.f
	xlf $(SWI)  isoph.f
