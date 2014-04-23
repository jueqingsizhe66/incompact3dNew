SWI= -c -O3   -qarch=pwr2
visuanun: vistrinun.o 
	xlf  -o visuanun vistrinun.o 
vistrinun.o: vistrinun.f param.f
	xlf $(SWI) vistrinun.f
 
