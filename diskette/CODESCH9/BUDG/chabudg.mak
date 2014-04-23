SWI= -C -c -O3 -g  -qarch=pwr2
chabudgbook: chabudgnn.o chabudgts.o chabudgio.o chabudgtu.o  
	xlf  -o chabudgbook chabudgnn.o chabudgts.o chabudgio.o  chabudgtu.o  
chabudgts.o: chabudgts.f param.f
	xlf $(SWI) chabudgts.f
chabudgnn.o: chabudgnn.f param.f
	xlf $(SWI) chabudgnn.f
chabudgio.o: chabudgio.f param.f
	xlf $(SWI) chabudgio.f
chabudgtu.o :chabudgtu.f param.f
	xlf $(SWI) chabudgtu.f 
 
