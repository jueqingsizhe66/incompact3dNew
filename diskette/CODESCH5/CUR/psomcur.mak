SWI= -c -O3   -qarch=pwr2
psomcur: psomcurco.o psomcurnn.o psomcurtn.o psomcurbc.o psomcurma.o psomcurmg.o psomcurit.o psomcurin.o
	xlf  -o psomcur psomcurma.o psomcurnn.o psomcurtn.o psomcurbc.o psomcurco.o   psomcurmg.o  psomcurit.o psomcurin.o
psomcurnn.o: psomcurnn.f param.f
	xlf $(SWI) psomcurnn.f
psomcurbc.o: psomcurbc.f param.f
	xlf $(SWI) psomcurbc.f
psomcurtn.o: psomcurtn.f param.f
	xlf $(SWI) psomcurtn.f
psomcurco.o :psomcurco.f param.f
	xlf $(SWI) psomcurco.f 
psomcurma.o :psomcurma.f param.f
	xlf $(SWI) psomcurma.f 
psomcurmg.o :psomcurmg.f param.f
	xlf $(SWI) psomcurmg.f 
psomcurit.o :psomcurit.f param.f
	xlf $(SWI) psomcurit.f 
psomcurin.o :psomcurin.f param.f
	xlf $(SWI) psomcurin.f 
 
