CFLAGS1 = -O
CFLAGS2 = -W
CC = g++

MC: main.o calc.o geometry.o xsection.o neutron.o
	$(CC) $(CFLAGS1) $(CFLAGS2) -o MC main.o calc.o geometry.o xsection.o neutron.o

main.o: main.cpp calc.h geometry.h xsection.h neutron.h
	$(CC) $(CFLAGS) -c main.cpp

xsection.o: xsection.cpp xsection.h neutron.h calc.h
	$(CC) $(CFLAGS) -c xsection.cpp

geometry.o: geometry.cpp neutron.h
	$(CC) $(CFLAGS) -c geometry.cpp

neutron.o: neutron.cpp neutron.h calc.h
	$(CC) $(CFLAGS) -c neutron.cpp

calc.o: calc.cpp calc.h
	$(CC) $(CFLAGS) -c calc.cpp

clean: 
	rm *.o