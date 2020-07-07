CC = g++
CFLAGS = -Wall -O3 --std=c++11

lab10.exe : lab10.o random.o classi.o
	mpicxx random.o lab10.o classi.o -o lab10.exe

lab10.o : lab10.cpp random.h classi.h
	mpicxx -c lab10.cpp -o lab10.o $(CFLAGS)

classi.o : classi.cpp classi.h random.h
	$(CC) -c classi.cpp -o classi.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o *.dat seed.out
