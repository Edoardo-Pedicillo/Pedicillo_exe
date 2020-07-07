CC = g++
CFLAGS = -Wall -O3 --std=c++11

lab9.exe : lab9.o random.o classi.o
	$(CC) random.o lab9.o classi.o -o lab9.exe

lab9.o : lab9.cpp random.h classi.h
	$(CC) -c lab9.cpp -o lab9.o $(CFLAGS)

classi.o : classi.cpp classi.h random.h
	$(CC) -c classi.cpp -o classi.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o *.dat main1.exe main2.exe main3.exe seed.out
