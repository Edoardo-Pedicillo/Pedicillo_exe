CC = g++
CFLAGS = -Wall -O3 --std=c++11



main1.exe : main1.o random.o funzioni.o classi.o
	$(CC) random.o main1.o funzioni.o classi.o -o main1.exe

main1.o : main1.cpp funzioni.h random.h classi.h
	$(CC) -c main1.cpp -o main1.o $(CFLAGS)

funzioni.o : funzioni.cpp funzioni.h random.h classi.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

classi.o : classi.cpp classi.h random.h
	$(CC) -c classi.cpp -o classi.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o *.dat main1.exe seed.out
