CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

%.o : %.cpp
		g++ -Wall -c $< ${INCS}

qmc1d: qmc1d.o
		g++ -O3 -Wall -o $@ $^ ${LIBS}

lab8.exe : lab8.o random.o
	$(CC) random.o lab8.o -o lab8.exe
lab8.o : lab8.cpp lab8.h
	$(CC) -c lab8.cpp -o lab8.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
		rm *.o qmc1d potential.dat kinetic.dat probability.dat lab8.exe
