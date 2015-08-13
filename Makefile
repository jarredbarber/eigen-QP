# Simple build script
CC=g++

CFLAGS=-O3 -I. -std=c++14

all: eigen-qp.o driver.o
	$(CC) $(CFLAGS) -o driver eigen-qp.o driver.o

eigen-qp.o: eigen-qp.hpp eigen-qp.cpp
	$(CC) $(CFLAGS) -c eigen-qp.cpp

driver.o: driver.cpp
	$(CC) $(CFLAGS) -c driver.cpp

clean:
	rm -f driver *.o