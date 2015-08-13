# Simple build script
CC=cc

CFLAGS=-O3 -I.

all: eigen-qp.o driver.o
	$(CC) $(CFLAGS) -o driver eigen-qp.o driver.o

eigen-qp.o: eigen-qp.hpp eigen-qp.cpp
	$(CC) $(CFLAGS) -c eigen-qp.cpp

driver.o: driver.cpp
	$(CC) $(CFLAGS) -c driver.cpp