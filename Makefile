# Simple build script
CC=g++

CFLAGS=-O3 -I. -std=c++14 -lboost_system -lboost_timer

all: driver.cpp eigen-qp.hpp
	$(CC) $(CFLAGS) -o driver driver.cpp

clean:
	rm -f driver *.o