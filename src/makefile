# Makefile for ClonalFrameML
CC = g++
CFLAGS += -O3
OBJECTS = main.o
HEADERS = main.h brent.h powell.h

.PHONY: clean 

all: ClonalFrameML

ClonalFrameML: $(OBJECTS)
	$(CC) $(LDFLAGS) -o ClonalFrameML $(OBJECTS)

main.o: main.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c -o main.o main.cpp

clean:
	rm -f $(OBJECTS)
