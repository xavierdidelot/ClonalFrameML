# Make file for ClonalFrameML
CC = g++
CFLAGS = -O3 -I./ -I./myutils -I./coalesce
LDFLAGS = 
OBJECTS = main.o
HEADERS = main.h brent.h powell.h version.h

.PHONY: clean version

all: version ClonalFrameML

ClonalFrameML: $(OBJECTS)
	$(CC) $(LDFLAGS) -o ClonalFrameML $(OBJECTS)

main.o: main.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c -o main.o main.cpp

version:
	echo "#define ClonalFrameML_SVNRevision \c" > version.h; svn info | grep "Revision: " | cut -f 2 -d' ' >> version.h

clean:
	rm $(OBJECTS)
