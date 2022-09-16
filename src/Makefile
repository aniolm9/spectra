CC=gcc
CFLAGS=-I./include/ -O3 -ffast-math -Wall
LIBS=

spectrum: spectrum.o windows.o
	$(CC) $(CFLAGS) -o spectrum spectrum.o windows.o $(LIBS) -lm -lkissfft-float -lgsl

clean:
	rm -f *.o spectrum
