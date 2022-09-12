CC=gcc
CFLAGS=-I./include/

spectrum: spectrum.o windows.o
	$(CC) -o spectrum spectrum.o windows.o -lm -lkissfft-float

clean:
	rm -f *.o spectrum
