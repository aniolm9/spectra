CC=gcc
CFLAGS=-I./include/

spectrum: spectrum.o windows.o
	$(CC) -o spectrum spectrum.o windows.o -lm -lfftw3

clean:
	rm -f *.o spectrum
