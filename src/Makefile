CC=gcc
CFLAGS=-I./include/ -O3 -ffast-math -Wall -fPIC
LIBS=

spectrum: spectrum.o windows.o
	$(CC) $(CFLAGS) -o spectrum spectrum.o windows.o $(LIBS) -lm -lkissfft-float -lgsl -lgslcblas

lib: spectrum.o windows.o
	$(CC) spectrum.o windows.o -shared -o libspectrum.so

clean:
	rm -f *.o *.so spectrum
