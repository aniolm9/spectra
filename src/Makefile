CC=gcc
CFLAGS=-I./include/
LIBS=

spectrum: spectrum.o windows.o
	$(CC) $(CFLAGS) -o spectrum spectrum.o windows.o $(LIBS) -lm -lkissfft-float

clean:
	rm -f *.o spectrum
