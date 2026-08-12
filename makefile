CC=g++ -std=c++11
CFLAGS=-g -Wall -O3 -pthread
IFLAGS=-I/opt/X11/include
LFLAGS=-L/opt/X11/lib -lX11
all: wiggle
wiggle: wiggle.cc graphics.cc polynomial.cc roots.cc figures.cc prune.cc certify.cc interface.cc
	$(CC) $(CFLAGS) $(IFLAGS) -o wiggle wiggle.cc $(LFLAGS) -lm

clean:
	rm -f wiggle

