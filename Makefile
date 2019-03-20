CFLAGS=-O3
LDLIBS=-lm

.PHONY: all clean

all: main
clean:
	-rm main

main: main.o matrix.o
