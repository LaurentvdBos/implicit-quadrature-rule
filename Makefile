CFLAGS=-O3 -march=native
LDLIBS=-lm

.PHONY: all clean

all: main
clean:
	-rm main

main: main.o matrix.o isort.o
