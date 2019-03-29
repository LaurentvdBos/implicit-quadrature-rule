CFLAGS=-O3 -march=native -g
LDLIBS=-lm

.PHONY: all clean

all: main
clean:
	-rm main

main: main.o matrix.o isort.o total_sequence.o trie.o
