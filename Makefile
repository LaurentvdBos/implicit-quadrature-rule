CFLAGS = -O3 -march=native -g -Wall
LDLIBS = -lm -llapack -llapacke

OBJECTS = implquad.o matrix.o isort.o total_sequence.o trie.o

.PHONY: all clean

all: implquad
clean:
	-rm $(OBJECTS)
	-rm implquad

implquad: $(OBJECTS)
