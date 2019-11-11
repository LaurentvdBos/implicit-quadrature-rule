CFLAGS = -O3 -g -Wall -pedantic -Wextra
LDLIBS = -lm

OBJECTS = implquad.o getopt.o matrix.o lu.o isort.o total_sequence.o tree.o

.PHONY: all clean

all: implquad
clean:
	-rm $(OBJECTS)
	-rm implquad

implquad: $(OBJECTS)
