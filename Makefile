CFLAGS = -O3 -march=native -g -Wall -pedantic

# If you want to statically link, use this:
LDFLAGS = -static
LDLIBS = -Wl,--start-group -lm -llapacke -llapack -lblas -lgfortran -lquadmath -Wl,--end-group

# Otherwise, simply use this:
LDFLAGS = 
LDLIBS = -lm -llapacke

OBJECTS = implquad.o matrix.o isort.o total_sequence.o trie.o

.PHONY: all clean

all: implquad
clean:
	-rm $(OBJECTS)
	-rm implquad

implquad: $(OBJECTS)
