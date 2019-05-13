CFLAGS = -O3 -march=native -g -Wall -pedantic -Wextra
LDLIBS = -lm -llapacke

# If you want to statically link, use this:
#LDFLAGS = -static
#LDLIBS = -Wl,--start-group -lm -llapacke -llapack -lblas -lgfortran -lquadmath -Wl,--end-group

OBJECTS = implquad.o getopt.o matrix.o isort.o total_sequence.o tree.o

.PHONY: all clean

all: implquad
clean:
	-rm $(OBJECTS)
	-rm implquad

implquad: $(OBJECTS)
