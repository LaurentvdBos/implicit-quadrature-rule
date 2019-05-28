CFLAGS = -O3 -march=native -g -Wall -pedantic -Wextra $(shell pkg-config --cflags lapacke)
LDLIBS = $(shell pkg-config --libs lapacke)

# Linking statically (necessary for invokation from MATLAB) can be done as follows:
# LDFLAGS = -static
# LDLIBS = -Wl,--start-group $(shell pkg-config --libs --static lapacke) -lgfortran -lquadmath  -Wl,--end-group

OBJECTS = implquad.o getopt.o matrix.o isort.o total_sequence.o tree.o

.PHONY: all clean

all: implquad
clean:
	-rm $(OBJECTS)
	-rm implquad

implquad: $(OBJECTS)
