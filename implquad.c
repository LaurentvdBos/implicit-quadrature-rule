#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>

int d = 0, n = 0;

void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s [-d dim] [-n number of nodes]\n\n", myname);
	fprintf(stderr, "This program start reading samples from standard input until eof and\n"
	                "prints the quadrature consisting on n nodes to standard output.\n");
}

int main(int argc, char **argv)
{
	if (argc < 1) {
		usage("[this-program]");
		return EXIT_FAILURE;
	} else if (argc == 1) {
		usage(argv[0]);
		return EXIT_FAILURE;
	}

	int opt;
	while ((opt = getopt(argc, argv, "d:n:h?")) != -1) {
		switch (opt) {
			case 'd':
				d = atoi(optarg);
				break;
			case 'n':
				n = atoi(optarg);
				break;
			case 'h':
			case '?':
			default:
				usage(argv[0]);
				return EXIT_FAILURE;
		}
	}

	if (d <= 0) {
		fprintf(stderr, "%s: invalid value for 'd': %d\n", argv[0], d);
		usage(argv[0]);
		return EXIT_FAILURE;
	}
	if (n <= 0) {
		fprintf(stderr, "%s: invalid value for 'n': %d\n", argv[0], n);
		usage(argv[0]);
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
