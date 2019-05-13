#include "getopt.h"

#include <string.h>
#include <stdio.h>

// Single character storing the option value
char opt;

// Value of the option
char *optarg;

// Index of the next option that will be parsed
int optind = 1;

// Index used if multi-char option is found
static char *scanchar = NULL;

int getopt(int argc, char **argv, const char *fmt)
{
	if (!scanchar || !*scanchar) {
		// Check if optind is still valid
		if (optind >= argc || argv[optind][0] != '-') {
			optind++;
			optarg = NULL;
			return -1;
		}
		scanchar = argv[optind];

		// Check if we found --
		if (scanchar[1] == '-') {
			scanchar = NULL;
			optind++;

			opt = '?'; optarg = NULL;
			return -1;
		}

		// Now scanchar points to the first character that must be scanned
		scanchar++;
	}

	// Find the option in fmt
	char *i = strchr(fmt, *scanchar);
	if (!i) {
		fprintf(stderr, "%s: illegal option '%c'\n", argv[0], *scanchar);
		scanchar++;
		
		opt = '?'; optarg = NULL;
	} else {
		// If this option requires an argument, see whether it is provided
		if (i[1] == ':') {
			if (optind < argc-1 && !scanchar[1]) {
				opt = *scanchar; optarg = argv[optind+1];
				optind += 2;
				scanchar = NULL;
			} else {
				fprintf(stderr, "%s: option '%c' requires argument\n", argv[0], *scanchar);
				scanchar++;

				// If scanchar is invalid, go to next option
				if (!scanchar || !*scanchar) {
					optind++;
				}

				opt = '?'; optarg = NULL;
			}
		} else {
			opt = *scanchar; optarg = NULL;
			scanchar++;

			// If scanchar is invalid, go to next option
			if (!scanchar || !*scanchar) {
				optind++;
			}
		}
	}


	return opt;
}
