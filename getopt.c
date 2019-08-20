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

// Parse command line options. See 'man 3 getopt' on a Linux distribution or
// google 'getopt C' for extensive information about this function.
//
// The format fmt should be a C-string which contains characters that are
// options. If a character is followed by :, it accepts an option. E.g. "a:b"
// accepts the options -a and -b (or -ab) and a should have an argument. Then
// "-b -a [bork]" or "-ba [bork]" is valid, but "-ab [bork]" is not. Repeatedly
// calling this function with the same argc/argv/fmt iterates through all
// options. It returns the character of the option and can be resetted by
// setting optind to 1.
//
// Parsing stops when a non-option is found, when -- is found, or when there
// are no options to parse anymore, whichever comes first. In that case it
// returns -1.
int getopt(int argc, char **argv, const char *fmt)
{
	if (!scanchar || !*scanchar) {
		// Check if optind is still valid
		if (optind >= argc || argv[optind][0] != '-') {
			opt = '?'; optarg = NULL;
			return -1;
		}
		scanchar = argv[optind];

		// Check if we found -- or -
		if (scanchar[1] == '-' || scanchar[1] == 0) {
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
				opt = *scanchar; optarg = argv[++optind];
				scanchar = NULL;
			} else {
				fprintf(stderr, "%s: option '%c' requires argument\n", argv[0], *scanchar);
				scanchar++;

				opt = '?'; optarg = NULL;
			}
		} else {
			opt = *scanchar; optarg = NULL;
			scanchar++;
		}
	}

	// If scanchar is invalid, go to next option
	if (!scanchar || !*scanchar) {
		optind++;
	}

	return opt;
}
