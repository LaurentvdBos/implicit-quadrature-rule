#ifndef GETOPT_H
#define GETOPT_H

// getopt is a straightforward approach to parse single character command line
// options. It is not available for all platforms, so we provide our own
// implementation.

extern char opt;
extern char *optarg;
extern int optind;

int getopt(int argc, char **argv, const char *fmt);

#endif
