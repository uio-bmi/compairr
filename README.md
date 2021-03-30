# vdjsearch

The tool vdjsearch can be used to compare two sets of immune receptor repertoires. Each set can contain many samples and each sample can contain many sequences.

For each of the two repertoires there must an input file of tab-separated values. The fields are:

1. Amino acid sequence (single letter code)
2. Frequency (floating point number)
3. V gene
4. J gene
5. Sample name

Example:

`CASSLALGNEQFF	0.0001	TCRBV07-06	TCRBJ02-01	ABC1`

The tool will find the sequences in the two sets that are similar. The user can specify whether 0, 1 or 2 differences are allowed, and whether indels (insertions or deletions) are allowed or only substitutions. By default just a single substitution is allowed.

The V and J genes must also match.

The similar sequences of each sample in each set are found and counted. Their frequency is taken into account and a matrix is output.

There might be bugs in the code, especially when 2 differences and indels are allowed.

The code is multi-threaded, but perhaps not very efficient on many threads yet.

The following command line options are accepted:

General options:
 -d, --differences INTEGER           number (0-2) of differences accepted (1)
 -i, --indels                        allow insertions or deletions (no)
 -h, --help                          display this help and exit
 -t, --threads INTEGER               number of threads to use (1)
 -v, --version                       display version information and exit

Input/output options:
 -l, --log FILENAME                  log to file, not to stderr
 -o, --output-file FILENAME          output results to file (stdout)
