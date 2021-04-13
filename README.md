# vdjsearch

The command line tool `vdjsearch` can be used to compare two sets of
immune receptor repertoires. Each set can contain many repertoires
(~1000) and each repertoire can contain many sequences (~100000).

The tool will find the sequences in the two sets that are similar and
output a matrix with results.

The user can specify whether 0, 1 or 2 differences are allowed when
comparing sequences, using the option `-d` or `--differences`. To allow
indels (insertions or deletions) the option `-i` or `--indels` may be
specified, otherwize only substitutions are allowed.  By default, no
differences are allowed.

The V and J genes specified for each sequence must also match, unless
the `-g` or `--ignore-genes` option is in effect.

The similar sequences of each sample in each set are found.  Their
frequencies is taken into account and a matrix is output containing a
value for each combination of samples in the two sets. The value is
the sum of the products of the frequency of the two sequences that
match in the two samples. If the option `-f` or `--ignore-frequency` is
specified, the frequency information is ignored and all frequencies
set to 1.0.

There is a bug in the code affecting the case when 2 differences and
indels are allowed at the same time. That combination is therefore
currently disabled and it is probably not useful anyway.

The code is multi-threaded. The number of threads may be specified
with the `-t` or `--threads` option.

Run the program with `-v` or `--version` for version information.

Use the `-h` or `--help` option to show some help information.

While the program is running it will print some status and progress
information to standard error (stderr) unless a log file has been
specified with the `-l` or `--log` option.


## Input files

For each of the two repertoires there must an input file of
tab-separated values. The two input files are specified on the command
line without any preceeding option letter.  The files must contain one
line per sequence. The five columns are:

1. Amino acid sequence (single letter code)
2. Frequency (floating point number, exponent allowed)
3. V gene name
4. J gene name
5. Sample name

See below for an example. The program is currently a bit picky on the
input.


## Output file

The output file is also a plain text file with tab-separated values,
containing the output matrix. The first line contains the text
"#sample" followed by the sample names in the second set. The
following lines contains the sample name in the first set, followed by
the values corresponding to the comparison of this sample with each of
the samples in the second set. The output is written to standard out
(stdout) unless a file name has been specified with the `-o` or
`--output-file` option.


## Command line options

```
vdjsearch 0.0.3 - Immune repertoire analysis

Usage: vdjsearch [OPTIONS] TSVFILE1 TSVFILE2

General options:
 -d, --differences INTEGER           number (0-2) of differences accepted (0)
 -i, --indels                        allow insertions or deletions (no)
 -f, --ignore-frequency              ignore frequency / read count information
 -g, --ignore-genes                  ignore V and J gene information
 -h, --help                          display this help and exit
 -t, --threads INTEGER               number of threads to use (1)
 -v, --version                       display version information and exit

Input/output options:
 -l, --log FILENAME                  log to file, not to stderr
 -o, --output-file FILENAME          output results to file (stdout)
```


## Example

Let's use two simple input files. The first is `seta.tsv`:

```
CASSTSHEQYF	0.01	TCRBV07-06	TCRBJ02-01	A1
CASSLRVGGYGYTF	0.03	TCRBV07-09	TCRBJ01-02	A2
```

The second is `setb.tsv`:

```
CASSLRVGGYGYTF	0.05	TCRBV07-09	TCRBJ01-02	B1
CASSLRVGGFGYTF	0.1	TCRBV07-09	TCRBJ01-02	B1
CASSTSHQQYF	0.07	TCRBV07-06	TCRBJ02-01	B2
```

We run the following command:

`vdjsearch -d 1 -o output.tsv seta.tsv setb.tsv`

Here is the output to the console:

```
vdjsearch 0.0.3 - Immune repertoire analysis

Repertoire set 1:  seta.tsv
Repertoire set 2:  setb.tsv
Ignore frequency:  No
Ignore genes:      No
Differences (d):   1
Indels (i):        No
Output file (o):   output.tsv
Threads: (t)       1

Immune receptor repertoire set 1
Reading sequences: 100%
Sequences: 2, residues: 25, shortest: 11, longest: 14, average: 12.5
Samples:           2
Indexing:          100%
Sorting:           100%

#no	seqs	freq	sample
1	      1	0.01000	A1
2	      1	0.03000	A2
Sum	2	0.04000

Immune receptor repertoire set 2
Reading sequences: 100%
Sequences: 3, residues: 39, shortest: 11, longest: 14, average: 13.0
Samples:           2
Indexing:          100%
Sorting:           100%

#no	seqs	freq	sample
1	      2	0.15000	B1
2	      1	0.07000	B2
Sum	3	0.22000

Unique v_genes:    2
Unique d_genes:    2
Computing hashes:  100%
Computing hashes:  100%
Hashing sequences: 100%
Analysing:         100%

```

The program gives some statistics on the input files after reading
them.

Here is the result in the `output.tsv` file:

```
#sample	B1	B2
A1	0.0e+00	7.0e-04
A2	4.5e-03	0.0e+00
```

Here, the sequence in sample A1 is similar to the sequence in sample
B2. The only difference is the E and Q in the 8th position. The gene
names are also the same. They have frequencies of 0.01 and 0.07,
respectively. The product is 0.0007 or 7.0e-04. That value is found in
the third column on the second line in the output file.

The sequence in sample A2 with frequency 0.03 is similar to both
sequences in sample B1, with frequencies 0.1 and 0.07. The first
sequence in B1 is identical, while the second sequence in B1 has an F
instead of a Y in the 10th position. The result is 0.03 * (0.05 + 0.1)
= 0.03 * 0.15 = 0.0045 = 4.5e-03. That values is found in the second
column on the third line.

Since there are no sequences from sample A1 similar to B1 or from A2
similar to B1, the other values are zero.


## Implementation

The program is written in C++. The strategy for finding similar
sequences is based on a similar concept developed in
[swarm](https://github.com/torognes/swarm). Basically, a 64-bit hash
is computed for all sequences in the sets. All hashes for one set are
stored in a Bloom filter and in a hash table. We then look for matches
to sequences in the second set by looking them up in the Bloom filter
and then, if there was a match, in the hash table. To find matches
with 1 or 2 substitutions or indels, the hashes of all these
'microvariants' are generated and looked up.


## Binaries

The code should compile easily by running `make` in the `src`
directory. Binaries for Linux and macOS are distributed with each
release.


## Performance

As a preliminary performance test, Cohort 2 ("Keck") of
[the dataset](https://s3-us-west-2.amazonaws.com/publishedproject-supplements/emerson-2017-natgen/emerson-2017-natgen.zip)
by Emerson et al. was compared to itself. It contains 120 repertoires
with a total of 24 544 336 sequences. The timing results are shown below.

Distance | Indels | Threads | Time (s) | Time (min) | Note
-------: | :----: | ------: | -------: | ---------: | ----
0 | no  | 1 |   74 |  1.14 | Almost all time was used to read files
1 | no  | 1 |  159 |  2.39
1 | no  | 4 |   97 |  1.37
1 | yes | 1 |  236 |  3.56
1 | yes | 4 |  123 |  2.03
2 | no  | 4 | 2871 | 48.51

Memory usage was 3.44GB, corresponding to an average of about 70 bytes
per sequence. Times are wall time measured by `/usr/bin/time`.

The analysis was performed on a Mac Mini M1.


## References

* Emerson RO, DeWitt WS, Vignali M, Gravley J, Hu JK, Osborne EJ, Desmarais C, Klinger M, Carlson CS, Hansen JA, Rieder M, Robins HS (2017)
**Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire.**
*Nature Genetics*, 49 (5): 659-665.
doi:[10.1038/ng.3822](https://doi.org/10.1038/ng.3822)


## Future plans

Some potential features:

* Given a set of sequences and a repertoire, produce a table of which
sequences are found in which sample of the repertoire.
