# CompAIRR

CompAIRR (`compairr`) is a command line tool to compare two sets of
adaptive immune receptor repertoires, or cluster the sequences in one
repertoire.


## Installation

The code is C++11 standard compliant and should compile easily using a
modern C++ compiler (e.g. GNU GCC or LLVM Clang). Run `make clean`,
`make`, `make test` or `make install` in the main folder to clean,
build, test or install the tool. It is no dependencies except for the
C standard library and the C++ standard libraries.

Binaries for Linux (x86_64) and macOS (x86_64 and ARM64) are also
distributed with each
[release](https://github.com/uio-bmi/compairr/releases/latest).

A `Dockerfile` is included if you want to make a Docker image.


## General options

Use the `-h` or `--help` option to show some help information.

Run the program with `-v` or `--version` for version information.

The code is multi-threaded. The number of threads may be specified
with the `-t` or `--threads` option.

The results will be written to standard out (stdout) unless a file
name has been specified with the `-o` or `--output-file` option.

While the program is running it will print some status and progress
information to standard error (stderr) unless a log file has been
specified with the `-l` or `--log` option. Error messages and warnings
will also be written here.

The user can specify whether 0, 1 or 2 differences are allowed when
comparing sequences, using the option `-d` or `--differences`. To allow
indels (insertions or deletions) the option `-i` or `--indels` may be
specified, otherwise only substitutions are allowed. By default, no
differences are allowed. The `-i` option is allowed only when d=1.

The V and J genes specified for each sequence must also match, unless
the `-g` or `--ignore-genes` option is in effect.


## Computing overlap between two repertoire sets

To compute the overlap between two repertoire sets, use the `-m` or
`--matrix` option.

For each of the two repertoire sets there must an input file of
tab-separated values. The two input files are specified on the command
line without any preceeding option letter. If only one filename is
specified on the command line, or the same filename is specified
twice, it is assumed that the set should be compared to itself. Each
file must contain the amino acid sequence of the rearrangement, the V
gene, the J gene, the duplicate count and the repertoire ID.

Each set can contain many repertoires and each repertoire can
contain many sequences.

The tool will find the sequences in the two sets that are similar and
output a matrix with results.

The similar sequences of each repertoire in each set are found by
comparing the sequences and their V and J genes.  The duplicate count
of each sequence is taken into account and a matrix is output
containing a value for each combination of repertoires in the two
sets. The value is the sum of the products of the counts of all pairs of
sequences in the two repertoires that match. If the option `-f` or
`--ignore-frequency` is specified, the duplicate count information is
ignored and all counts are set to 1.

The output will be a matrix of values in a tab-separated plain text
file. Two different formats can be selected. In the default format,
the first line contains the character '#' followed by the repertoire
ID's from the second set. The following lines contains the repertoire
ID from the first set, followed by the values corresponding to the
comparison of this repertoire with each of the repertoires in the
second set.

An alternative output format is used when the `-a` or `--alternative`
option is specified. It will write the results in a three column
format with the repertoire ID from set 1 and set 2 in the two first
columns, respectively, and the value in the third column. There will
be one line for each combination of repertoires in the sets.


## Clustering the sequences in a repertoire

To cluster the sequences in one repertoire, use the `-c` or
`--cluster` option.

One input file in the tab-separated format must be specified on the
command line.

The tool will cluster the sequences using single linkage hierarchical
clustering, according to the specified distance and indel options
(`-d`, `--distance`, `-i`, `--indels`). The V and J genes will be
taken into account unless the `-g` or `--ignore-genes` option is
specified.

The output will be in a similar format as the input file, but preceeded
with two additional columns. The first column will contain a cluster
number, starting at 1. The second column will contain the size of the
cluster. The clusters are sorted by size, in descending order.


## Input files

The input files must be in tab-separated value (TSV) format accoring
to the
[Rearrangement Schema](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)
of the
[AIRR standards 1.3 documentation](https://docs.airr-community.org/en/stable/).

The first line must contain the header. The rest of the file must
contain one line per sequence. The following fields are required:

* junction_aa: amino acid sequence (single letter code)
* duplicate_count: number of identical copies of the same rearrangement
* v_call: V gene name with allele
* j_call: J gene name with allele
* repertoire_id: identifier of the repertoire

See below for an example. Other fields may be included, but will be
ignored.


## Command line options

```
CompAIRR 1.0.0 - Compare Adaptive Immune Receptor Repertoires
https://github.com/uio-bmi/compairr

Usage: compairr [OPTIONS] TSVFILE1 [TSVFILE2]

Commands:
 -h, --help                  display this help and exit
 -v, --version               display version information
 -m, --matrix                compute overlap matrix between two sets
 -c, --cluster               cluster sequences in one repertoire

General options:
 -d, --differences INTEGER   number (0-2) of differences accepted (0)
 -i, --indels                allow insertions or deletions
 -f, --ignore-frequency      ignore duplicate count information
 -g, --ignore-genes          ignore V and J gene information
 -t, --threads INTEGER       number of threads to use (1)

Input/output options:
 -a, --alternative           output overlap results in column format
 -l, --log FILENAME          log to file (stderr)
 -o, --output FILENAME       output results to file (stdout)
```


## Example - Repertoire overlap

Let's use two simple input files. The first is `seta.tsv`:

```tsv
junction_aa	duplicate_count	v_call	j_call	repertoire_id
CASSTSHEQYF	1	TCRBV07-06	TCRBJ02-01	A1
CASSLRVGGYGYTF	3	TCRBV07-09	TCRBJ01-02	A2
```

The second is `setb.tsv`:

```tsv
junction_aa	duplicate_count	v_call	j_call	repertoire_id
CASSLRVGGYGYTF	5	TCRBV07-09	TCRBJ01-02	B1
CASSLRVGGFGYTF	10	TCRBV07-09	TCRBJ01-02	B1
CASSTSHQQYF	7	TCRBV07-06	TCRBJ02-01	B2
```

We run the following command:

`compairr -m -d 1 -o output.tsv seta.tsv setb.tsv`

Here is the output to the console:

```
CompAIRR 1.0.0 - Immune repertoire analysis
https://github.com/uio-bmi/compairr

Command:           Overlap
Repertoire set 1:  seta.tsv
Repertoire set 2:  setb.tsv
Differences (d):   1
Indels (i):        No
Ignore counts (f): No
Ignore genes (g):  No
Output file (o):   output.tsv
Output format (a): Matrix
Threads: (t)       1

Immune receptor repertoire set 1

Reading sequences: 100%
Sequences:         2
Residues:          25
Shortest:          11
Longest:           14
Average length:    12.5
Repertoires:       2
Indexing:          100%
Sorting:           100%

Repertoires:
  # Sequences Count Repertoire ID
  1         1     1 A1
  2         1     3 A2
Sum         2     4

Immune receptor repertoire set 2

Reading sequences: 100%
Sequences:         3
Residues:          39
Shortest:          11
Longest:           14
Average length:    13.0
Repertoires:       2
Indexing:          100%
Sorting:           100%

Repertoires:
  # Sequences Count Repertoire ID
  1         2    15 B1
  2         1     7 B2
Sum         3    22

Unique V genes:    2
Unique J genes:    2
Computing hashes:  100%
Computing hashes:  100%
Hashing sequences: 100%
Analysing:         100%
Writing results:   100%
```

The program gives some statistics on the input files after reading
them.

Here is the result in the `output.tsv` file:

```tsv
#	B1	B2
A1	0	7
A2	45	0
```

Here, the sequence in repertoire A1 is similar to the sequence in
repertoire B2. The only difference is the E and Q in the 8th
position. The gene names are also the same. They have duplicate counts
of 1 and 7, respectively. The product is 7. That value is found in the
third column on the second line in the output file.

The sequence in repertoire A2 with frequency 3 is similar to both
sequences in repertoire B1, with frequencies 10 and 5. The first
sequence in B1 is identical, while the second sequence in B1 has an F
instead of a Y in the 10th position. The result is 3 * (5 + 10)
= 3 * 15 = 45. That values is found in the second
column on the third line.

Since there are no sequences from repertoire A1 similar to B1 or from
A2 similar to B1, the other values are zero.

This small dataset is included in the test folder and the tool can
automatically be tested by running `make test`.


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


## Performance

As a preliminary performance test, Cohort 2 ("Keck") of
[the dataset](https://s3-us-west-2.amazonaws.com/publishedproject-supplements/emerson-2017-natgen/emerson-2017-natgen.zip)
by Emerson et al. was compared to itself. It contains 120 repertoires
with a total of 24 205 557 extracted sequences. The test was performed
with CompAIRR version 1.0.0. The timing results are shown below.

Distance | Indels | Threads | Time (s) | Time (mm:ss)
-------: | :----: | ------: | -------: | -----------:
0 | no | 1 | 37 | 0:37
0 | no | 4 | 32 | 0:32
1 | no | 1 | 244 | 4:04
1 | no | 4 | 89 | 1:29
1 | yes | 1 | 402 | 6:42
1 | yes | 4 | 133 | 2:13
2 | no | 4 | 3437 | 57:17

When the distance is zero almost all of the time was used to read
files.

Memory usage was 2.5GB, corresponding to an average of about 100 bytes
per sequence.

Since this is a comparison of a repertoire set to itself, the dataset
is only read once, and the memory needed is also reduced as compared
to a situation were two different repertoire sets are compared.

Wall time and memory usage was measured by `/usr/bin/time`. The
analysis was performed on a Mac Mini M1.


## Development team

The code has been developed by Torbjørn Rognes based on code from
Swarm where Frédéric Mahé and Lucas Czeck made important
contributions. Geir Kjetil Sandve had the idea of developing a tool
for rapid repertoire set comparison. Lonneke Scheffer has tested and
benchmarked the tool, and suggested new features. Milena Pavlovic and
Victor Greiff have also contributed to the project.


## Citing CompAIRR

We are preparing a manuscript about CompAIRR, but it is not yet
available. For the time being, please cite the following if you use
CompAIRR in any published work:

* Rognes T, Scheffer L, Greiff V, Sandve GK (2021)
**CompAIRR: Efficient computation of adaptive immune receptor repertoire-similarity allowing non-exact sequence matching.**
(in prep.)


## Support

We will prioritize fixing important bugs. We will also try to answer
questions, improve documentation and implement suggested enhancements
as time permits. As we have no dedicated funding for this project we
cannot make any guarantees on the level of support.

To report a potential bug, suggest enhancements or ask questions,
please use one of the following means:

* [Submit an issue on GitHub](https://github.com/uio-bmi/compairr/issues) (preferred)

* Send an email to [`torognes@ifi.uio.no`](mailto:torognes@ifi.uio.no)

If you would like to contribute with code you are most welcome to
[submit a pull request](https://github.com/uio-bmi/compairr/pulls).


## References

* Emerson RO, DeWitt WS, Vignali M, Gravley J, Hu JK, Osborne EJ, Desmarais C, Klinger M, Carlson CS, Hansen JA, Rieder M, Robins HS (2017)
**Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire.**
*Nature Genetics*, 49 (5): 659-665.
doi:[10.1038/ng.3822](https://doi.org/10.1038/ng.3822)
