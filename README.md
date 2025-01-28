[![](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)

# CompAIRR

CompAIRR (`compairr`) is a command line tool to compare two sets of
adaptive immune receptor repertoires and compute their overlap. It can
also identify which sequences are present in which repertoires.
Furthermore, CompAIRR can cluster the sequences in a repertoire
set. Sequence comparisons can be exact or approximate. CompAIRR has
been shown to be very fast and to have a small memory footprint
compared to similar tools, when up to 2 differences are allowed.


## Installation

The code is C++11 standard compliant and should compile easily using
`make` and a modern C++ compiler (e.g. GNU GCC or LLVM Clang). Run
`make clean`, `make`, `make test` and `make install` in the main
folder to clean, build, test and install the tool. There are no
dependencies except for the C and C++ standard libraries.

Binaries for Linux (x86_64) and macOS (x86_64 and Arm64) are also
distributed with each
[release](https://github.com/uio-bmi/compairr/releases/latest).

A `Dockerfile` is included if you want to make a Docker image.  A
docker image may be built with the following command:

```sh
docker build -t compairr .
```

Ready-made Docker images for CompAIRR can be found on the
[Docker Hub](https://hub.docker.com/r/torognes/compairr).

CompAIRR can be installed on macOS using homebrew with
`brew install torognes/bioinf/compairr`.


## Tutorial

For an introduction to how to use CompAIRR, please have a look at the
[CompAIRR tutorial](https://github.com/LonnekeScheffer/compairr-tutorial).


## General options

Use the `-h` or `--help` option to show some help information.

Run the program with `-v` or `--version` for version information.

The type of operation that should be performed is specified with one
of the options `-m`, `-x`, `-c` or `-z` (or the corresponding long option
forms `--matrix`, `--existence`, `--cluster`, or `--deduplicate`).

The code is multi-threaded. The number of threads may be specified
with the `-t` or `--threads` option.

The results will be written to standard out (stdout) unless a file
name has been specified with the `-o` or `--output-file` option.

While the program is running it will print some status and progress
information to standard error (stderr) unless a log file has been
specified with the `-l` or `--log` option. Error messages and warnings
will also be written here.

The default is to compare amino acid sequences, but nucleotide
sequences are compared if the `-n` or `--nucleotides` option is given.
The accepted amino acid symbols are `ACDEFGHIKLMNPQRSTVWY`, while the
accepted nucleotide symbols are `ACGTU`. Lower case letters are also
accepted. The program will abort with an error message if any other
symbol is encountered in a sequence, unless one specifies the `-u` or
`--ignore-unknown` option, in which case CompAIRR will simply ignore
that sequence. If the program encounters an empty sequence it will
also abort with an error message, unless the `-e` or `--ignore-empty`
option is given.

By default, the sequences should be given in the `junction` or
`junction_aa` column of the input file, for nucleotide and amino acid
sequences, respectively. Alternatively, the sequences may be present
in the `cdr3` or `cdr3_aa` column, if the `--cdr3` option is given.

The user can specify how many differences are allowed when comparing
sequences, using the option `-d` or `--differences`. To allow indels
(insertions or deletions) the option `-i` or `--indels` may be
specified, otherwise only substitutions are allowed. By default, no
differences are allowed. The `-i` option is allowed only when d=1. The
number of differences allowed strongly influences the speed of
CompAIRR. The program will be slower as more differences
are allowed. When d=0 or d=1 it is very fast, but it will be relatively
slow with d=2 and even slower when d>2. See the section on performance
below for an example.

The V and J gene alleles specified for each sequence must also match,
unless the `-g` or `--ignore-genes` option is in effect.


## Computing overlap between two repertoire sets

To compute the overlap between two repertoire sets, use the `-m` or
`--matrix` option.

For each of the two repertoire sets there must an input file of
tab-separated values formatted according to [the AIRR standard for
rearrangements](https://docs.airr-community.org/en/stable/datarep/rearrangements.html).
The two input files are specified on the command line without any
preceding option letter. If only one filename is specified on the
command line, or the same filename is specified twice, it is assumed
that the set should be compared to itself. Each file must contain the
repertoire ID and either the nucleotide or the amino acid sequence of
the rearrangement. If the repertoire ID column is missing, all
sequences are assumed to belong to the same repertoire (with ID 1 or
2, respectively, for the two sets). A sequence ID may also be
included. Unless they should be ignored, the V gene, the J gene, and
the duplicate count is also needed.

Each set can contain many repertoires and each repertoire can contain
many sequences. The tool will find the sequences in the two sets that
are similar and output a matrix with results.

CompAIRR assumes that all sequences within each repertoire are
distinct, and that the abundance of each sequence is indicated in the
`duplicate_count` field in the input file. Duplicated sequences,
i.e. identical sequences (with the same V and J genes) within the same
repertoire, may lead to unexpected results. CompAIRR will warn if it
detects duplicates. Duplicates may be merged with the `--deduplicate`
command.

The similar sequences of each repertoire in each set are found by
comparing the sequences and their V and J genes.  The duplicate count
of each sequence is taken into account and a matrix is output
containing a value for each combination of repertoires in the two
sets. The value is usually the sum of the products of the duplicate
counts of all pairs of sequences in the two repertoires that match. If
the option `-f` or `--ignore-counts` is specified, the duplicate count
information is ignored and all counts are treated as 1. Instead of
summing the product of the counts, the ratio, min, max, or mean may be
used if specified with the `-s` or `--score` option. The Morisita-Horn
index or Jaccard index will be calculated if `MH` or `Jaccard` is
specified with the `-s` option. These indices can only be computed
when d=0.

The output will be a matrix of values in a tab-separated plain text
file. Two different formats can be selected. In the default format,
the first line contains the hash character (`#`) followed by the
repertoire ID's from the second set. The following lines contains the
repertoire ID from the first set, followed by the values corresponding
to the comparison of this repertoire with each of the repertoires in
the second set.

An alternative output format is used when the `-a` or `--alternative`
option is specified. It will write the results in a three column
format with the repertoire ID from set 1 and set 2 in the two first
columns, respectively, and the value in the third column. There will
be one line for each combination of repertoires in the sets. The very
first line will contain a hash character (`#`) followed by the field
names separated by tabs.

If the `-p` or `--pairs` option is specified, CompAIRR will write
information about all pairs of matching sequences to a specified TSV
file. Please note that such files may grow very large when there are
many matches. Use of multithreading may be of little use in this
case. The order of the lines in the file is unspecified. The following
columns from both input files will be included in the output:
`repertoire_id`, `sequence_id`, `duplicate_count`, `v_call`, `j_call`,
and `junction`. The term `junction` will be replaced with
`junction_aa`, `cdr3`, or `cdr3_aa` as appropriate. Additional columns
from the input files may be copied to the pairs file using the `-k` or
`--keep-columns` option. Multiple columns, separated by commas (but no
spaces), may be given. A warning will be given if any of the specified
columns are missing. In the header, columns from the first and second
input file will be suffixed by `_1` and `_2`, respectively. The
distance between the sequences will be included if the `--distance`
option is included. This is usually the Hamming distance (minimum
number of substitutions), unless the `--indel` (or `-i´) option is
specified, in which case the distance is the Levenshtein distance
(minimum number of substitutions or indels).


## Analysing in which repertoires a set of sequences are present

Use the option `-x` or `--existence` to analyse in which repertoires a
set of sequences are present, and create a sequence presence matrix.

Two input files with repertoire sets in standard format must be
specified on the command line. The first file should contain the
different sequences to analyse. The `sequence_id` column must be
present in this file. If the optional `repertoire_id` column is
present, all those identifiers must be identical. The second file must
contain the repertoires to match. The `repertoire_id` column must be
present in the second file, otherwise the ID will be set to 2 for all
sequences.

CompAIRR will identify in which repertoires each sequence is present
and will output the results either as a matrix or as a three-column
table (if the `-a` option is specified). The options `-d`, `-i`, `-g`,
and `-n` (and the corresponding long option names `--differences`,
`--indels`, `--ignore-genes`, and `--nucleotides`) will be taken into
account when comparing sequences.

The output will be in a similar format as when computing the overlap
(above), but the first column will contain the `sequence_id` from the
first file instead of the `repertoire_id`.

The `-p` or `--pairs` option may be specified to output all pairs of
matching sequences in the same way as for the overlap computation.


## Clustering the sequences in a repertoire

To cluster the sequences in one repertoire, use the `-c` or
`--cluster` option.

One input file in tab-separated format must be specified on the
command line.

The tool will cluster the sequences using single linkage hierarchical
clustering, according to the specified distance and indel options
(`-d`, `--distance`, `-i`, `--indels`). The V and J gene alleles will
be taken into account unless the `-g` or `--ignore-genes` option is
specified. The options `-n` or `--nucleotides` indicate that the
comparison should be performed with nucleotide sequences, not amino
acid sequences. If the repertoire ID column is missing, all
sequences are assumed to belong to the same repertoire (with ID 1).

The output will be in a similar TSV format as the input file, but
preceded with two additional columns. The first column will contain a
cluster number, starting at 1. The second column will contain the size
of the cluster. The subsequent columns are `repertoire_id`,
`sequence_id`, `duplicate_count`, `v_call`, `j_call`, and `junction`
(or `junction_aa`, `cdr3` or `cdr3_aa`, as appropriate).

The clusters are sorted by size, in descending order.


## Deduplication

The `--deduplicate` command may be used to deduplicate a data set by
merging entries in the same repertoire with identical sequences and
identical V and J genes. This may be necessary to get correct results
when computing overlaps between repertoires. Duplicates may be present
for instance in cases were the data set contains both nucleotide and
amino acid sequences from the same rearrangement, where the nucleotide
sequences may be distinct while the amino acid sequences may not be,
due to the degeneracy of the genetic code.

One input file in TSV format must be specified on the command line.

Strictly identical sequences in the same repertoire will be merged and
their counts will be added together. If the `-g` or `--ignore_genes`
option is specified, the V and J genes are ignored. The `-n` or
`--nucleotides` option may be specified if the input is nucleotide
sequences, otherwise amino acid sequences will be assumed. If the `-f`
or `--ignore_counts` option is specified, the counts in the input file
will be ignored, and just the number of identical sequences will be
counted. If the repertoire ID column is missing, all sequences are
assumed to belong to the same repertoire (with ID 1).

The output will be in a similar TSV format as the input file, with the
following columns: `repertoire_id`, `duplicate_count`, `v_call`,
`j_call`, and `junction` (or `junction_aa`, `cdr3` or `cdr3_aa`, as
appropriate). If the `-g` or `--ignore_genes` option is specified, the
`v_call` and `j_call` columns will not be included.


## Input files

The input files must be in tab-separated value (TSV) format accoring
to the [Rearrangement
Schema](https://docs.airr-community.org/en/stable/datarep/rearrangements.html)
of the [AIRR standards 1.3
documentation](https://docs.airr-community.org/en/stable/).

The first line must contain the header. The rest of the file must
contain one line per sequence. The following fields should be included:

* `repertoire_id`: identifier of the repertoire
* `sequence_id`: identifier of the sequence (optional except for for first file when using `-x` or `--existence`)
* `duplicate_count`: number of identical copies of the same rearrangement (required unless `-f` option given)
* `v_call`: V gene name with allele (required unless `-g` option given)
* `j_call`: J gene name with allele (required unless `-g` option given)
* `junction`: nucleotide sequence (required if `-n` option given and `--cdr3` option not given)
* `junction_aa`: amino acid sequence (single letter code) (required unless `-n` or `--cdr3` options given)
* `cdr3`: nucleotide sequence (required if both `-n` and `--cdr3` options given)
* `cdr_aa`: amino acid sequence (single letter code) (required if `--cdr3` option given and `-n` option not given)

See below for an example. Other fields may be included, but will be
ignored.


## Command line option overview

The command line should look like this:

```
compairr OPTIONS TSVFILE1 [TSVFILE2]
```

Exactly one of the command options `-m`, `-x` or `-c` (or their long forms) must be specified. Other options as indicated in the table below could also be included. With the `-m` and `-x` command options, the names of two tab-separated value files with repertoires must also be specified on the command line, with the `-c` command option, only one such file should be specified.

Short | Long               | Argument | Default  | Description
------|--------------------|----------|----------|-------------
`-a`  | `--alternative`    |          |          | Output results in three-column format, not matrix
`  `  | `--cdr3`           |          |          | Use the `cdr3` or `cdr3_aa` column instead of `junction` or `junction_aa`
`-c`  | `--cluster`        |          |          | Cluster sequences in one repertoire
`-d`  | `--differences`    | INTEGER  | 0        | Number of differences accepted
`  `  | `--distance`       |          |          | Include sequence distance in pairs file
`-e`  | `--ignore-empty`   |          |          | Ignore empty sequences
`-f`  | `--ignore-counts`  |          |          | Ignore duplicate count information
`-g`  | `--ignore-genes`   |          |          | Ignore V and J gene information
`-h`  | `--help`           |          |          | Display help text and exit
`-i`  | `--indels`         |          |          | Allow insertions or deletions
`-k`  | `--keep-columns`   | STRING   |          | Copy given comma-separated columns to pairs file
`-l`  | `--log`            | FILENAME | (stderr) | Log to specified file instead of stderr
`-m`  | `--matrix`         |          |          | Compute overlap matrix between two sets
`-n`  | `--nucleotides`    |          |          | Compare nucleotides, not amino acids
`-o`  | `--output`         | FILENAME | (stdout) | Output results to specified file instead of stdout
`-p`  | `--pairs`          | FILENAME | (none)   | Output matching pairs to specified file
`-s`  | `--score`          | STRING   | product  | Sum `product`, `ratio`, `min`, `max`, or `mean`; or compute `MH` or `Jaccard` index
`-t`  | `--threads`        | INTEGER  | 1        | Number of threads to use (1-256)
`-u`  | `--ignore-unknown` |          |          | Ignore sequences including unknown residue symbols
`-v`  | `--version`        |          |          | Display version information
`-x`  | `--existence`      |          |          | Check existence of sequences in repertoires
`-z`  | `--deduplicate`    |          |          | Deduplicate sequences


## Example 1: Repertoire overlap

In this example we will compute the overlap of two repertoire sets.

Let's use two simple input files. The first is `seta.tsv`:

```tsv
repertoire_id	sequence_id	duplicate_count	v_call	j_call	junction	junction_aa	sequence	rev_comp	productive	d_call	sequence_alignment	germline_alignment	v_cigar	d_cigar	j_cigar
A1	R	1	TCRBV07-06	TCRBJ02-01	tgcgcgagcagcaccagccatgaacagtatttt	CASSTSHEQYF									
A2	S	3	TCRBV07-09	TCRBJ01-02	tgcgcgagcagcctgcgcgtgggcggctatggctataccttt	CASSLRVGGYGYTF									
```


The second is `setb.tsv`:

```tsv
repertoire_id	sequence_id	duplicate_count	v_call	j_call	junction	junction_aa	sequence	rev_comp	productive	d_call	sequence_alignment	germline_alignment	v_cigar	d_cigar	j_cigar
B1	T	5	TCRBV07-09	TCRBJ01-02	tgcgcgagcagcctgcgcgtgggcggctatggctataccttt	CASSLRVGGYGYTF									
B1	U	10	TCRBV07-09	TCRBJ01-02	tgcgcgagcagcctgcgcgtgggcggctttggctataccttt	CASSLRVGGFGYTF									
B2	V	7	TCRBV07-06	TCRBJ02-01	tgcgcgagcagcaccagccatcagcagtatttt	CASSTSHQQYF									
```

We run the following command:

`compairr -m seta.tsv setb.tsv -d 1 -o output.tsv -p pairs.tsv`

Here is the output to the console:

```
CompAIRR 1.7.0 - Comparison of Adaptive Immune Receptor Repertoires
https://github.com/uio-bmi/compairr

Start time:        Thu Mar 03 12:29:32 CET 2022
Command (m/c/x):   Overlap (-m)
Repertoire set 1:  seta.tsv
Repertoire set 2:  setb.tsv
Nucleotides (n):   No
Differences (d):   1
Indels (i):        No
Ignore counts (f): No
Ignore genes (g):  No
Ign. unknown (u):  No
Threads (t):       1
Output file (o):   output.tsv
Output format (a): Matrix
Score (s):         Sum of products of counts
Pairs file (p):    pairs.tsv
Log file (l):      (stderr)

Immune receptor repertoire set 1

Reading sequences: 100% (0s)
Repertoires:       2
Sequences:         2
Residues:          25
Shortest:          11
Longest:           14
Average length:    12.5
Total dupl. count: 4
Indexing:          100% (0s)

Repertoires in set:
# Sequences Count Repertoire ID
1         1     1 A1
2         1     3 A2

Immune receptor repertoire set 2

Reading sequences: 100% (0s)
Repertoires:       2
Sequences:         3
Residues:          39
Shortest:          11
Longest:           14
Average length:    13.0
Total dupl. count: 22
Indexing:          100% (0s)

Repertoires in set:
# Sequences Count Repertoire ID
1         2    15 B1
2         1     7 B2

Unique V genes:    2
Unique J genes:    2
Computing hashes:  100% (0s)
Computing hashes:  100% (0s)
Hashing sequences: 100% (0s)
Analysing:         100% (0s)
Writing results:   100% (0s)

End time:          Thu Mar 03 12:29:32 CET 2022
```

Repertoires will be sorted alphabetically by ID. The program gives some
statistics on the input files after reading them.

Here is the result in the `output.tsv` file:

```tsv
#	B1	B2
A1	0	7
A2	45	0
```

And here is the result in the `pairs.tsv` file:

```tsv
#repertoire_id_1	sequence_id_1	duplicate_count_1	v_call_1	j_call_1	junction_aa_1	repertoire_id_2	sequence_id_2	duplicate_count_2	v_call_2	j_call_2	junction_aa_2
A1	R	1	TCRBV07-06	TCRBJ02-01	CASSTSHEQYF	B2	V	7	TCRBV07-06	TCRBJ02-01	CASSTSHQQYF
A2	S	3	TCRBV07-09	TCRBJ01-02	CASSLRVGGYGYTF	B1	T	5	TCRBV07-09	TCRBJ01-02	CASSLRVGGYGYTF
A2	S	3	TCRBV07-09	TCRBJ01-02	CASSLRVGGYGYTF	B1	U	10	TCRBV07-09	TCRBJ01-02	CASSLRVGGFGYTF
```

Here, sequence R in repertoire A1 is similar to sequence V in
repertoire B2. The only difference is the E and Q in the 8th
position. The gene allele names are also the same. They have duplicate
counts of 1 and 7, respectively. The product is 7. That value is found
in the third column on the second line in the main output file.

Sequence S in repertoire A2 with duplicate count 3 is similar to both
sequence T and U in repertoire B1, with duplicate counts of 5 and 10,
respectively. Sequence T in B1 is identical, while sequence U in B1
has an F instead of a Y in the 10th position. The result is 3 * (5 +
10) = 3 * 15 = 45. That value is found in the second column on the
third line of the main output file.

Since there are no sequences from repertoire A1 similar to B1 or from
A2 similar to B1, the other values in the matrix are zero.

This small dataset is included in the test folder and the tool can
automatically be tested by running `make test`.


## Example 2: Sequence existence

In this example we will use the `-x` or `--existence` option to find
out in which repertoires a set of sequences are present.

The file `setc.tsv` contains the sequences that we will analyse:

```tsv
repertoire_id	sequence_id	duplicate_count	v_call	j_call	junction	junction_aa	sequence	rev_comp	productive	d_call	sequence_alignment	germline_alignment	v_cigar	d_cigar	j_cigar
C	X	1	TCRBV07-09	TCRBJ01-02	tgcgcgagcagcctgcgcgtgggcggctttggctataccttt	CASSLRVGGFGYTF									
C	Y	1	TCRBV07-06	TCRBJ02-01	tgcgcgagcagcaccagccatcagcagtatttt	CASSTSHQQYF									
```

The file above is included in the folder `test` in the distribution.

We will compare it to repertoire sets in the file `setb.tsv` described
earlier.

We run the following command:

`compairr -x setc.tsv setb.tsv -d 1 -f -o output.tsv -p pairs.tsv`

Here is the output to the console:

```
CompAIRR 1.7.0 - Comparison of Adaptive Immune Receptor Repertoires
https://github.com/uio-bmi/compairr

Start time:        Thu Mar 03 12:31:16 CET 2022
Command (m/c/x):   Existence (-x)
Repertoire:        setc.tsv
Repertoire set:    setb.tsv
Nucleotides (n):   No
Differences (d):   1
Indels (i):        No
Ignore counts (f): Yes
Ignore genes (g):  No
Ign. unknown (u):  No
Threads (t):       1
Output file (o):   output.tsv
Output format (a): Matrix
Score (s):         Sum of products of counts
Pairs file (p):    pairs.tsv
Log file (l):      (stderr)

Immune receptor repertoire set 1

Reading sequences: 100% (0s)
Repertoires:       1
Sequences:         2
Residues:          25
Shortest:          11
Longest:           14
Average length:    12.5
Total dupl. count: 2
Indexing:          100% (0s)

Repertoires in set:
# Sequences Count Repertoire ID
1         2     2 C

Immune receptor repertoire set 2

Reading sequences: 100% (0s)
Repertoires:       2
Sequences:         3
Residues:          39
Shortest:          11
Longest:           14
Average length:    13.0
Total dupl. count: 22
Indexing:          100% (0s)

Repertoires in set:
# Sequences Count Repertoire ID
1         2    15 B1
2         1     7 B2

Unique V genes:    2
Unique J genes:    2
Computing hashes:  100% (0s)
Computing hashes:  100% (0s)
Hashing sequences: 100% (0s)
Analysing:         100% (0s)
Writing results:   100% (0s)

End time:          Thu Mar 03 12:31:16 CET 2022
```

Here is the result in the `output.tsv` file:

```tsv
#	B1	B2
X	2	0
Y	0	1
```

Please note that the `-f` option was used to ignore the duplicate
counts.

And here is the result in the `pairs.tsv` file:

```tsv
#repertoire_id_1	sequence_id_1	duplicate_count_1	v_call_1	j_call_1	junction_aa_1	repertoire_id_2	sequence_id_2	duplicate_count_2	v_call_2	j_call_2	junction_aa_2
C	X	1	TCRBV07-09	TCRBJ01-02	CASSLRVGGFGYTF	B1	U	10	TCRBV07-09	TCRBJ01-02	CASSLRVGGFGYTF
C	X	1	TCRBV07-09	TCRBJ01-02	CASSLRVGGFGYTF	B1	T	5	TCRBV07-09	TCRBJ01-02	CASSLRVGGYGYTF
C	Y	1	TCRBV07-06	TCRBJ02-01	CASSTSHQQYF	B2	V	7	TCRBV07-06	TCRBJ02-01	CASSTSHQQYF
```

The results indicate that sequence X was found (twice) in repertoire
B1 (matching sequences T and U) and that sequence Y was found in
repertoire B2 (matching sequence V).


## Example 3: Clustering sequences

This time we will cluster the nucleotide sequences in the file
`setb.tsv` using the `-c` or `--cluster` option.

The command line to run is:

`compairr -c setb.tsv -d 1 -n -o output.tsv`

The output during the clustering is as follows:

```
CompAIRR 1.7.0 - Comparison of Adaptive Immune Receptor Repertoires
https://github.com/uio-bmi/compairr

Start time:        Thu Mar 03 12:33:05 CET 2022
Command (m/c/x):   Cluster (-c)
Repertoire:        setb.tsv
Nucleotides (n):   Yes
Differences (d):   1
Indels (i):        No
Ignore counts (f): No
Ignore genes (g):  No
Ign. unknown (u):  No
Threads (t):       1
Output file (o):   output.tsv
Log file (l):      (stderr)

Immune receptor repertoire clustering

Reading sequences: 100% (0s)
Repertoires:       2
Sequences:         3
Residues:          117
Shortest:          33
Longest:           42
Average length:    39.0
Total dupl. count: 22
Indexing:          100% (0s)

Unique V genes:    2
Unique J genes:    2

Computing hashes:  100% (0s)
Hashing sequences: 100% (0s)
Building network:  100% (0s)
Clustering:        100% (0s)
Sorting clusters:  100% (0s)
Writing clusters:  100% (0s)

Clusters:          2
End time:          Thu Mar 03 12:33:05 CET 2022
```

The result in the file `output.tsv` looks like this:

```tsv
#cluster_no	cluster_size	repertoire_id	sequence_id	duplicate_count	v_call	j_call	junction
1	2	B1	T	5	TCRBV07-09	TCRBJ01-02	tgcgcgagcagcctgcgcgtgggcggctatggctataccttt
1	2	B1	U	10	TCRBV07-09	TCRBJ01-02	tgcgcgagcagcctgcgcgtgggcggctttggctataccttt
2	1	B2	V	7	TCRBV07-06	TCRBJ02-01	tgcgcgagcagcaccagccatcagcagtatttt
```

In this case, there are 2 clusters. The first contains 2 sequences (T
and U from B1), while the second cluster contains 1 sequence (V from
B2). The sequences are clustered across repertoires.


## Example 4: Deduplication

This time we will deduplicate the amino acid sequences in the file
`setb.tsv` using the `-z` or `--deduplicate` option.

The command line to run is:

`compairr -z setb.tsv -o output.tsv`

The output will look like this:

```
CompAIRR 1.8.0 - Comparison of Adaptive Immune Receptor Repertoires
https://github.com/uio-bmi/compairr

Start time:        Thu Sep 15 17:10:51 CEST 2022
Command:           Deduplicate (--deduplicate)
Repertoire:        setb.tsv
Nucleotides (n):   No
Differences (d):   0
Indels (i):        No
Ignore counts (f): No
Ignore genes (g):  No
Ign. unknown (u):  No
Threads (t):       1
Output file (o):   output.tsv
Log file (l):      (stderr)

Reading sequences: 100% (0s)
Repertoires:       2
Sequences:         3
Residues:          39
Shortest:          11
Longest:           14
Average length:    13.0
Total dupl. count: 22
Indexing:          100% (0s)
Unique V genes:    2
Unique J genes:    2
Computing hashes:  100% (0s)
Deduplicating:     100% (0s)
Duplicates merged: 0
Writing output:    100% (0s)

End time:          Thu Sep 15 17:10:51 CEST 2022
```

The result in the file `output.tsv` looks like this:

```tsv
repertoire_id	duplicate_count	v_call	j_call	junction_aa
B1	5	TCRBV07-09	TCRBJ01-02	CASSLRVGGYGYTF
B1	10	TCRBV07-09	TCRBJ01-02	CASSLRVGGFGYTF
B2	7	TCRBV07-06	TCRBJ02-01	CASSTSHQQYF
```

There were no duplicates in this dataset so the output is essentially
identical to the input data, but does not include all the original
columns. If the two sequences in repertoire B1 had been identical, the
two lines would have been merged and the new `duplicate_count` would
have been 15.


## Implementation

The program is written in C++. The strategy for finding similar
sequences is based on a similar concept developed for the tool
[Swarm](https://github.com/torognes/swarm) (Mahé et al.
2021). Basically, a 64-bit hash is computed for all sequences in the
sets. All hashes for one set are stored in a Bloom filter and in a
hash table. We then look for matches to sequences in the second set by
looking them up in the Bloom filter and then, if there was a match, in
the hash table. To find matches with 1 or 2 substitutions or indels,
the hashes of all these variant sequences are generated and looked
up. When d>2, a different strategy is used where all sequences are
compared against each other and the number of differences is found.


## Performance

As a preliminary performance test, Cohort 2 ("Keck") of [the
dataset](https://s3-us-west-2.amazonaws.com/publishedproject-supplements/emerson-2017-natgen/emerson-2017-natgen.zip)
by Emerson et al. (2017) was compared to itself. It contains 120 repertoires
with a total of 24 205 557 extracted sequences. The test was performed
with CompAIRR version 1.3.1. The timing results are shown below.

Distance | Indels | Threads | Time (s) | Time (mm:ss)
-------: | :----: | ------: | -------: | -----------:
0 | no | 1 | 18 | 0:18
0 | no | 4 | 12 | 0:12
1 | no | 1 | 224 | 3:44
1 | no | 4 | 72 | 1:12
1 | yes | 1 | 367 | 6:07
1 | yes | 4 | 111 | 1:51
2 | no | 4 | 3200 | 53:20

When the distance is zero almost all of the time was used to read
files.

Memory usage was 2.5GB, corresponding to an average of about 100 bytes
per sequence.

Since this is a comparison of a repertoire set to itself, the dataset
is only read once, and the memory needed is also reduced as compared
to a situation were two different repertoire sets are compared.

Wall time and memory usage was measured by `/usr/bin/time`. The
analysis was performed on an Apple Mac Mini M1 (2020) with 16GB RAM.


## Benchmarking

The AIRR overlap functionality of CompAIRR has been thoroughly
benchmarked against similar tools. All data, scripts, and results are
available in a separate [CompAIRR benchmarking
repository](https://github.com/uio-bmi/compairr-benchmarking).


## Tips

If computer memory is limited, the dataset may be split into blocks
before running CompAIRR on each block separately. Results then needs
to be merged together again afterwards. This may be achieved with a
simple script. We will consider providing such a script.


## Development team

The code has been developed by Torbjørn Rognes based on code from
Swarm where Frédéric Mahé and Lucas Czech made important
contributions. Geir Kjetil Sandve had the idea of developing a tool
for rapid repertoire set comparison. Lonneke Scheffer has tested and
benchmarked the tool, and suggested new features. Milena Pavlovic and
Victor Greiff have also contributed to the project.


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


## Citing CompAIRR

Please cite the following if you use CompAIRR in any published work:

* Rognes T, Scheffer L, Greiff V, Sandve GK (2021) **CompAIRR: ultra-fast comparison of adaptive immune receptor repertoires by exact and approximate sequence matching.** *Bioinformatics*, btac505. doi: [10.1093/bioinformatics/btac505](https://doi.org/10.1093/bioinformatics/btac505)

The article is also available in preprint form:

* Rognes T, Scheffer L, Greiff V, Sandve GK (2021) **CompAIRR: ultra-fast comparison of adaptive immune receptor repertoires by exact and approximate sequence matching.** *bioRxiv*, 2021.10.30.466600. doi: [10.1101/2021.10.30.466600](https://doi.org/10.1101/2021.10.30.466600)


## References

* Emerson RO, DeWitt WS, Vignali M, Gravley J, Hu JK, Osborne EJ, Desmarais C, Klinger M, Carlson CS, Hansen JA, Rieder M, Robins HS (2017) **Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire.** *Nature Genetics*, 49 (5): 659-665. doi: [10.1038/ng.3822](https://doi.org/10.1038/ng.3822)

* Mahé F, Czech L, Stamatakis A, Quince C, de Vargas C, Dunthorn M, Rognes T (2021) **Swarm v3: Towards Tera-Scale Amplicon Clustering.** *Bioinformatics*, btab493. doi: [10.1093/bioinformatics/btab493](https://doi.org/10.1093/bioinformatics/btab493)
