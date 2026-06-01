# Genome sequencing Longest Increasing Subsequence algorithm

Place .fasta data in the dataset directory

Program build:

```cmake
mkdir -p build && cmake .. && make
```

.fasta file and program executable have to be within the same directory

Launch program with

```bash
./fasta_lis_sequencing -a local -m 2 -n -2 -g -2 -k 15 -w 10 -f 0.0001 > dataset/alignments.paf
```
