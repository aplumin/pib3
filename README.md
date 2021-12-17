# Programming in Bioinformatics 3 - Aaron Plumin

For now, I left all the STAR output in the repository.

## TODO update environment.yml
## TODO complete run_me.sh
## TODO check .gitignore
## TODO write tests
## TODO check flake8
## TODO add type hints (Add type hints to your custom Python function/method signatures. It will be enough to only add type hints for all input arguments and the return values, although you are of course welcome to add them for any local variables as well)
## TODO check docstrings


## TODO instructions on how to deploy (e.g., set up a Conda environment or build a Docker image; see below)

## TODO instructions to run code

## Answers to exercises
### Session 1, Exercise 1.4
### Output:

Nucleotide fractions queries:

A: 0.2

C: 0.2

G: 0.3

T: 0.2

Nucleotide fractions references:

A: 0.2

C: 0.3

G: 0.3

T: 0.2

{'sequence1': {'chr1': [759], 'chr2': [], 'chr3': [], 'chr4': []}, 'sequence2': {'chr1': [], 'chr2': [1422], 'chr3': [], 'chr4': []}, 'sequence4': {'chr1': [], 'chr2': [1039], 'chr3': [1422], 'chr4': [1455]}}

### Interpretation:

- sequence3 contains non-DNA letters
- sequence1 and sequence 2 both are only present in one of the chromosomes (could be a gene, for example)
- sequence4 is much shorter and present once in each of chromosomes 2, 3, and 4

### Session 2, Exercise 2.1
How many alignments were reported?

18068

How many reads were uniquely mapped?

17672

How many reads were mapped to multiple loci?

180

How many reads could not be mapped?

2148 (actually, 2144 were unmapped and 4 mapped to too many loci - I left the max at the default of 10 - according to the log file Log.final.out, where all the numbers of unique, multi-mapping, and unmapped reads can also be found).

Compare the sum of uniquely mapped, multi-mapped and unmapped reads to the total number of reads in the FASTQ input files. Do the numbers match?

Yes (2*10000 - 17672 - 180 - 2148 = 0).
