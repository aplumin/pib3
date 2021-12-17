# Programming in Bioinformatics, Part 3 - Aaron Plumin
These sequencing exercises involve some fundamental manipulations and analyses of DNA/RNA sequencing data.

## Prerequisites
The script requires bash on a Unix terminal.

It is easiest to use a `conda` environment (see https://docs.conda.io/en/latest/miniconda.html for details).
The environment dependencies are listed in the `environment.yml` file. 

Sequence alignments are done by the software `STAR` (see https://github.com/alexdobin/STAR), which is available via 
the conda package manager. To install `STAR`, run: `conda install -c bioconda star`, or alternatively refer to the 
distribution.

## Run the code
All the code for the exercises is included in the shell script `run_me.sh`,
which can be run with `bash run_me.sh` on a terminal with bash.


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

{'sequence1': {'chr1': [759], 'chr2': [], 'chr3': [], 'chr4': []}, 'sequence2': {'chr1': [], 'chr2': [1422], 
'chr3': [], 'chr4': []}, 'sequence4': {'chr1': [], 'chr2': [1039], 'chr3': [1422], 'chr4': [1455]}}

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

2148 (actually - according to the log file Log.final.out, where all the numbers of unique, multi-mapping, and unmapped 
reads can also be found - 2144 were unmapped and 4 mapped to too many loci; I left the max at the default of 10).

Compare the sum of uniquely mapped, multi-mapped and unmapped reads to the total number of reads in the FASTQ input 
files. Do the numbers match?

Yes (2*10000 - 17672 - 180 - 2148 = 0).
