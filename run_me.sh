# Create index
STAR \
  --runMode genomeGenerate \
  --genomeFastaFiles data/Mus_musculus.GRCm38.dna_rm.chr19.fa \
  --sjdbGTFfile data/Mus_musculus.GRCm38.88.chr19.gtf \
  --genomeDir results/genomeDir/ \
  --genomeSAindexNbases 11

# Align
STAR \
  --genomeDir results/genomeDir/ \
  --readFilesIn data/reads.mate_1.fq data/reads.mate_2.fq \
  --outFileNamePrefix results/

# How many alignments were reported?
grep -v "^@" results/Aligned.out.sam | grep -c -v "^$"

# How many reads were uniquely mapped?
grep -v "^@" results/Aligned.out.sam | grep -v "^$" | cut -f12- -d$'\t' | grep -c "NH:i:1"

# How many reads were mapped to multiple loci?
str=$(grep -v "^@" results/Aligned.out.sam | grep -v "^$" | cut -f12- -d$'\t' | grep -v "NH:i:1" | sed 's/^.*\(NH:i:\)//' | sort)
last_line="${str##*$'\n'}"
declare -i highest_multimap=${last_line:0:1}
declare -i multimaps=0
for (( i=2; i<=highest_multimap; i++ ))
do
  multimaps=$(((((($(grep -v "^@" results/Aligned.out.sam | grep -v "^$" | cut -f12- -d$'\t' | grep -v "NH:i:1" | sed 's/^.*\(NH:i:\)//' | sort | grep -c "^$i") / i)) + multimaps))))
done
echo $multimaps

# How many reads could not be mapped? 
declare -i total_reads=$(($(grep -c "^@" data/reads.mate_1.fq) + $(grep -c "^@" data/reads.mate_2.fq)))
echo $((total_reads - 17672 - 180 ))

# Convert SAM to FASTA file
python src/sam_to_fasta.py results/Aligned.out.sam results/Aligned.out.fa

# Apply functions to the FASTA file
python src/process_fasta.py results/Aligned.out.fa data/genome.fasta

# flake8
flake8

# pytest
coverage run --source src/ -m pytest && coverage report -m
