cd src/resources

# Create index
STAR  --runMode genomeGenerate --genomeFastaFiles Mus_musculus.GRCm38.dna_rm.chr19.fa --sjdbGTFfile Mus_musculus.GRCm38.88.chr19.gtf --genomeSAindexNbases 11

# Align
STAR --readFilesIn reads.mate_1.fq,reads.mate_2.fq

# How many alignments were reported?
grep -v "^@" Aligned.out.sam | grep -v "^$" | wc -l

# How many reads were uniquely mapped?
grep -v "^@" Aligned.out.sam | grep -v "^$" | cut -f12- -d$'\t' | grep "NH:i:1" | wc -l

# How many reads were mapped to multiple loci?
str=$(grep -v "^@" Aligned.out.sam | grep -v "^$" | cut -f12- -d$'\t' | grep -v "NH:i:1" | sed 's/^.*\(NH:i:\)//' | sort)
last_line="${str##*$'\n'}"
declare -i highest_multimap=${last_line:0:1}
declare -i multimaps=0
for (( i=2; i<=$highest_multimap; i++ ))
do
  multimaps=$(((((($(grep -v "^@" Aligned.out.sam | grep -v "^$" | cut -f12- -d$'\t' | grep -v "NH:i:1" | sed 's/^.*\(NH:i:\)//' | sort | grep "^$i" | wc -l) / $i)) + $multimaps))))
done
echo $multimaps

# How many reads could not be mapped? 
declare -i total_reads=$(($(grep "^@" reads.mate_1.fq | wc -l) + $(grep "^@" reads.mate_2.fq | wc -l)))
echo $(($total_reads - 17672 - 180 ))

# Convert SAM to FASTA file
cd ..
python sam_to_fasta.py resources/Aligned.out.sam resources/Aligned.out.fa

# Apply functions to the FASTA file
python process_fasta.py resources/sequences.fasta resources/genome.fasta

# flake8
cd ..
flake8

# pytest
# TODO