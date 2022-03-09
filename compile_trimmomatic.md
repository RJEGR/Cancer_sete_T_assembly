1) Make a quick view of the meta data:

```bash
ls *R2* | grep fastq | sort > R2.tmp && ls *R1* | grep fastq | sort > R1.tmp && \
cut -d "_" -f1 R1.tmp > factors.tmp && \
cut -d "_" -f1 R1.tmp > codes.tmp && \
paste factors.tmp codes.tmp R1.tmp R2.tmp | awk '{gsub(/\-/,"_",$2); \
print $1,$2,$3,$4}' | column -t > samples.file && rm *tmp
```

Columns 1 and 2 are used as factors between samples or experimental design)

```bash
head samples.file
```

Use gedit or other text editor to make file changes

2) Run Trimmomatic

```bash
#SBATCH --error=slurm-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00
#SBATCH -p cicese


TRIMMOMATIC=/LUSTRE/bioinformatica_data/genomica_funcional/bin/Trimmomatic-0.36
TRUSEQ=/home/rgomez


# /// loading Java current version
module load jdk1.8.0_60
for file in $(ls *R1*gz | grep fastq)
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_R*.fastq.gz}"
java -jar $TRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 \
    ${base}_R1_001.fastq.gz ${base}_R2_001.fastq.gz \
    ${base}_R1.P.qtrim.fq.gz ${base}_R1.UP.qtrim.fq.gz \
    ${base}_R2.P.qtrim.fq.gz ${base}_R2.UP.qtrim.fq.gz \
    ILLUMINACLIP:$TRUSEQ/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36 LEADING:5 TRAILING:5
done

exit 
```

Ejecutar `fastqc` de las bibliotecas que pasaron el filtro de calidad

```bash
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

mkdir -p fastqc
fastqc *.gz -t 24 --nogroup -o ./fastqc 

exit

#export PATH=/LUSTRE/apps/Anaconda/conda2/bin:$PATH
#source activate multiqc_py2.7
#multiqc ./fastqc/*zip -o ./multiqc
```

Which Illumina adapter would you to consider? Check TruSeq3-PE-2.fa file

2.1) Running cutadapter (if you know your adapter sequences. Example:)

```bash
#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
 
for file in $(ls *R1*gz | grep -e 'fastq' -e 'fq' -e 'gz')
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_R*.fastq.gz}"
./cutadapt -a ^TTGTACACACCGCCC...GTAGGTGAACCTGCRGAAGG -A ^CCTTCYGCAGGTTCACCTAC...GGGCGGTGTGTACAA --discard-untrimmed --cores 24 -o ${base}_R1.clipped.fastq.gz -p ${base}_R2.clipped.fastq.gz ${base}_R1_001.fastq.gz ${base}_R2_001.fastq.gz 
done

mkdir -p clipped

mv *clipped.fastq.gz ./clipped

exit 

```





