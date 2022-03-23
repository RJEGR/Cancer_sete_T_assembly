[TOC]



### 1) Make a quick view of the meta data:

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

Count the number of reads per lib

```bash
for file in $(ls *_L001_R1_001*gz | grep -e 'fastq' -e 'fq' -e 'gz')
do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_R*.fastq.gz}"
#zcat ${base}_R1_001.fastq.gz | grep -c "^@" ;
#zcat ${base}_R1.clipped.fastq.gz | grep -c "^@" ; 
echo $base'
done

# | awk 'BEGIN { FS = "-" } ; { print $3"."$2 }

```

Use gedit or other text editor to make file changes

### 2) Run Trimmomatic

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

2.2) Fastp tool (as alternative to trimmomatic)

```bash
# A tool designed to provide fast all-in-one preprocessing for FastQ files
# https://github.com/OpenGene/fastp#fastp

TOOL=/LUSTRE/bioinformatica_data/genomica_funcional/bin/

export PATH=$TOOL:$PATH

# for paired end data (gzip compressed)
fastp -i P1283_S30_L001_R1.P.qtrim.fq -I P1283_S30_L001_R2.P.qtrim.fq -o out.R1.fq.gz -O out.R2.fq.gz
```

Then , report in multiqc document the result

https://multiqc.info/docs/

Prepare dataset for dataviz

```bash
ls *R2.P.qtrim.fq | sort > R2.tmp && ls *R1.P.qtrim.fq | sort > R1.tmp && \
cut -d "_" -f1 R1.tmp > factors.tmp && \
cut -d "_" -f1 R1.tmp > codes.tmp && \
paste factors.tmp codes.tmp R1.tmp R2.tmp | awk '{gsub(/\-/,"_",$2); \
print $1,$2,$3,$4}' | column -t > samples.file && rm *tmp
```

```bash

grep 'Input Read Pairs' slurm-171018.err | awk '{print $4, $7, $12, $17, $20}' > trimmomatic.log

#. or


grep phred33 slurm-171018.err | awk '{gsub(/\_/," ", $2); print $2}' | awk '{print $1}' > trimmomatic_id.log

paste trimmomatic_id.log trimmomatic.log > trimmomatic.csv

# Download and enter to Rstudios


```

```R
c("ID",	"Input Read Pairs",	"Both Surviving",	"Forward Only",	"Reverse Only",	"Dropped")

```

### 3) Test genome-guide assembly (hisat2)

3.1) Concatenate forward and reverse reads in a single batch

```bash
cat *.R1.P.qtrim.fq > reads_f.fq
cat *.R2.P.qtrim.fq > reads_r.fq
```

3.2) Build the hg (v 38) index (For human genome take 4.5 Gb of memory, > 2 hours)

```bash
# https://nbisweden.github.io/workshop-RNAseq/1906/lab_assembly.html#22_hisat2
# https://github.com/DaehwanKimLab/hisat2

TOOL=/LUSTRE/apps/bioinformatica/hisat2-2.1.0/

export PATH=$TOOL:$PATH

# Build the hg (v 38) index (Already done, 22/03/22)
# https://www.ncbi.nlm.nih.gov/genome/51?genome_assembly_id=582967

hg=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Genomes/human/hg38ome/hg38.fa

hg_name=`basename ${hg%.fa}`

# hisat2-build -p $SLURM_NPROCS $hg ${hg%.fa}

# ll -h /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Genomes/human/hg38ome/*ht2

# Run on paired-end reads

hisat2 --phred33 -p $SLURM_NPROCS -x ./hisat2-build/$hg_name -1 reads_f.fq -2 reads_r.fq -S output.sam

# test single at:
## srun hisat2 --phred33 -p 24 -x ./hisat2-build/$hg_name -1 P1283_S30_L001_R1.P.qtrim.fq  -2 P1283_S30_L001_R2.P.qtrim.fq -S output.sam 2> hisat2.log &

mkdir hisat2
mv output.sam hisat2
cd hisat2

# convert from sam to bam as StringTie requiere the bam format to assemble the reads into a single transcript

samtools view -bS -o output.bam output.sam

samtools sort -o output.sorted.bam output.bam

```

Integrate hisat to multiqc

https://multiqc.info/docs/#hisat2

### 4) StringTie

https://ccb.jhu.edu/software/stringtie/

```bash
cd ..
mkdir stringtie

# run StringTie

stringtie hisat2/accepted_hits.sorted.bam -o stringtie/transcripts.gtf


```

For the guided assembly results, you need first to extract the transcript sequences from the gtf transcript file :

```bash
# Option 1
wget https://raw.githubusercontent.com/NBISweden/AGAT/bf48b0d7bc18ab204d5880545acc4b66c1ef13b7/bin/agat_sp_extract_sequences.pl .

chmod +x agat_sp_extract_sequences.pl

# Can't locate Clone.pm in @INC (you may need to install the Clone module)

# Option 2

# get the chosen AGAT container version
singularity pull docker://quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0 
# run the container
singularity run agat_0.8.0--pl5262hdfd78af_0.sif
# 
agat_convert_sp_gxf2gxf.pl --help
```



#### 4.1) Assembly quality

```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=BUSCOpy
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

module load python-2.7-anaconda

fasta=$1
out=${fasta%.fasta}

BUSCO=/home/rgomez/bin/busco-master/scripts/
export PATH=$BUSCO:$PATH

eukaryota_odb9=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oyster-rawdata/method_v2/ANNOTATE/BUSCOdb/eukaryota_odb9
mammalia_odb9=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oyster-rawdata/method_v2/ANNOTATE/BUSCOdb/mammalia_odb9
bacteria_odb9=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oyster-rawdata/method_v2/ANNOTATE/BUSCOdb/bacteria_odb9

run_BUSCO.py -i $fasta -l $eukaryota_odb9 -m transcriptome -o ${out}_eukaryota_odb9 -c 24
run_BUSCO.py -i $fasta -l $mammalia_odb9 -m transcriptome -o ${out}_mammalia_odb9 -c 24
run_BUSCO.py -i $fasta -l $bacteria_odb9 -m transcriptome -o ${out}_bacteria_odb9 -c 24

exit

# https://busco.ezlab.org/busco_userguide.html#interpreting-the-results
```

#### 4.2)

```bash

```



### 5) Annotation

```bash

```

### 6) Quantification

```bash

```





Trinity (Failed to run)

```bash
module load gcc-7.2

export PATH=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/:$PATH
export PATH=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/trinity-plugins/jellyfish-2.2.6/bin/:$PATH

which Trinity

srun Trinity --seqType fq --max_memory 100G --samples_file samples.file --no_normalize_reads --CPU 24 --output trinity_out --no_salmon & 2> trinity.log &

 /LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util/insilico_rea
d_normalization.pl --seqType fq --JM 100G  --max_cov 200 --min_cov 1 --CPU 24 --output /LUSTRE/bioinformatica_data/ecol
ogia_molecular/Francesco/Cancer/Trimmomatic/trinity_outputs/insilico_read_normalization --max_CV 10000  --left FILES ---right FILES --pairs_together --PARALLEL_STATS


--pairs_together --PARALLEL_STATS   died with ret 6400 at /LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/Trinity line 2745
        main::process_cmd('/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util...') called at /LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/Trinity line 3295
        main::normalize('/LUSTRE/bioinformatica_data/ecologia_molecular/Francesco/Canc...', 200, 'ARRAY(0x2884cf0)', 'ARRAY(0x2884d20)') called at /LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/Trinity line 3238
        main::run_normalization(200, 'ARRAY(0x2884cf0)', 'ARRAY(0x2884d20)') called at /LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/Trinity line 1350
        
        
        
        
        
        
CMD: /LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util/..//Inchworm/bin/fastaToKmerCoverageStats --reads left.fa --kmers jellyfish.K25.min2.kmers.fa --kmer_size 25  --num_threads 12  --DS  > left.fa.K25.stats
CMD: /LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.8.5/util/..//Inchworm/bin/fastaToKmerCoverageStats --reads right.fa --kmers jellyfish.K25.min2.kmers.fa --kmer_size 25  --num_threads 12  --DS  > right.fa.K25.stats
-reading Kmer occurrences...-reading Kmer occurrences...

terminate called after throwing an instance of 'std::bad_alloc'
```

