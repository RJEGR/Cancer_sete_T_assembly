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

3.1) Concatenate forward and reverse reads in a single bulk `bulk_reads.sh`

**Note**: if try to quantify using stringtie (instead of RSEM), needs to run hisat2 in a per-sample-aligment mode instead of concatenate all the reads in a single bulk. After it you can use samtools to merge and combine bam files after alignmente Ej: https://rnabio.org/module-02-alignment/0002/03/01/Alignment/

```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=cat
#SBATCH --error=slurm-%j.err
#SBATCH -N 1
#SBATCH --mem=12GB
#SBATCH --ntasks-per-node=20

cat *R1.P.qtrim.fq > reads_f.fq
cat *R2.P.qtrim.fq > reads_r.fq

exit 0
```

3.2) Build the hg (v 38) index (For human genome take 4.5 Gb of memory, > 2 hours)

```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=assembly
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 15-00:00:00
#SBATCH -p d15

# https://nbisweden.github.io/workshop-RNAseq/1906/lab_assembly.html#22_hisat2
# https://github.com/DaehwanKimLab/hisat2

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"
cd $SLURM_SUBMIT_DIR

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

3.3) testing different indexes:

Biostring [question](https://www.biostars.org/p/251741/): What is the difference between `genome`, `genome_tran` and `genome_snp_tran`

```bash
Genome is the basic index of the genome. genome_tran additionally includes annotated splicing boundaries. genome_snp_tranadditionally includes a number of SNPs, so you can (theoretically) get better alignment around them.

# cite: https://www.biostars.org/p/251741/

# vi urls:

# genome	
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# genome_snp	
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_snp.tar.gz
# genome_tran	
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz
# genome_snp_tran	
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_snptran.tar.gz
# genome_rep(above 2.2.0)	
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_rep.tar.gz
# genome_snp_rep(above 2.2.0)	
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_snprep.tar.gz

```

```bash
# cd /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/hisat2-build/snp_trans

srun cat urls | sh &> download.log &

```

Run in a single batch

```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=GRTA
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 15-00:00:00
#SBATCH -p d15

# Guide-Reference-Transcriptiome-Assembly
# M.C Ricardo Gomez-Reyes


echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"
cd $SLURM_SUBMIT_DIR

TOOL=/LUSTRE/apps/bioinformatica/hisat2-2.1.0/
export PATH=$TOOL:$PATH

STOOL=/LUSTRE/bioinformatica_data/genomica_funcional/bin/stringtie-2.2.1.Linux_x86_64
export PATH=$STOOL:$PATH

snp_trans=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Genomes/human/hg38ome/hisat2-build/snp_trans/grch38_snp_tran/

TOOL=/LUSTRE/bioinformatica_data/genomica_funcional/bin/gffread-0.12.7.Linux_x86_64
export PATH=$TOOL:$PATH

SAMTOOL=/LUSTRE/bioinformatica_data/genomica_funcional/cei/bin/samtools
export PATH=$SAMTOOL:$PATH

FPATH=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Genomes/human/hg38ome/hisat2-build/snp_trans/

snp_trans=${FPATH}/grch38_snp_tran


F=Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Outputs

out_sam=genome_snp_tran_vs_reads_f_and_reads_r_output.sam

# 1) Alignment

hisat2 --phred33 -p $SLURM_NPROCS -x $snp_trans/genome_snp_tran -1 reads_f.fq -2 reads_r.fq -S $out_sam

# 1.1 sam to bam conversion and sorting

samtools view -bS -o ${out_sam%.sam}.bam $out_sam

samtools sort -o ${out_sam%.sam}.sorted.bam ${out_sam%.sam}.bam

# 2) Assembly

stringtie ${out_sam%.sam}.sorted.bam -o transcripts.gtf

# 2.1 Reference annotation transcripts (-G)
# Flag -G is for include reference annotation to use for guiding the assembly process (GTF/GFF)

stringtie ${out_sam%.sam}.sorted.bam -B -G ${FPATH}/Homo_sapiens.GRCh38.84.gtf -o Homo_sapiens.GRCh38.84_transcripts.gtf

# 3) Generate a FASTA file with the DNA sequences for all transcripts in a GFF file.


# -w flag write a fasta file with spliced exons for each transcript

gffread -w transcripts.fa -g $FPATH/$F transcripts.gtf

echo "Fecha de paro: `date`"

exit

```



Integrate hisat to multiqc

https://multiqc.info/docs/#hisat2

### 4) StringTie

Assembly tools such as StringTie (Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. *Nature protocols*, *11*(9), 1650-1667.) ... use the gaps identified in the alignments to derive exon boundaries and possible splice sites. These de novo transcript assembly tools are particularly useful when the reference genome annotation may be missing or incomplete, or where aberrant transcripts (for example, in tumour tissue) are of interest ( Stark, R., Grzelak, M., & Hadfield, J. ([2019](https://www.nature.com/articles/s41576-019-0150-2)). RNA sequencing: the teenage years. *Nature Reviews Genetics*, *20*(11), 631-656) 

``` /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/StringTie```

https://ccb.jhu.edu/software/stringtie/

https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=StringTie
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

# mkdir StringTie
# cd StringTie

TOOL=/LUSTRE/bioinformatica_data/genomica_funcional/bin/stringtie-2.2.1.Linux_x86_64

export PATH=$TOOL:$PATH

# 

# which stringtie
# ln -s /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/output.sorted.bam .

# run StringTie

stringtie output.sorted.bam -o transcripts.gtf

exit 0
```

4.1) For the guided assembly results, you need first to extract the transcript sequences from the gtf transcript file. First load the compiled version of gffread tool:

```bash
TOOL=/LUSTRE/bioinformatica_data/genomica_funcional/bin/gffread-0.12.7.Linux_x86_64
export PATH=$TOOL:$PATH

SAMTOOL=/LUSTRE/bioinformatica_data/genomica_funcional/cei/bin/samtools
export PATH=$SAMTOOL:$PATH
```

4.2) Many bioinformatics programs represent genes and transcripts in GFF format (**G**eneral **F**eature **F**ormat) which simply describes the locations and the attributes of gene and transcript features on the genome (chromosome or scaffolds/contigs). GFF has many versions, but the two most popular that are **GTF2** and GTF3 (Gene Transfer Format). 

```bash
gffread -E transcripts.gtf | head
gffread -E transcripts.gtf -T -o- | head
```



4.3) `gffread` can also be used to generate a FASTA file with the DNA sequences for all transcripts in a GFF file. For this operation a fasta file with the genomic sequences have to be provided as well. Note that the retrieval of the transcript sequences this way is going to be much faster if a fasta index file (genome.fa.fai in this example) is found in the same directory with the genomic fasta file.  Such an index file can be created with the [**samtools**](http://samtools.sourceforge.net/) utility prior to running gffread, like this: `samtools faidx genome.fa`. By default gffread 0.12.7 version creates an index at first.

```bash
hg=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Genomes/human/hg38ome/hg38.fa

# take a few minutes

# -w flag write a fasta file with spliced exons for each transcript

srun gffread -w transcripts.fa -g $hg transcripts.gtf &>gffread.log
```

### 5) Assembly quality

#### 5.1) gff compare:  Evaluating transcript discovery accuracy

Gffcompare can be used to evaluate and compare the accuracy of transcript assemblers - in terms of their structural correctness (exon/intron coordinates). This assessment can even be performed in case of more generic "transcript discovery" programs like gene finders. The best way to do this would be to use a simulated data set (where the "reference annotation" is also the set of the expressed transcripts being simulated), but for well annotated reference genomes (human, mouse etc.), gffcompare can be used to evaluate and compare the general accuracy of isoform discovery programs on a real data set, using just the known (reference) annotation of that genome ([cite](http://ccb.jhu.edu/software/stringtie/gffcompare.shtml)).  Assume the StringTie's output is in `transcript.gtf`, while the reference annotation would be in a file called `Homo_sapiens.GRCh38.84.gtf`, the gffcompare commands would be:

```bash
# For interpretation Continue w/ http://ccb.jhu.edu/software/stringtie/gffcompare.shtml
TOOL=/LUSTRE/bioinformatica_data/genomica_funcional/bin/gffcompare-0.12.6.Linux_x86_64
export PATH=$TOOL:$PATH

# -r reference annotation file (GTF/GFF)
# -R for -r option, consider only the reference transcripts that overlap any of the input transfrags (Sn correction)
# -o is for prefix to add in the output files

# Ex. gffcompare -R -r Homo_sapiens.GRCh38.84.gtf -o gffcompare transcripts.gtf
gffcompare -R -r Homo_sapiens.GRCh38.84.gtf -o gffcompare transcripts.gtf

# A Sensitivity | Precision stats will print out in the *.stats file

# Download References

ENSEMBL_RELEASE=84

ENSEMBL_GRCh38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna

ENSEMBL_GRCh38_GTF=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}

F=Homo_sapiens.GRCh38.dna.primary_assembly.fa

wget ${ENSEMBL_GRCh38_BASE}/$F.gz

wget ${ENSEMBL_GRCh38_GTF}.gtf.gz
```



#### 5.2) BUSCO (Completeness)

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
out=${fasta%.fa}

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

After running busco, 

```bash
mkdir summaries
cp ./run_transcripts*_odb*/short_summary* summaries

module load R-3.3.1

python2.7 /home/rgomez/bin/busco-master/scripts/generate_plot.py --working_directory ./summaries/
cp summaries/busco_figure.R .
sed -i 's/_odb9//g' summaries/busco_figure.R
Rscript summaries/busco_figure.R

# Additionally, prepare full_Table
mkdir full_tables
cp ./run_transcripts*_odb*/full_table_transcripts* full_tables

```

scp

```bash
scp -r rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/snp_trans/summaries .
scp -r rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/StringTie/summaries summaries_genome
```



#### 5.3) Transrate

`sbatch transrate.sh transcripts.fa reads_f.fq reads_r.fq`	

Test the oyster River protocol versio  `git clone https://github.com/macmanes-lab/Oyster_River_Protocol.git`

```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=ORvrP
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

#tool=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/transrate-1.0.3-linux-x86_64
tool=/LUSTRE/bioinformatica_data/genomica_funcional/bin/Oyster_River_Protocol/software/orp-transrate
export PATH=$tool:$PATH

fasta=$1
left=$2
right=$3
out=${fasta%.fa}



transrate --assembly $fasta --left $left --right $right --threads $SLURM_NPROCS --output ${out}_transrate

exit
# 609,125,172
# -mcp 10,000,000 to 700000000
# testig solution to https://github.com/blahah/transrate/issues/204
# You need to add this line between lines 39 and 40 cmd << " -mcp 700000000" # maximum candidate pool size
# /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/transrate-1.0.3-linux-x86_64/lib/app/lib/transrate/snap.rb

```

/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/snp_trans/transrate

5.4) Strand spec

https://github.com/trinityrnaseq/trinityrnaseq/wiki/Examine-Strand-Specificity

### 6) Annotation

https://rjegr.github.io/Transcriptomics/markdown/trinotate

/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/annotation/dataBase

Because a problem during the `Build_Trinotate_Boilerplate_SQLite_db`  step (error at the SNOGG download) We will use a copy of an emtpy Trinotate.sqlite file from `/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/oyster_full_assembly/annotation`

```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=transrate
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24

# error en el batch con la descarga,

SOURCE=/LUSTRE/bioinformatica_data/genomica_funcional/scripts/exports

source $SOURCE

$NOTATE_AD/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate

$BLAST/makeblastdb -in uniprot_sprot.pep -dbtype prot
gunzip Pfam-A.hmm.gz
$HMM/hmmpress Pfam-A.hmm &> hmmpress.log

exit
```

prepare isoform genes relation 

```bash
grep '^>' transcripts.fa |  sed 's/>//g' > transcript.map
cat transcript.map | awk '{gsub(/\.[0-9]$/,"",$1); print $1}' > genes.map
paste genes.map transcript.map > genes_trans_map
rm *.map
```

Run trinotate `sbatch annot.sh transcript.fa`

`chmod +x annot.sh`

```bash
#!/bin/bash
### Directivas
#SBATCH -p cicese
#SBATCH --job-name=trinotate
#SBATCH --ntasks-per-node=24
#SBATCH -N 2
#SBATCH --output=trinotate-%j.log
#SBATCH --error=trinotate-%j.err
#SBATCH -t 6-00:00:00

fasta=$1

RUN=/LUSTRE/apps/bioinformatica/Trinotate/auto/

$RUN/autoTrinotate.pl --Trinotate_sqlite Trinotate.sqlite --transcripts $fasta --gene_to_trans_map genes_trans_map --conf conf.txt --CPU $SLURM_NPROCS

exit
```

At the end

```bash
TRINO=/LUSTRE/apps/bioinformatica/Trinotate/
export PATH=$TRINO:$PATH 

Trinotate Trinotate.sqlite init --gene_trans_map genes_trans_map --transcript_fasta transcripts.fa --transdecoder_pep transcripts.fa.transdecoder.pep

# 1) Resultados a nivel proteina:

Trinotate Trinotate.sqlite LOAD_swissprot_blastp swissprot.blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
Trinotate Trinotate.sqlite LOAD_signalp signalp.out
# Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out

# 2) Resultados a nivel RNA
Trinotate Trinotate.sqlite  LOAD_swissprot_blastx swissprot.blastx.outfmt6
# Trinotate Trinotate.sqlite LOAD_rnammer good.Trinity.fasta.rnammer.gff
Trinotate Trinotate.sqlite report > Trinotate.xls 

```



Then download

```bash
scp rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/annotation/Trinotate.xls .
```

### 7) Quantification

7.1) RSEM

```bash
ls ../*R2* | grep fq | sort > R2.tmp && ls ../*R1* | grep fq | sort > R1.tmp && \
cut -d "_" -f1 R1.tmp > factors.tmp && \
cut -d "_" -f1 R1.tmp > codes.tmp && \
paste factors.tmp codes.tmp R1.tmp R2.tmp | awk '{gsub(/\-/,"_",$2); \
print $1,$2,$3,$4}' | column -t > samples.file && rm *tmp

# sed 's/[.,/]//g' samples.file
```

`vi RSEM.sh`
`chmod +x RSEM.sh`
`sbatch RSEM.sh transcripts.fa samples.file`

/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/quantification

```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=RSEM
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 15-00:00:00
#SBATCH -p d15

#Defining paths for RSEM
BOWTIE2=/LUSTRE/apps/bioinformatica/bowtie2
UTILS=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.5.1/util
RSEM=/LUSTRE/bioinformatica_data/genomica_funcional/bin/RSEM
SAM2LS=/LUSTRE/apps/bioinformatica/samtools-1.7/bin

#exporting paths:
export PATH=$UTILS:$PATH
export PATH=$BOWTIE2:$PATH
export PATH=$RSEM:$PATH
export PATH=$SAM2LS:$PATH

fasta=$1
sam_file=$2

align_and_estimate_abundance.pl --transcripts $fasta --seqType fq \
    --est_method RSEM \
    --aln_method bowtie2 \
    --prep_reference \
    --samples_file $sam_file \
    --thread_count=$SLURM_NPROCS
    # --gene_trans_mapls

exit


```

Then, prepare matrix for Differential Expression Analysis

```bash
mkdir DiffExp
cd DiffExp

ls -d ../quantification/*/*.isoforms.results > isoforms.results

module load R-3.5.0_bio

UTILS=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.5.1/util

$UTILS/abundance_estimates_to_matrix.pl --gene_trans_map genes_trans_map --est_method RSEM --out_prefix abundance --quant_files isoforms.results --name_sample_by_basedir
```

And download

```bash
scp -r rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/DiffExp .
```

Additionally, prepera Exploratory data analysis

```bash
PTR=/LUSTRE/apps/bioinformatica/trinityrnaseq-Trinity-v2.5.1/Analysis/DifferentialExpression
export PATH=$PTR:$PATH

module load R-3.5.0_bio

sbatch PtR.sh abundance.isoform.counts.matrix sam_file
```

Then

```bash
scp -r rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/DiffExp/figures/ .
```

7.) Stringtie to Ballow (Para implementar estar version, es neceasrio generar el alineamiento por biblioteca en el paso de `hisat2 > sam2bam > merge sorted.bam` para generar la tabla de abundancia sample-to-sample) Por ejemplo: `/Users/cigom/Documents/DOCTORADO/stringTie`

https://rnabio.org/module-03-expression/0003/02/01/Expression/

https://www.bioconductor.org/packages/release/bioc/vignettes/ballgown/inst/doc/ballgown.html

```bash
/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/human_cancer/snp_trans/transcripts.gtf

# ‘-B’ enable output of Ballgown table files which will be created in the same directory as the output GTF (requires -G, -o recommended)
# ‘-A’ output path/file name for gene abundance estimates

#!/bin/sh
## Directivas
#SBATCH --job-name=StringTie
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

# mkdir StringTie
# cd StringTie

TOOL=/LUSTRE/bioinformatica_data/genomica_funcional/bin/stringtie-2.2.1.Linux_x86_64

export PATH=$TOOL:$PATH

stringtie genome_snp_tran_vs_reads_f_and_reads_r_output.sorted.bam -B -G Homo_sapiens.GRCh38.84.gtf -o Homo_sapiens.GRCh38.84_transcripts.gtf -A transcripts_abundance.tsv


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

Error during `/LUSTRE/bioinformatica_data/genomica_funcional/bin/Trinotate/util/rnammer_support/RnammerTranscriptome.pl --transcriptome transcripts.fa --path_to_rnammer /LUSTRE/bioinformatica_data/genomica_funcional/bin/rnammer/rnammer`



```bash
#!/bin/sh
## Directivas
#SBATCH --job-name=rnamer
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 15-00:00:00
#SBATCH -p d15


HMMER=/LUSTRE/bioinformatica_data/genomica_funcional/bin/RNAMMER2_TRINOTATE/hmmer-2.3

#exporting paths:
export PATH=$HMMER:$PATH

/LUSTRE/bioinformatica_data/genomica_funcional/bin/Trinotate/util/rnammer_support/RnammerTranscriptome.pl --transcriptome transcripts.fa --path_to_rnammer /LUSTRE/bioinformatica_data/genomica_funcional/bin/rnammer/rnammer
# /LUSTRE/bioinformatica_data/genomica_funcional/bin/RNAMMER2_TRINOTATE/hmmer-2.3

exit

```



```bash
SuperScaffold 100
acc: STRG.106564.1

Done.

CMD: perl /LUSTRE/bioinformatica_data/genomica_funcional/bin/rnammer/rnammer -S euk -m tsu,lsu,ssu -gff tmp.superscaff.rnammer.gff < transcriptSuperScaffold.fasta
$VAR1 = {
          'sequence' => './temp.2598.fsa',
          'global_fullmodel_score' => '0',
          'flankStop' => '4500',
          'mode' => 'silent',
          'flankBegin' => '4500',
          'feature' => 'rRNA',
          'domain_fullmodel_evalue' => '1e-05',
          'domain_fullmodel_score' => '0',
          'postmodel' => '/LUSTRE/bioinformatica_data/genomica_funcional/bin/rnammer/lib/euk.lsu.rnammer.hmm',
          'tempdir' => '.',
          'id' => '100000006187698955',
          'global_spottermodel_evalue' => '1e-05',
          'spottermodel' => '/LUSTRE/bioinformatica_data/genomica_funcional/bin/rnammer/lib/euk.lsu.rnammer.initial.hmm',
          'global_spottermodel_score' => '0',
          'description' => '28s_rRNA',
          'xml_output' => './2598.lsu.xml',
          'global_fullmodel_evalue' => '1e-05',
          'domain_spottermodel_score' => '0',
          'domain_spottermodel_evalue' => '1e-05',
          'hmmsearch' => '/LUSTRE/apps/bioinformatica/hmmer-2.3.2/bin/hmmsearch',
          'config' => './2598.lsu.cf'
        };
apply_revcompl(): entry 1 (transcriptSuperScaffold): 1223061233 bp
fasta entry: 1
STRAND: pos
./100000006187698955.1.pos.fsa
build_jobs(): running ./100000006187698955.1.pos.fsa on spotter model
write_fasta(): Writing sequence to './100000006187698955.1.pos.fsa'...(1223061233 bases)
CMD: /LUSTRE/apps/bioinformatica/hmmer-2.3.2/bin/hmmsearch --compat --domE 1e-05 --domT 0 -E 1e-05 -T 0 /LUSTRE/bioinformatica_data/genomica_funcional/bin/rnammer/lib/euk.lsu.rnammer.initial.hmm "./100000006187698955.1.pos.fsa" > "./100000006187698955.1.pos.fsa.hmmsearchresult"
```

