Ricardo Gómez-Reyes

### Tools

Esearch Installation: [Here](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

SRA Tool Kit Installation (For prefetch and fastq-dump) : [Here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)

### Get all SRA runs for a BioProject based on an SRA Run ID

```bash
# PRJNA450372, Single Cell RNA sequencing of Adult Human Breast Epithelial Cells
esearch -db sra -query "PRJNA450372" |  efetch -format docsum | xtract -pattern Runs -ACC @acc  -element "&ACC" > SraAccList.txt
```

### A list of Runs:

Proof-read: [Here](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump#how-to-use-prefetch-and-fasterq-dump-to-extract-fastq-files-from-sra-accessions)

```bash
# Test
head -n1 SraAccList.txt | xargs -n 1 -P 12 fastq-dump --split-files --gzip --skip-technical --outdir .
# run All downloads
fastq-dump --outdir . --split-files SraAccList.txt
```

