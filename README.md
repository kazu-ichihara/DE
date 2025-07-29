# README for DE

2025-7-29 Kazuya Ichihara

#### Python packages

* Python 3.7.0
* Numpy 1.17.3
* Pysam 0.15.3
* Scipy 1.4.0rc1
* Pandas 0.25.3
* Biopython 1.75
* OpenCV 3.4.2

#### Obtain genome annotation files

* Download human genome and transcriptome annotation files from GENCODE (Harrow et al., 2012).

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.transcripts.fa.gz
```

#### Prepare annotation files

```bash
python CreateAnnotation.py
```

#### Create metagene plots

```
python metagene.py Ribo.bam ./Ribo
python readdist.py Ribo.bam ./Ribo
python readend.py Ribo.bam ./Ribo
python regioncount.py Ribo.bam ./Ribo
python regioncount_DE.py Ribo.bam ./Ribo
```