### Quick start

Download [binary](https://github.com/philres/ngmlr/releases/tag/v0.2.6) from github and unzip or [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/ngmlr/README.html) or pull docker [![Docker Automated buil](https://img.shields.io/docker/automated/jrottenberg/ffmpeg.svg)](https://hub.docker.com/r/philres/ngmlr/). For updates follow [![Twitter URL](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&style=plastic)](https://twitter.com/philres1)

Download precompiled version:
```bash
wget https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz
tar xvzf ngmlr-0.2.7-linux-x86_64.tar.gz
cd ngmlr-0.2.7/
```

For PacBio data run:
```bash
ngmlr -t 4 -r reference.fasta -q reads.fastq -o test.sam
```
For Oxford Nanopore run:
```bash
ngmlr -t 4 -r reference.fasta -q reads.fastq -o test.sam -x ont
```

### Introduction
 
CoNvex Gap-cost alignMents for Long Reads (ngmlr) is a long-read mapper designed to sensitively align PacBio or Oxford Nanopore to (large) reference genomes. It was designed to quickly and correctly align the reads, including those spanning (complex) structural variations. Ngmlr uses an SV aware k-mer search to find approximate mapping locations for a read and then a banded Smith-Waterman alignment algorithm to compute the final alignment. Ngmlr uses a convex gap cost model that penalizes gap extensions for longer gaps less than for shorter ones to compute precise alignments. The gap model allows ngmlr to account for both the sequencing error and real genomic variations at the same time and makes it especially effective at more precisely identifying the position of breakpoints stemming from structural variations. The k-mer search helps to detect and split reads that cannot be aligned linearly, enabling ngmlr to reliably align reads to a wide range of different structural variations including nested SVs (e.g. inversions flanked by deletions).

With 10 cores (AMD Opteron 6348), ngmlr currently takes about 90 minutes and 10 GB RAM for aligning 3Gbp (~ 1x human data) of PacBio reads.


### Citation:
Please see and cite our paper:
https://www.nature.com/articles/s41592-018-0001-7


**Poster & Talks:**

[Accurate and fast detection of complex and nested structural variations using long read technologies](http://schatzlab.cshl.edu/presentations/2016/2016.10.28.BIODATA.PacBioSV.pdf)<br>
Biological Data Science, Cold Spring Harbor Laboratory, Cold Spring Harbor, NY, 26 - 29.10.2016

[NGMLR: Highly accurate read mapping of third generation sequencing reads for improved structural variation analysis](http://www.cibiv.at/~philipp_/files/gi2016_poster_phr.pdf)<br> 
Genome Informatics 2016, Wellcome Genome Campus Conference Centre, Hinxton, Cambridge, UK, 19.09.-2.09.2016

### Parameters

```
Usage: ngmlr [options] -r <reference> -q <reads> [-o <output>]

Input/Output:
    -r <file>,  --reference <file>
        (required)  Path to the reference genome (FASTA/Q, can be gzipped)
    -q <file>,  --query <file>
        Path to the read file (FASTA/Q) [/dev/stdin]
    -o <string>,  --output <string>
        Path to output file [stdout]
    --skip-write
        Don't write reference index to disk [false]
    --bam-fix
        Report reads with > 64k CIGAR operations as unmapped. Required to be compatible with the BAM format [false]
    --rg-id <string>
        Adds RG:Z:<string> to all alignments in SAM/BAM [none]
    --rg-sm <string>
        RG header: Sample [none]
    --rg-lb <string>
        RG header: Library [none]
    --rg-pl <string>
        RG header: Platform [none]
    --rg-ds <string>
        RG header: Description [none]
    --rg-dt <string>
        RG header: Date (format: YYYY-MM-DD) [none]
    --rg-pu <string>
        RG header: Platform unit [none]
    --rg-pi <string>
        RG header: Median insert size [none]
    --rg-pg <string>
        RG header: Programs [none]
    --rg-cn <string>
        RG header: sequencing center [none]
    --rg-fo <string>
        RG header: Flow order [none]
    --rg-ks <string>
        RG header: Key sequence [none]

General:
    -t <int>,  --threads <int>
        Number of threads [1]
    -x <pacbio, ont>,  --presets <pacbio, ont>
        Parameter presets for different sequencing technologies [pacbio]
    -i <0-1>,  --min-identity <0-1>
        Alignments with an identity lower than this threshold will be discarded [0.65]
    -R <int/float>,  --min-residues <int/float>
        Alignments containing less than <int> or (<float> * read length) residues will be discarded [0.25]
    --no-smallinv
        Don't detect small inversions [false]
    --no-lowqualitysplit
        Split alignments with poor quality [false]
    --verbose
        Debug output [false]
    --no-progress
        Don't print progress info while mapping [false]

Advanced:
    --match <float>
        Match score [2]
    --mismatch <float>
        Mismatch score [-5]
    --gap-open <float>
        Gap open score [-5]
    --gap-extend-max <float>
        Gap open extend max [-5]
    --gap-extend-min <float>
        Gap open extend min [-1]
    --gap-decay <float>
        Gap extend decay [0.15]
    -k <10-15>,  --kmer-length <10-15>
        K-mer length in bases [13]
    --kmer-skip <int>
        Number of k-mers to skip when building the lookup table from the reference [2]
    --bin-size <int>
        Sets the size of the grid used during candidate search [4]
    --max-segments <int>
        Max number of segments allowed for a read per kb [1]
    --subread-length <int>
        Length of fragments reads are split into [256]
    --subread-corridor <int>
        Length of corridor sub-reads are aligned with [40]
```

### Running with docker
```bash
docker run -ti -v /home/user/data/:/home/user/data/ philres/ngmlr ngmlr -r /home/user/data/ref.fa -q /home/user/data/reads.fasta -o /home/user/data/output.sam
```

### Building ngmlr from source
OS: Linux and Mac OSX:
Requirements: zlib-dev, cmake, gcc/g++ (>=4.8.2)

```bash
git clone https://github.com/philres/ngmlr.git
cd ngmlr/
mkdir -p build
cd build/
cmake ..
make

cd ../bin/ngmlr-*/
./ngmlr
```

### Building ngmlr for linux with docker
```bash
git clone https://github.com/philres/ngmlr.git
mkdir -p ngmlr/build
docker run -v `pwd`/ngmlr:/ngmlr philres/nextgenmaplr-buildenv bash -c "cd /ngmlr/build && cmake .. &&  make"
`pwd`/ngmlr/bin/ngmlr-*/ngmlr
```

### NGMLR progress information
Example:
```
Processed: 92198 (0.66), R/S: 37.44, RL: 8857, Time: 2.00 5.00 11.62, Align: 0.96, 490, 0.81
```

92198 reads were processed so far
66 % of the 92198 reads were mapped (with > 25 % of their bp mapped)
37.44 are mapped on average per second
8857 is the average read length so far

"Time" and "Align" are for debugging purpose and will be removed.

### Datasets used in the mansucript:
We provide the NGMLR aligned reads and the Sniffles calls for the data sets used:  

Arabidopsis trio: [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Arabidopsis_trio](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Arabidopsis_trio) . 

Genome in a Bottle trio: 
+ Mappings: [ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/) . 

+ SV calls: [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/GiaB/](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/GiaB/)

NA12878: [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/NA12878/](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/NA12878/) .  

SKBR3: [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Skbr3/](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Skbr3/) . 

