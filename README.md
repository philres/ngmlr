### Quick start

Download [binary](https://github.com/philres/ngmlr/releases/tag/v0.2.3) from github and unzip or [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/ngmlr/README.html) or pull docker [![Docker Repository on Quay](https://quay.io/repository/philres/ngmlr/status "Docker Repository on Quay")](https://quay.io/repository/philres/ngmlr)

For Pacbio data run:
```bash
ngmlr -t 4 -r reference.fasta -q reads.fastq -o test.sam
```
For Oxford Nanopore run:
```bash
ngmlr -t 4 -r reference.fasta -q reads.fastq -o test.sam -x ont
```

### Intorduction
 
Ngmlr is a long-read mapper desigend to sensitively align PacBilo or Oxford Nanopore to (large) reference genomes. It was desigend to correctly align reads stemming from (complex) structural variations. Ngmlr uses an SV aware k-mer search to find approximate mapping locations for a read and a banded Smith-Waterman alignment algorithm with a non-affine gap model that penalizes gap extensions for longer gaps less than for shorter ones to compute precise alignments. The gap model allows ngmlr to account for both the sequencing error and real genomic variations at the same time and makes it especially effective at more precisely identifying the position of breakpoints stemming from (complex) structural variations. The k-mer search helps to detect and split reads that cannot be aligned linearly, enabling ngmlr to reliably align reads to a wide range of different structural variations including nested SVs (e.g. inversions flanked by deletions).
Currently ngmlr takes about 60 minutes (on a AMD Opteron 6348) and 10 GB RAM for aligning 1Gbp of Pacbio Reads when using 10 threads.

**Poster & Talks:**

[Accurate and fast detection of complex and nested structural variations using long read technologies](http://schatzlab.cshl.edu/presentations/2016/2016.10.28.BIODATA.PacBioSV.pdf)
Biological Data Science, Cold Spring Harbor Laboratory, Cold Spring Harbor, NY, 26 - 29.10.2016

[NGMLR: Highly accurate read mapping of third generation sequencing reads for improved structural variation analysis](http://www.cibiv.at/~philipp_/files/gi2016_poster_phr.pdf) 
Genome Informatics 2016, Wellcome Genome Campus Conference Centre, Hinxton, Cambridge, UK, 19.09.-2.09.2016

### Parameters

```
Usage: ngmlr [options] -r <reference> -q <reads> [-o <output>]

Input/Output:
    -r <file>,  --reference <file>
        (required)  Path to the reference genome (FASTA/Q, can be gzipped)
    -q <file>,  --query <file>
        (required)  Path to the read file (FASTA/Q)
    -o <file>,  --output <file>
        Path to output file [stdout]

General:
    -t <int>,  --threads <int>
        Number of threads [1]
    -x <pacbio, ont>,  --presets <pacbio, ont>
        Parameter presets for different sequencing technologies [pacbio]
    -i <0-1>,  --min-identity <0-1>
        Alignments with an identity lower than this threshold will be discarded [0.65]
    -R <int/float>,  --min-residues <int/float>
        Alignments containing less than <int> or (<float> * read length) residues will be discarded [50]
    --no-smallinv
        Do not detect small inversions [false]
    --no-lowqualitysplit
        Do not split alignments with poor quality [false]
    --verbose
        Debug output [false]
    --no-progress
        Do not print progress info while mapping [false]

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
    --subread-length <int>
        Length of fragments reads are split into [256]
    --subread-corridor <int>
        Length of corridor sub-reads are aligned with [40]
```

### Running with docker
```bash
docker run -ti -v /home/user/data/:/home/user/data/ quay.io/philres/ngmlr ngmlr -r /home/user/data/ref.fa -q /home/user/data/reads.fasta -o /home/user/data/output.sam
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
