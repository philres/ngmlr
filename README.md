### Quick start

Download [binary](https://github.com/philres/nextgenmap-lr/releases/tag/0.1.5) from github
```bash
tar xvzf ngmlr-0.1.5.tar.gz
ngmlr-0.1.5/ngmlr -t 4 -r reference.fasta -q reads.fastq -o test.sam
```

### Intorduction
 
NextGenMap-LR is a long-read mapper desigend to sensitively align PacBilo or Oxford Nanopore to (large) reference genomes. It was desigend to correctly align reads stemming from (complex) structural variations. NextGenMap-LR uses an SV aware k-mer search to find approximate mapping locations for a read and a banded Smith-Waterman alignment algorithm with a non-affine gap model that penalizes gap extensions for longer gaps less than for shorter ones to compute precise alignments. The gap model allows NextGenMap-LR to account for both the sequencing error and real genomic variations at the same time and makes it especially effective at more precisely identifying the position of breakpoints stemming from (complex) structural variations. The k-mer search helps to detect and split reads that cannot be aligned linearly, enabling NextGenMap-LR to reliably align reads to a wide range of different structural variations including nested SVs (e.g. inversions flanked by deletions).
Currently NextGenMap-LR takes about 60 minutes (on a AMD Opteron 6348) and 10 GB RAM for aligning 1Gbp of Pacbio Reads when using 10 threads.

##### Poster
[NextGenMap-LR: Highly accurate read mapping of third generation sequencing reads for improved structural variation analysis](http://www.cibiv.at/~philipp_/files/gi2016_poster_phr.pdf) 
Genome Informatics 2016, Wellcome Genome Campus Conference Centre, Hinxton, Cambridge, UK, 19.09.-2.09.2016.

### Parameters

```
Usage: ngmlr [options] -r <reference> -q <reads> [-o <output>]

Input/Output:
    -r <file>,  --reference <file>
        (required)  Path to the reference genome (FASTA/Q, can be gzipped)
    -q <file>,  --query <file>
        (required)  Path to the read file (FASTA/Q, SAM/BAM)
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
    --match <int>
        Match score [1]
    --mismatch <int>
        Mismatch score [-4]
    --gap-open <int>
        Gap open score [-2]
    --gap-extend <int>
        Gap open extend [-2]
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
    --vcf <file>
        SNPs will be taken into account when building reference index [none]
    --bed-filter <file>
        Only reads in the regions specified by the BED file are read from the input file (requires BAM input) [none]
```

### Building NextGenMap-LR from source
OS: Linux and Mac OSX (experimental):
Requirements: zlib-dev, cmake, gcc/g++ (>=5.1)

```bash
git clone https://github.com/philres/nextgenmap-lr.git
cd nextgenmap-lr/
mkdir -p build
cd build/
cmake ..
make

cd ../bin/ngmlr-*/
./ngmlr
```