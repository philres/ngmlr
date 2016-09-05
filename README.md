### Install

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

### Quick start

```bash
wget http://www.cibiv.at/~philipp_/files/ngmlr-testseq.tgz
tar xvzf ngm-testseq.tgz
ngmlr -r dh10b_ecoli.fasta -q dh10b_ecoli.fasta_pacbio.fastq -o test.sam
```