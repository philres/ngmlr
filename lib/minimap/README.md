## Introduction

Minimap is an *experimental* tool to efficiently find multiple approximate
mapping positions between two sets of long sequences, such as between reads and
reference genomes, between genomes and between long noisy reads. By default, it
is tuned to have high sensitivity to 2kb matches around 20% divergence but with
low specificity. Minimap does not generate alignments as of now and because of
this, it is usually tens of times faster than mainstream *aligners*. With four
CPU cores, minimap can map 1.6Gbp PacBio reads to human in 2.5 minutes, 1Gbp
PacBio E. coli reads to pre-indexed 9.6Gbp bacterial genomes in 3 minutes, to
pre-indexed >100Gbp nt database in ~1 hour (of which ~20 minutes are spent on
loading index from the network filesystem; peak RAM: 10GB), map 2800 bacteria
to themselves in 1 hour, and map 1Gbp E. coli reads against themselves in a
couple of minutes.

Minimap does not replace mainstream aligners, but it can be useful when you
want to quickly identify long approximate matches at moderate divergence among
a huge collection of sequences. For this task, it is much faster than most
existing tools.

## Usage

* Map two sets of long sequences:
  ```sh
  minimap target.fa.gz query.fa.gz > out.mini
  ```
  The output is TAB-delimited with each line consisting of query name, length,
  0-based start, end, strand, target name, length, start, end, the number of
  matching bases, the number of co-linear minimizers in the match and the
  fraction of matching bases.

* All-vs-all PacBio read self-mapping for [miniasm][miniasm]:
  ```sh
  minimap -Sw5 -L100 -m0 reads.fa reads.fa | gzip -1 > reads.paf.gz
  ```

* Prebuild index and then map:
  ```sh
  minimap -d target.mmi target.fa.gz
  minimap -l target.mmi query.fa.gz > out.mini
  ```
  Minimap indexing is very fast (1 minute for human genome; 50 minutes for >100Gbp
  nt database retrieved on 2015-09-30), but for huge
  repeatedly used databases, prebuilding index is still preferred.

* Map sequences against themselve without diagnal matches:
  ```sh
  minimap -S sequences.fa sequences.fa > self-match.mini
  ```
  The output may still contain overlapping matches in repetitive regions.

## Algorithm Overview

1. Indexing. Collect all [(*w*,*k*)-minimizers][mini] in a batch (**-I**=4
   billion bp) of target sequences and store them in a hash table. Mark top
   **-f**=0.1% of most frequent minimizers as repeats. Minimap
   uses [invertible hash function][invhash] to avoid taking ploy-A as
   minimizers.

2. For each query, collect all (*w*,*k*)-minimizers and look up the hash table for
   matches (*q<sub>i</sub>*,*t<sub>i</sub>*,*s<sub>i</sub>*), where
   *q<sub>i</sub>* is the query position, *t<sub>i</sub>* the target position
   and *s<sub>i</sub>* indicates whether the minimizer match is on the same
   strand.

3. For matches on the same strand, sort by {*q<sub>i</sub>*-*t<sub>i</sub>*}
   and then cluster matches within a **-r**=500bp window. Minimap merges
   two windows if **-m**=50% of minimizer matches overlap. For matches on different
   strands, sort {*q<sub>i</sub>*+*t<sub>i</sub>*} and apply a similar
   clustering procedure. This is inspired by the [Hough transformation][hough].

4. For each cluster, sort (*q<sub>i</sub>*,*t<sub>i</sub>*) by *q<sub>i</sub>*
   and solve a [longest increasing sequence problem][lis] for *t<sub>i</sub>*. This
   finds the longest co-linear matching chain. Break the chain whenever there
   is a gap longer than **-g**=10000.

5. Output the start and end of the chain if it contains **-c**=4 or more
   minimizer matches and the matching length is no less than **-L**=40.

6. Go to 1 and rewind to the first record of query if there are more target
   sequences; otherwise stop.

To increase sensitivity, we may decrease **-w** to index more minimizers;
we may also decrease **-k**, though this may greatly impact performance for
mammalian genomes.

Also note that by default, if the total length of target sequences is less than
4Gbp (1G=1 billion; controlled by **-I**), minimap creates one index and stream
all the query sequences in one go. The multiple hits of a query sequence is
adjacent to each other in the output. If the total length is greater than
4Gbp, minimap needs to read query sequences multiple times. The multiple hits
of a query may not be adjacent.

[mini]: http://bioinformatics.oxfordjournals.org/content/20/18/3363.abstract
[lis]: https://en.wikipedia.org/wiki/Longest_increasing_subsequence
[hough]: https://en.wikipedia.org/wiki/Hough_transform
[invhash]: https://gist.github.com/lh3/974ced188be2f90422cc
[miniasm]: https://github.com/lh3/miniasm
