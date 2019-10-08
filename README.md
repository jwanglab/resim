resim
-----
An efficient and flexible simulator for high-throughput optical maps

All code is distributed under the MIT license.

Requirements
------------

zlib (https://zlib.net/)

To install on Ubuntu:

    sudo apt install zlib1g-dev

Installation
------------

    git clone https://github.com/txje/resim && cd resim
    git clone https://github.com/attractivechaos/klib
    make

Usage
-----

    Usage: resim -f <ref.fasta> -r <CTTAAG> -x <100> -s <ground_truth.tsv> [options]
    Options:
      -f: fasta: Reference sequence to simulate from
      -r: cutseq: Recognition/label site sequence (default: CTTAAG)
      -x: Simulated molecule coverage
      --break-rate: Probability of genome fragmentation per locus (default: 0.000005)
      --fn: Probability of missed label at true restriction site (default: 0.09893)
      --fp: Probability of false-positive label (default: 0.07558)
      --stretch-mean: Fragment stretch mean (default: 0.991385)
      --stretch-std: Fragment stretch standard deviation (default: 0.033733)
      --min-frag: Minimum detectable fragment size (default: 500)
      -s, --source-output: Output the reference positions of the simulated molecules to the given file

Example
-------

To simulate 100x coverage from a reference genome (fasta) using the DLE-1 recognition site:

    resim -f <fasta> -r CTTAAG -x 100 -s <output_truth> > <output_bnx>

Simulate 10x coverage from the human reference genome with DLE-1 (should take <1 minute):

    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
    resim -f GCF_000001405.39_GRCh38.p13_genomic.fna.gz -r CTTAAG -x 10 -s GRCh38_rekit_10x_truth.tsv > GRCh38_rekit_10x.bnx
