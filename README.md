Resim
-----
An efficient and flexible simulator for high-throughput optical maps

All code is distributed under the MIT license - see LICENSE.

    --------*------------*---*------------*----*----
        --------*------*---------*--------*-----*--
      ----*--------*-------------*---*--
                 ---------*-------**-------*-

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

Output
------

Resim outputs simulated molecules in BNX format (https://bionanogenomics.com/wp-content/uploads/2018/04/30038-BNX-File-Format-Specification-Sheet.pdf) to standard output. In most cases, stdout should be redirected to the desired output (.bnx) file:

    resim [options] > output.bnx

A short runtime log and any errors are output to stderr and can also be redirected to a file if desired:

    resim [options] > output.bnx 2> run.log

To generate the "ground truth" of simulated molecules, include the -s flag and a file name:

    resim [options] -s ground_truth.tsv > output.bnx

The ground truth is output in a tab-delimited format:

    ref_id	start_pos
    19	    45090529
    13	    128801843
    10	    139531341
    43	    18835145
    13	    163174323
    ...

Where ref\_id is the source reference sequence, 0-indexed in the order they appeared in the input FASTA file and start\_pos
is the proximal end (closer to 0) of the simulated molecule, regardless of the orientation of the simulated molecule.

Example
-------

To simulate 100x coverage from a reference genome (fasta) using the DLE-1 recognition site:

    resim -f <fasta> -r CTTAAG -x 100 -s <output_truth> > <output_bnx>

Simulate 10x coverage from the human reference genome with DLE-1 (should take <1 minute):

    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
    resim -f GCF_000001405.39_GRCh38.p13_genomic.fna.gz -r CTTAAG -x 10 -s GRCh38_rekit_10x_truth.tsv > GRCh38_rekit_10x.bnx
