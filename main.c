/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2019 Jeremy Wang
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <getopt.h>

#include "bnx.h"
#include "sim.h"
#include "digest.h"

void usage() {
  fprintf(stderr, "Usage: resim -f <ref.fasta> -r <CTTAAG> -x <100> -s <ground_truth.tsv> [options]\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -f: fasta: Reference sequence to simulate from\n");
  fprintf(stderr, "  -r: cutseq: Recognition/label site sequence (default: CTTAAG)\n");
  fprintf(stderr, "  -x: Simulated molecule coverage\n");
  fprintf(stderr, "  --break-rate: Probability of genome fragmentation per locus (default: 0.000005)\n");
  fprintf(stderr, "  --fn: Probability of missed label at true restriction site (default: 0.09893)\n");
  fprintf(stderr, "  --fp: Probability of false-positive label (default: 0.07558)\n");
  fprintf(stderr, "  --stretch-mean: Fragment stretch mean (default: 0.991385)\n");
  fprintf(stderr, "  --stretch-std: Fragment stretch standard deviation (default: 0.033733)\n");
  fprintf(stderr, "  --min-frag: Minimum detectable fragment size (default: 500)\n");
  fprintf(stderr, "  -s, --source-output: Output the reference positions of the simulated molecules to the given file\n");
}

static struct option long_options[] = {
// if these are the same as a single-character option, put that character in the 4th field instead of 0
  { "break-rate",             required_argument, 0, 0 },
  { "fn",                     required_argument, 0, 0 },
  { "fp",                     required_argument, 0, 0 },
  { "min-frag",               required_argument, 0, 0 },
  { "stretch-mean",           required_argument, 0, 0 },
  { "stretch-std",            required_argument, 0, 0 },
  { "help",                   no_argument,       0, 0 },
  { "source-output",          required_argument, 0, 0 },
  { 0, 0, 0, 0}
};

int main(int argc, char *argv[]) {
  srand(time(NULL));

  char* fasta_file = NULL; // .fasta file path/name
  char* restriction_seq = NULL; // restriction enzyme or label recognition sequence (must also be reverse complemented if not symmetrical)
  char* source_outfile = NULL; // output file for the truth/source positions
  int verbose = 0;

  float coverage = 0.0;
  int covg_threshold = 10;
  float break_rate = 0.000005; // one every 200Kb
  float fn = 0.09893; // based on empirical data for Saphyr DLE1
  float fp = 0.07558; // based on empirical data
  float min_frag = 500; // minimum reported gap between labels
  float stretch_mean = 0.991385; // these are based empirically on Cauchy distribution of NA12878 DLE1 data
  float stretch_std = 0.033733;

  int opt, long_idx;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "hf:r:vx:s:", long_options, &long_idx)) != -1) {
    switch (opt) {
      case 'h':
        usage();
        return 0;
        break;
      case 'f':
        fasta_file = optarg;
        break;
      case 'r':
        restriction_seq = optarg;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'x':
        coverage = atof(optarg);
        break;
      case 's':
        source_outfile = optarg;
        break;
      case '?':
        if (optopt == 'r' || optopt == 'f' || optopt == 'x' || optopt == 's')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      case 0:
        //fprintf(stderr, "option %s\n", long_options[long_idx].name);
        if (long_idx == 0) break_rate = atof(optarg); // --break-rate
        else if (long_idx == 1) fn = atof(optarg); // --fn
        else if (long_idx == 2) fp = atof(optarg); // --fp
        else if (long_idx == 3) min_frag = atoi(optarg); // --min-frag
        else if (long_idx == 4) stretch_mean = atof(optarg); // --stretch-mean
        else if (long_idx == 5) stretch_std = atof(optarg); // --stretch-std
        else if (long_idx == 7) {usage(); return 0;} // --help
        else if (long_idx == 8) source_outfile = optarg; // --source-output
        break;
      default:
        usage();
        return 1;
    }
  }

  int i = 0;
  int ret = 1;

  if(fasta_file == NULL) {
    fprintf(stderr, "FASTA file required (-f)\n\n");
    usage();
    return 1;
  }
  if(restriction_seq == NULL) {
    fprintf(stderr, "Restriction sequence is required (-r)\n\n");
    usage();
    return 1;
  } else {
    if(strcmp(restriction_seq, "DLE1") == 0 || strcmp(restriction_seq, "DLE-1") == 0) {
      restriction_seq = "CTTAAG";
    }
  }
  if(coverage < FLT_EPSILON) {
    fprintf(stderr, "Coverage is required (-x)\n\n");
    usage();
    return 1;
  }
  // make a list of restriction seqs - that's what digest wants
  char** rseqs = malloc(1 * sizeof(char*));
  rseqs[0] = restriction_seq;

  fprintf(stderr, "-- Running optical mapping simulation --\n");
  cmap c = simulate_bnx(fasta_file, rseqs, 1, break_rate, fn, fp, stretch_mean, stretch_std, min_frag, coverage);
  fprintf(stderr, "Done simulating, writing to BNX...\n");
  ret = write_bnx(&c, stdout);

  if(source_outfile != NULL) {
    fprintf(stderr, "Writing truth/source positions to '%s'\n", source_outfile);
    FILE *fp = fopen(source_outfile, "w");
    fprintf(fp, "ref_id\tstart_pos\n");
    for(i = 0; i < kv_size(c.source); i++) {
      fprintf(fp, "%u\t%u\n", kv_A(c.source, i).ref_id, kv_A(c.source, i).pos);
    }
    fclose(fp);
  }

  // free everything
  return ret;
}
