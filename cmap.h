#ifndef _cmap_h_
#define _cmap_h_

#include <stdint.h>
#include "klib/kvec.h"
#include "klib/kstring.h"

#define string_begins_with(s, pre) (strlen(pre) <= strlen(s) && strncmp(pre, s, strlen(pre)) == 0)

// various vectors to handle fragment/label data
typedef kvec_t(uint8_t) byteVec;

typedef kvec_t(uint32_t) u32Vec;
typedef kvec_t(u32Vec*) fragVec;

typedef kvec_t(char*) seqVec;
typedef kvec_t(char*) cstrVec;

typedef struct ref_pos {
  uint32_t ref_id;
  uint32_t pos;
} ref_pos;

typedef kvec_t(ref_pos) posVec;

/*

# CMAP File Version:	0.1
# Label Channels:	1
# Nickase Recognition Site 1:	CTTAAG
# Number of Consensus Nanomaps:	66
#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence
#f int	float	int	int	int	float	float	int	int
1	195471971.0	44559	1	1	3004108.0	1.0	1	1
1	195471971.0	44559	2	1	3010894.0	1.0	1	1
...

*/

typedef struct label {
  uint32_t position;
  float stdev;
  uint16_t coverage;
  uint8_t channel;
  uint16_t occurrence;
} label;

typedef struct molecule {
  uint32_t id;
  size_t n_labels;
  size_t length;
  label* labels;
} molecule;

typedef struct cmap {
  molecule* molecules;
  uint32_t n_maps;
  char** rec_seqs;
  uint32_t n_rec_seqs;
  posVec source;
} cmap;

// cmap and associated IO functions
char* get_val(char* buf);
void next_line(FILE *fp, char *buf, size_t bufsize);

int write_cmap(cmap *c, FILE* fp);
cmap read_cmap(const char* fn);
int add_map(cmap* c, uint32_t molid, uint32_t* positions, uint32_t n_pos, uint8_t channel);
void init_cmap(cmap* c);

size_t filter_labels(label* labels, size_t n_labels, label* filtered_labels, int resolution_min);

#endif
