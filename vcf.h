#ifndef __VCF_H__
#define __VCF_H__


#include <zlib.h>

#include "snp.h"
#include "chrom.h"

#define VCF_MAX_QUAL 1024
#define VCF_MAX_FILTER 1024
#define VCF_MAX_FORMAT 1024
#define VCF_N_CHROM_INIT 25

typedef struct {
  long n_samples;
  long n_geno_prob_col;
  long n_haplo_col;
  long n_header_lines;


  /* records true length of ref / alt alleles, which can be
   * truncated by limited buffer size of SNP datastructure
   */
  size_t ref_len;
  size_t alt_len;

  char qual[VCF_MAX_QUAL];
  char filter[VCF_MAX_FILTER];
  char info[VCF_MAX_FILTER];
  char format[VCF_MAX_FORMAT];

  long n_chrom;
  long max_chrom;
  Chromosome *chrom;

  /* used for reading lines */
  size_t buf_size;
  char *buf;
  
  /* could store lots of header info here */
} VCFInfo;



VCFInfo *vcf_info_new();
void vcf_info_free();

void vcf_read_header(gzFile vcf_fh, VCFInfo *vcf_info);

int vcf_read_line(gzFile vcf_fh, VCFInfo *vcf_info, SNP *snp);


#endif
