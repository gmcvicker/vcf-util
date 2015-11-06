
#include <zlib.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <stdint.h>

#include "vcf.h"
#include "snp.h"
#include "util.h"
#include "memutil.h"



void usage(char **argv) {
  fprintf(stderr, "\nusage: %s VCF1 VCF2 ... > MERGED_VCF\n"
	  "\n"
	  "Description:\n"
	  "  This program merges VCF files. Input VCF files must be sorted\n"
	  "\n", argv[0]);
}




/**
 * Find subset of chromosomes that are present in all VCFs
 */
Chromosome *chrom_table_intersect(VCFInfo **vcf, int n_vcf) {
  int i, j, k, n_chrom, n_intersect;
  int *counts;
  Chromosome *intersect;

  if(n_vcf < 1) {
    my_err("expected at least 1 vcf");
  }

  counts = my_malloc(sizeof(int) * vcf[0]->n_chrom);
  for(i = 0; i < vcf[0]->n_chrom; i++) {
    counts[i] = 1;
  }

  n_intersect = 0;

  /* count how many other VCFs each chrom occurs in */
  for(i = 0; i < vcf[0]->n_chrom; i++) {
    for(j = 1; j < n_vcf; j++) {
      for(k = 0; k < vcf[j]->n_chrom; k++) {
	if(strcmp(vcf[0]->chrom[i].name, vcf[j]->chrom[k].name) == 0) {
	  counts[i] += 1;
	  if(counts[i] == n_vcf) {
	    /* this chrom found in all VCFs */
	    n_intersect += 1;
	  }
	  break;
	}
      }
    }
  }

  intersect = my_malloc(sizeof(Chromosome) * n_intersect);
  j = 0;
  for(i = 0; i < vcf[0]->n_chrom; i++) {
    if(counts[i] == n_vcf) {
      /* copy chrom found in all VCFs */
      intersect[j].id = j;
      intersect[j].name = util_str_dup(vcf[0]->chrom[i].name);
      intersect[j].assembly = util_str_dup(vcf[0]->chrom[i].assembly);
      intersect[j].len = vcf[0]->chrom[i].len;
      j += 1;
    } else {
      fprintf(stderr, "skipping chromosome %s because only found "
	      "in %d/%d VCFs\n", vcf[0]->chrom[i].name, counts[i], n_vcf);
    }
  }
  
  return intersect;
}





void merge_vcf(int n_vcf, char **vcf_filenames) {
  VCFInfo **vcf;
  SNP *cur_snps;
  SNP snp;
  int i;
  gzFile *gzf;
  char *cur_chrom, *is_done;
  Chromosome *chrom_tab;
  

  vcf = my_malloc(sizeof(VCFInfo *) * n_vcf);
  cur_snps = my_malloc(sizeof(SNP) * n_vcf);
  gzf = my_malloc(sizeof(SNP) * n_vcf);
  is_done = my_malloc(sizeof(char) * n_vcf);
  
  cur_chrom = NULL;

  fprintf(stderr, "n_vcf_files=%d\n", n_vcf);
  
  /* open all files and read header/sample information from them */
  for(i = 0; i < n_vcf; i++) {
    vcf[i] = vcf_info_new();

    fprintf(stderr, "reading VCF header from %s\n", vcf_filenames[i]);
    gzf[i] = util_must_gzopen(vcf_filenames[i], "rb");
    vcf_read_header(gzf[i], vcf[i]);
    fprintf(stderr, "  VCF header lines: %ld\n", vcf[i]->n_header_lines);
    
    is_done[i] = FALSE;
    
    /* initialize memory for current SNPs */
    cur_snps[i].geno_probs = my_malloc(sizeof(float) * vcf[i]->n_geno_prob_col);
    cur_snps[i].haplotypes = my_malloc(sizeof(char) * vcf[i]->n_haplo_col);
    
    /** TODO: need to do renaming of samples when there are duplicate
     ** sample names (postfix -1, -2, etc.) 
     **/
  }

  /* find chromosomes that are present in ALL VCFs */
  chrom_tab = chrom_table_intersect(vcf, n_vcf);

  
  /* read first SNP from all files, check that chromosomes are the same */
  int ret;
  for(i = 0; i < n_vcf; i++) {
    ret = vcf_read_line(gzf[i], vcf[i], &cur_snps[i]);
    if(ret == -1) {
      is_done[i] = TRUE;
    } else {
      if(cur_chrom == NULL) {
	cur_chrom = util_str_dup(cur_snps[i].chrom);
      } else {
	if(strcmp(cur_chrom, cur_snps[i].chrom) != 0) {
	  my_err("%s:%d: VCFs do not all start with same chromosome. "
		 "%s != %s", __FILE__, __LINE__, cur_chrom, cur_snps[i].chrom);
	}
      }
    }
  }

  fprintf(stderr, "parsing files\n");

  /** TODO! **/
  while(vcf_read_line(gzf, vcf[0], &snp) != -1) {
    /** TODO **/
  }
  
  fprintf(stderr, "\n");
  
  gzclose(gzf);

  for(i = 0; i < n_vcf; i++) {
    vcf_info_free(vcf[i]);
  }
  my_free(vcf);
  
}




int main(int argc, char **argv) {
  int n_vcf;
  char **vcf_filenames;

  n_vcf = argc-1;

  if(n_vcf < 2) {
    usage(argv);
    exit(255);
  }

  vcf_filenames = &argv[1];
  
  merge_vcf(n_vcf, vcf_filenames);

  fprintf(stderr, "done\n");
  
  return 0;
}
