
#include <zlib.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <stdint.h>

#include "vcf.h"
#include "snp.h"
#include "util.h"
#include "memutil.h"

typedef struct  {
  gzFile gzf;
  Chromosome *cur_chrom;
  VCFInfo *vcf;
  char is_done;
  SNP cur_snp;
} FileInfo;



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
Chromosome *chrom_table_intersect(VCFInfo **vcf, int n_vcf, int *n_intersect) {
  int i, j, k, n_chrom;
  int *counts;
  Chromosome *intersect;

  if(n_vcf < 1) {
    my_err("expected at least 1 vcf");
  }

  counts = my_malloc(sizeof(int) * vcf[0]->n_chrom);
  for(i = 0; i < vcf[0]->n_chrom; i++) {
    counts[i] = 1;
  }

  *n_intersect = 0;

  /* count how many other VCFs each chrom occurs in */
  for(i = 0; i < vcf[0]->n_chrom; i++) {
    for(j = 1; j < n_vcf; j++) {
      for(k = 0; k < vcf[j]->n_chrom; k++) {
	if(strcmp(vcf[0]->chrom[i].name, vcf[j]->chrom[k].name) == 0) {
	  counts[i] += 1;
	  if(counts[i] == n_vcf) {
	    /* this chrom found in all VCFs */
	    *n_intersect += 1;
	  }
	  break;
	}
      }
    }
  }

  intersect = my_malloc(sizeof(Chromosome) * *n_intersect);
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


  my_free(counts);
  
  return intersect;
}




void init_file_info(int n_vcf, char **vcf_filenames) {
  FileInfo *f_info;
  int i, ret;

  f_info = my_malloc(sizeof(FileInfo *) * n_vcf);

  for(i = 0; i < n_vcf; i++) {
    f_info[i].vcf = vcf_info_new();
    fprintf(stderr, "reading VCF header from %s\n", vcf_filenames[i]);
    f_info[i].gzf = util_must_gzopen(vcf_filenames[i], "rb");
    vcf_read_header(f_info[i].gzf[i], vcf[i]);
    fprintf(stderr, "  VCF header lines: %ld\n", f_info[i].vcf->n_header_lines);

    f_info[i].is_done = FALSE;

    
    /* initialize memory for current SNPs */
    f_info[i].cur_snp.geno_probs =
      my_malloc(sizeof(float) * file_info[i].vcf->n_geno_prob_col);
    f_info[i].cur_snp.haplotypes =
      my_malloc(sizeof(char) * file_info[i].vcf->n_haplo_col);
    
    /* read first SNP from all files */
    ret = vcf_read_line(f_info[i].gzf, f_info[i].vcf, &f_info[i].cur_snp);
    if(ret == -1) {
      my_err("file %s contains no SNPs\n", vcf_filename[i]);
    }
  }

  return f_info;
}


void free_file_info(FileInfo *f_info, int n) {
  for(i = 0; i < n; i++) {
    vcf_info_free(f_info[i].vcf);
    
    gzclose(f_info.gzf[i]);

    if(f_info[i].cur_snps.haplotypes) {
      my_free(f_info[i].cur_snps.haplotypes);
    }
    if(f_info[i].cur_snps.geno_probs) {
      my_free(f_info[i]..geno_probs);
    }
  }
  my_free(f_info);
}



void merge_vcf(int n_vcf, char **vcf_filenames) {
  FileInfo *f_info;
  int n_done;

  f_info = init_file_info(n_vcf, vcf_filenames):
  
  /* find chromosomes that are present in ALL VCFs */
  chrom_tab = chrom_table_intersect(f_info, &n_chrom);
  n_done = 0;

  
  fprintf(stderr, "parsing files\n");
  while(n_done < n_vcf) {
  
  while(vcf_read_line(gzf[0], f_info.vcf[0], &cur_snps[0]) != -1) {
    /** TODO **/
    break;
  }
  
  fprintf(stderr, "\n");

  free_file_info(f_info, n_vcf);
  
  for(i = 0; i < n_vcf; i++) {
    vcf_info_free(vcf[i]);
    gzclose(gzf[i]);

    if(cur_snps[i].haplotypes) {
      my_free(cur_snps[i].haplotypes);
    }
    if(cur_snps[i].geno_probs) {
      my_free(cur_snps[i].geno_probs);
    }
  }
  for(i = 0; i < n_chrom; i++) {
    my_free(chrom_tab[i].name);
    my_free(chrom_tab[i].assembly);
  }
  my_free(chrom_tab);
  my_free(vcf);
  my_free(cur_snps);
  my_free(is_done);
  my_free(gzf);
  
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
