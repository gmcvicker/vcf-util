
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
Chromosome *chrom_table_intersect(FileInfo *f_info, int n_vcf, int *n_intersect) {
  int i, j, k, n_chrom;
  int *counts;
  Chromosome *intersect;

  if(n_vcf < 1) {
    my_err("expected at least 1 vcf");
  }

  counts = my_malloc(sizeof(int) * f_info[0].vcf->n_chrom);
  for(i = 0; i < f_info[0].vcf->n_chrom; i++) {
    counts[i] = 1;
  }

  *n_intersect = 0;

  /* count how many other VCFs each chrom occurs in */
  for(i = 0; i < f_info[0].vcf->n_chrom; i++) {
    for(j = 1; j < n_vcf; j++) {
      for(k = 0; k < f_info[j].vcf->n_chrom; k++) {
	if(strcmp(f_info[0].vcf->chrom[i].name,
		  f_info[j].vcf->chrom[k].name) == 0) {
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
  for(i = 0; i < f_info[0].vcf->n_chrom; i++) {
    if(counts[i] == n_vcf) {
      /* copy chrom found in all VCFs */
      intersect[j].id = j;
      intersect[j].name = util_str_dup(f_info[0].vcf->chrom[i].name);
      intersect[j].assembly = util_str_dup(f_info[0].vcf->chrom[i].assembly);
      intersect[j].len = f_info[0].vcf->chrom[i].len;
      j += 1;
    } else {
      fprintf(stderr, "skipping chromosome %s because only found "
	      "in %d/%d VCFs\n", f_info[0].vcf->chrom[i].name, counts[i], n_vcf);
    }
  }


  my_free(counts);
  
  return intersect;
}




FileInfo *init_file_info(int n_vcf, char **vcf_filenames) {
  FileInfo *f_info;
  int i, ret;

  f_info = my_malloc(sizeof(FileInfo) * n_vcf);

  for(i = 0; i < n_vcf; i++) {
    f_info[i].vcf = vcf_info_new();
    fprintf(stderr, "reading VCF header from %s\n", vcf_filenames[i]);
    f_info[i].gzf = util_must_gzopen(vcf_filenames[i], "rb");
    vcf_read_header(f_info[i].gzf, f_info[i].vcf);
    fprintf(stderr, "  VCF header lines: %ld\n", f_info[i].vcf->n_header_lines);

    f_info[i].is_done = FALSE;

    /* initialize memory for current SNPs */
    f_info[i].cur_snp.geno_probs =
      my_malloc(sizeof(float) * f_info[i].vcf->n_geno_prob_col);
    f_info[i].cur_snp.haplotypes =
      my_malloc(sizeof(char) * f_info[i].vcf->n_haplo_col);

    f_info[i].cur_snp.has_geno_probs = FALSE;
    f_info[i].cur_snp.has_haplotypes = FALSE;

    f_info[i].cur_chrom = NULL;
  }

  return f_info;
}


void free_file_info(FileInfo *f_info, int n) {
  int i;
  
  for(i = 0; i < n; i++) {
    vcf_info_free(f_info[i].vcf);
    
    gzclose(f_info[i].gzf);

    if(f_info[i].cur_snp.haplotypes) {
      my_free(f_info[i].cur_snp.haplotypes);
    }
    if(f_info[i].cur_snp.geno_probs) {
      my_free(f_info[i].cur_snp.geno_probs);
    }
  }
  my_free(f_info);
}



void set_cur_chrom(FileInfo *f_info, Chromosome *chrom_tab, int n_chrom) {
  int i;
  
  if(f_info->cur_chrom != NULL) {
    if(strcmp(f_info->cur_snp.chrom_name, f_info->cur_chrom->name) == 0) {
      /* current chromosome matches name in SNP */
      return;
    }
  }
  for(i = 0; i < n_chrom; i++) {
    if(strcmp(f_info->cur_snp.chrom_name, chrom_tab[i].name) == 0) {
      f_info->cur_chrom = &chrom_tab[i];
    }
  }
}



/**
 * compare (chrom, pos) from SNPs from two files, return -1 if SNP 1 is lower,
 * return +1 if SNP 2 is lower, return 0 if they have equal position.
 * 
 */
int f_cmp(FileInfo *f1, FileInfo *f2) {
  if(f1->cur_chrom->id < f2->cur_chrom->id) {
    return -1;
  }
  if(f1->cur_chrom->id > f2->cur_chrom->id) {
    return 1;
  }
  if(f1->cur_snp.pos < f2->cur_snp.pos) {
    return -1;
  }
  if(f1->cur_snp.pos > f2->cur_snp.pos) {
    return 1;
  }
  return 0;
}


/**
 * Finds file(s) with current SNP(s) with lowest coordinates.
 * Sets values in is_lowest and lowest arrays:
 * - is_lowest is of length n_vcf and has TRUE or FALSE flags.
 * - lowest has *n_lowest values which are indices pointing
 *   to elements in f_info array.
 */
void find_lowest(FileInfo *f_info, int n_vcf,
		 int *is_lowest, int *lowest, int *n_lowest) {
  int i, j;
  FileInfo *cur_lowest;

  *n_lowest = 0;
  
  for(i = 0; i < n_vcf; i++) {
    is_lowest[i] = FALSE;
  }
  for(i = 0; i < n_vcf; i++) {
    if(f_info[i].is_done) {
      /* at end of this file */
      continue;
    }
    if(*n_lowest == 0) {
      /* first file that is not at end */
      lowest[0] = i;
      *n_lowest = 1;
      cur_lowest = &f_info[i];
    } else {
      /* is this file's SNP lower than current lowest? */
      int c = f_cmp(&f_info[i], cur_lowest);
      if(c < 0) {
	/* new lowest SNP */
	cur_lowest = &f_info[i];
	lowest[0] = i;
	*n_lowest = 1;
      } else if(c == 0) {
	  /* another SNP that matches lowest */
	  lowest[*n_lowest] = i;
	  *n_lowest += 1;
      }
    }
  }

  for(i = 0; i < *n_lowest; i++) {
    is_lowest[lowest[i]] = TRUE;
  }
}


void write_output(FILE *f, FileInfo *f_info, int n_vcf, int *is_lowest,
		  int *lowest, int write_geno_probs, int write_haplotypes) {
  SNP *s;
  char *format_str, *filter_str;
  int qual;
  

  /* TODO: BUILD NEW INFO field */
  /* TODO: NOT SURE WHAT TO DO ABOUT QUAL, FILTER */
  /* TODO: check that alleles match! */

  if(write_geno_probs && write_haplotypes) {
    format_str = "GL;GT";
  }
  else if(write_haplotypes) {
    format_str = "GT";
  }
  else if(write_geno_probs) {
    format_str = "GL";
  }
  
  qual = 100;
  filter_str = "PASS";
  
  /* obtain SNP info from first of SNPs that is in group of lowest SNPs */
  s = &f_info[lowest[0]].cur_snp;
  fprintf(f, "%s\t%ld\t%s\t%s\t%s\t%d\t%s\t%s", s->chrom_name, s->pos, s->name,
	  s->allele1, s->allele2, qual, filter_str, format_str);

  /* TODO: write out genotype information for every file... */
  
  fprintf(f, "\n");
  

}


void merge_vcf(int n_vcf, char **vcf_filenames) {
  FileInfo *f_info;
  int n_done, n_chrom, i, *is_lowest, *lowest, n_lowest;
  int ret, use_geno_probs, use_haplotypes;
  Chromosome *chrom_tab;

  f_info = init_file_info(n_vcf, vcf_filenames);
  
  /* find chromosomes that are present in ALL VCFs */
  chrom_tab = chrom_table_intersect(f_info, n_vcf, &n_chrom);
  n_done = 0;
  is_lowest = my_malloc(sizeof(int) * n_vcf);
  lowest = my_malloc(sizeof(int) * n_vcf);

  /* only use genotypes and haplotypes if they are present in ALL files */
  use_geno_probs = TRUE;
  use_haplotypes = TRUE;
  
  /* read first SNP from all files */
  for(i = 0; i < n_vcf; i++) {
    ret = vcf_read_line(f_info[i].gzf, f_info[i].vcf, &f_info[i].cur_snp);
    if(ret == -1) {
      /* file is over */
      n_done += 1;
      f_info[i].is_done = TRUE;
      my_warn("file %s contains no SNPs\n", vcf_filenames[i]);
      f_info[i].cur_chrom = NULL;
    } else {
      set_cur_chrom(&f_info[i], chrom_tab, n_chrom);

      if(!f_info[i].cur_snp.has_geno_probs) {
	if(use_geno_probs) {
	  fprintf(stderr, "Not using genotype likelihoods (GL) because "
		  "not present in file %s\n", vcf_filenames[i]);
	}
	use_geno_probs = FALSE;
      }
      if(!f_info[i].cur_snp.has_haplotypes) {
	if(use_haplotypes) {
	  fprintf(stderr, "Not using genotypes (GT) because "
		  "not present in file %s\n", vcf_filenames[i]);
	}
	use_haplotypes = FALSE;
      }
    }
  }
  
  fprintf(stderr, "parsing files\n");

  while(n_done < n_vcf) {
    /* find SNP(s) with lowest (chrom, pos) */
    find_lowest(f_info, n_vcf, is_lowest, lowest, &n_lowest);

    /* merge counts and write line for these SNPs */
    write_output(stdout, f_info, n_vcf, is_lowest, lowest,
		 use_geno_probs, use_haplotypes);
    
    /* advance files with lowest SNPs */
    for(i = 0; i < n_vcf; i++) {
      if(!f_info[i].is_done && is_lowest[i]) {
	if(vcf_read_line(f_info[i].gzf, f_info[i].vcf, &f_info[i].cur_snp) == -1) {
	  /* have reached end of this file */
	  n_done += 1;
	  f_info[i].is_done = TRUE;
	}
      }
    }
  }
  
  fprintf(stderr, "done!\n");

  free_file_info(f_info, n_vcf);
  for(i = 0; i < n_chrom; i++) {
    my_free(chrom_tab[i].name);
    my_free(chrom_tab[i].assembly);
  }
  my_free(chrom_tab);
  my_free(is_lowest);
  
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
