
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





void parse_vcf(int n_vcf_files, char **vcf_filenames) {
  VCFInfo **vcf;
  SNP *cur_snps;
  FileInfo file_info;
  SNP snp;
  int n_vcf;
  gzFile gzf;
  float *geno_probs;
  char *haplotypes;

  Chromosome *chrom;

  vcf = my_malloc(sizeof(VCFInfo *) * n_vcf);
  cur_snps = my_malloc(sizeof(SNP) * n_vcf);

  for(i = 0; i < n_vcf; i++) {
    vcf[i] = vcf_info_new();

    /* read information from header lines */
    fprintf(stderr, "reading VCF header from %s\n", vcf_filenames[i]);
    gzf = util_must_gzopen(vcf_filenames[i], "rb");
    vcf_read_header(gzf, vcf[i]);
    fprintf(stderr, "  VCF header lines: %ld\n", vcf[i]->n_header_lines);

    /** TODO: need to do renaming of samples when there are duplicate
     ** sample names (postfix -1, -2, etc.) 
     **/
    
    row = 0;

    fprintf(stderr, "parsing files\n");
    
    while(vcf_read_line(gzf, vcf, &snp, geno_probs, haplotypes) != -1) {

      if(geno_probs) {
	write_h5matrix_row(gprob_info, row, geno_probs);
      }
      if(haplotypes) {
	write_h5matrix_row(haplotype_info, row, haplotypes);
      }

      /*  set snp_index element at this chromosome position
       * to point to row in matrices / SNP table 
       */
      if(snp.pos > chrom->len || snp.pos < 1) {
	my_err("%s:%d: SNP position (%ld) is outside of "
	       "chromomosome %s range:1-%ld", __FILE__, __LINE__,
	       snp.pos, chrom->len);
      }

      if(snp_index) {
	/* set value in snp_index array to point to current row */
	snp_index[snp.pos-1] = row;
      }

      /* append row to SNP table */
      if(snp_tab) {
	snp_tab_append_row(snp_tab, &snp);
      }
      
      row++;
      if((row % 1000) == 0) {
	fprintf(stderr, ".");
      }
    }
    fprintf(stderr, "\n");

    if(row != file_info.n_row) {
      my_warn("%s:%d: expected %ld data rows, but only read %ld\n",
	      __FILE__, __LINE__, file_info.n_row, row);
    }

    gzclose(gzf);
  }

  vcf_info_free(vcf);
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
