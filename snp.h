#ifndef __SNP_H__
#define __SNP_H__


#define SNP_MAX_ALLELE 1024
#define SNP_MAX_CHROM 1024
#define SNP_MAX_NAME 1024


typedef struct {
  char name[SNP_MAX_NAME];
  char chrom_name[SNP_MAX_CHROM];
  long pos;
  char allele1[SNP_MAX_ALLELE];
  char allele2[SNP_MAX_ALLELE];
  
  char has_geno_probs;
  char has_haplotypes;
  float *geno_probs;
  char *haplotypes;
} SNP;



#endif
