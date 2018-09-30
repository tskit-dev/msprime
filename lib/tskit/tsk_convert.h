#ifndef TSK_CONVERT_H
#define TSK_CONVERT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "tsk_genotypes.h"

/* TODO do we really need to expose this or would a simpler function be
 * more appropriate? Depends on how we use it at the Python level probably. */

typedef struct {
    size_t num_samples;
    size_t num_vcf_samples;
    unsigned int ploidy;
    char *genotypes;
    char *header;
    char *record;
    char *vcf_genotypes;
    size_t vcf_genotypes_size;
    size_t contig_id_size;
    size_t record_size;
    size_t num_sites;
    unsigned long contig_length;
    unsigned long *positions;
    tsk_vargen_t *vargen;
} tsk_vcf_converter_t;

int tsk_vcf_converter_alloc(tsk_vcf_converter_t *self,
        tsk_treeseq_t *tree_sequence, unsigned int ploidy, const char *chrom);
int tsk_vcf_converter_get_header(tsk_vcf_converter_t *self, char **header);
int tsk_vcf_converter_next(tsk_vcf_converter_t *self, char **record);
int tsk_vcf_converter_free(tsk_vcf_converter_t *self);
void tsk_vcf_converter_print_state(tsk_vcf_converter_t *self, FILE *out);


int tsk_convert_newick(tsk_tree_t *tree, tsk_id_t root, size_t precision, int flags,
        size_t buffer_size, char *buffer);

#ifdef __cplusplus
}
#endif
#endif
