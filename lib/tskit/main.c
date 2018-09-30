#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdarg.h>
#include <float.h>
#include <stdlib.h>

#include <regex.h>
#include "argtable3.h"

#include "tskit.h"

/* This file defines a crude CLI for tskit. It is intended for development
 * use only.
 */

typedef struct {
    int alphabet;
    double mutation_rate;
} mutation_params_t;

static void
fatal_error(const char *msg, ...)
{
    va_list argp;
    fprintf(stderr, "main:");
    va_start(argp, msg);
    vfprintf(stderr, msg, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}

static void
fatal_library_error(int err, const char *msg, ...)
{
    va_list argp;
    fprintf(stderr, "error:");
    va_start(argp, msg);
    vfprintf(stderr, msg, argp);
    va_end(argp);
    fprintf(stderr, ":%d:'%s'\n", err, tsk_strerror(err));
    exit(EXIT_FAILURE);
}

static void
load_tree_sequence(tsk_treeseq_t *ts, const char *filename)
{
    int ret = tsk_treeseq_load(ts, filename, 0);
    if (ret != 0) {
        fatal_library_error(ret, "Load error");
    }
}

static void
print_variants(tsk_treeseq_t *ts)
{
    int ret = 0;
    tsk_vargen_t vg;
    uint32_t j, k;
    tsk_variant_t* var;

    printf("variants (%d) \n", (int) tsk_treeseq_get_num_sites(ts));
    ret = tsk_vargen_alloc(&vg, ts, NULL, 0, 0);
    if (ret != 0) {
        fatal_library_error(ret, "tsk_vargen_alloc");
    }
    j = 0;
    while ((ret = tsk_vargen_next(&vg, &var)) == 1) {
        printf("%.2f\t", var->site->position);
        for (j = 0; j < var->num_alleles; j++) {
            for (k = 0; k < var->allele_lengths[j]; k++) {
                printf("%c", var->alleles[j][k]);
            }
            if (j < var->num_alleles - 1) {
                printf(",");
            }
        }
        printf("\t");
        for (k = 0; k < ts->num_samples; k++) {
            printf("%d\t", var->genotypes.u8[k]);
        }
        printf("\n");
    }
    if (ret != 0) {
        fatal_library_error(ret, "tsk_vargen_next");
    }
    tsk_vargen_free(&vg);
}

static void
print_haplotypes(tsk_treeseq_t *ts)
{
    int ret = 0;
    tsk_hapgen_t hg;
    uint32_t j;
    char *haplotype;

    printf("haplotypes \n");
    ret = tsk_hapgen_alloc(&hg, ts);
    if (ret != 0) {
        fatal_library_error(ret, "tsk_hapgen_alloc");
    }
    for (j = 0; j < ts->num_samples; j++) {
        ret = tsk_hapgen_get_haplotype(&hg, (tsk_id_t) j, &haplotype);
        if (ret < 0) {
            fatal_library_error(ret, "tsk_hapgen_get_haplotype");
        }
        printf("%d\t%s\n", j, haplotype);
    }
    tsk_hapgen_free(&hg);
}

static void
print_ld_matrix(tsk_treeseq_t *ts)
{
    int ret;
    size_t num_sites = tsk_treeseq_get_num_sites(ts);
    tsk_site_t sA, sB;
    double *r2 = malloc(num_sites * sizeof(double));
    size_t j, k, num_r2_values;
    tsk_ld_calc_t ld_calc;

    if (r2 == NULL) {
        fatal_error("no memory");
    }
    ret = tsk_ld_calc_alloc(&ld_calc, ts);
    printf("alloc: ret = %d\n", ret);
    if (ret != 0) {
        fatal_library_error(ret, "tsk_ld_calc_alloc");
    }
    tsk_ld_calc_print_state(&ld_calc, stdout);
    for (j = 0; j < num_sites; j++) {
        ret = tsk_ld_calc_get_r2_array(&ld_calc, j, TSK_DIR_FORWARD, num_sites,
                DBL_MAX, r2, &num_r2_values);
        if (ret != 0) {
            fatal_library_error(ret, "tsk_ld_calc_get_r2_array");
        }
        for (k = 0; k < num_r2_values; k++) {
            ret = tsk_treeseq_get_site(ts, j, &sA);
            if (ret != 0) {
                fatal_library_error(ret, "get_site");
            }
            ret = tsk_treeseq_get_site(ts, (j + k + 1), &sB);
            if (ret != 0) {
                fatal_library_error(ret, "get_site");
            }
            printf("%d\t%f\t%d\t%f\t%.3f\n",
                (int) sA.id, sA.position, (int) sB.id, sB.position, r2[k]);
        }
    }
    free(r2);
    ret = tsk_ld_calc_free(&ld_calc);
    if (ret != 0) {
        fatal_library_error(ret, "tsk_ld_calc_write_table");
    }
}

static void
print_stats(tsk_treeseq_t *ts)
{
    int ret = 0;
    uint32_t j;
    size_t num_samples = tsk_treeseq_get_num_samples(ts) / 2;
    tsk_id_t *sample = malloc(num_samples * sizeof(tsk_id_t));
    double pi;

    if (sample == NULL) {
        fatal_error("no memory");
    }
    for (j = 0; j < num_samples; j++) {
        sample[j] = (tsk_id_t) j;
    }
    ret = tsk_treeseq_get_pairwise_diversity(ts, sample, num_samples, &pi);
    if (ret != 0) {
        fatal_library_error(ret, "get_pairwise_diversity");
    }
    printf("pi = %f\n", pi);
    free(sample);
}

static void
print_vcf(tsk_treeseq_t *ts, unsigned int ploidy, const char *chrom, int verbose)
{
    int ret = 0;
    char *record = NULL;
    char *header = NULL;
    tsk_vcf_converter_t vc;

    ret = tsk_vcf_converter_alloc(&vc, ts, ploidy, chrom);
    if (ret != 0) {
        fatal_library_error(ret, "vcf alloc");
    }
    if (verbose > 0) {
        tsk_vcf_converter_print_state(&vc, stdout);
        printf("START VCF\n");
    }
    ret = tsk_vcf_converter_get_header(&vc, &header);
    if (ret != 0) {
        fatal_library_error(ret, "vcf get header");
    }
    printf("%s", header);
    while ((ret = tsk_vcf_converter_next(&vc, &record)) == 1) {
        printf("%s", record);
    }
    if (ret != 0) {
        fatal_library_error(ret, "vcf next");
    }
    tsk_vcf_converter_free(&vc);
}

static void
print_newick_trees(tsk_treeseq_t *ts)
{
    int ret = 0;
    char *newick = NULL;
    size_t precision = 8;
    size_t newick_buffer_size = (precision + 5) * tsk_treeseq_get_num_nodes(ts);
    tsk_tree_t tree;

    newick = malloc(newick_buffer_size);
    if (newick == NULL) {
        fatal_error("No memory\n");
    }

    ret = tsk_tree_alloc(&tree, ts, 0);
    if (ret != 0) {
        fatal_error("ERROR: %d: %s\n", ret, tsk_strerror(ret));
    }
    for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
        ret = tsk_convert_newick(&tree, tree.left_root, precision,
                0, newick_buffer_size, newick);
        if (ret != 0) {
            fatal_library_error(ret ,"newick");
        }
        printf("%d:\t%s\n", (int) tree.index, newick);
    }
    if (ret < 0) {
        fatal_error("ERROR: %d: %s\n", ret, tsk_strerror(ret));
    }
    tsk_tree_free(&tree);
    free(newick);
}

static void
print_tree_sequence(tsk_treeseq_t *ts, int verbose)
{
    int ret = 0;
    tsk_tree_t tree;

    tsk_treeseq_print_state(ts, stdout);
    if (verbose > 0) {
        printf("========================\n");
        printf("trees\n");
        printf("========================\n");
        ret = tsk_tree_alloc(&tree, ts, TSK_SAMPLE_COUNTS|TSK_SAMPLE_LISTS);
        if (ret != 0) {
            fatal_error("ERROR: %d: %s\n", ret, tsk_strerror(ret));
        }
        for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
            printf("-------------------------\n");
            printf("New tree: %d: %f (%d)\n", (int) tree.index,
                    tree.right - tree.left, (int) tree.num_nodes);
            printf("-------------------------\n");
            tsk_tree_print_state(&tree, stdout);
        }
        if (ret < 0) {
            fatal_error("ERROR: %d: %s\n", ret, tsk_strerror(ret));
        }
        tsk_tree_free(&tree);
    }
}

static void
run_ld(const char *filename, int TSK_UNUSED(verbose))
{
    tsk_treeseq_t ts;

    load_tree_sequence(&ts, filename);
    print_ld_matrix(&ts);
    tsk_treeseq_free(&ts);
}

static void
run_haplotypes(const char *filename, int TSK_UNUSED(verbose))
{
    tsk_treeseq_t ts;

    load_tree_sequence(&ts, filename);
    print_haplotypes(&ts);
    tsk_treeseq_free(&ts);
}

static void
run_variants(const char *filename, int TSK_UNUSED(verbose))
{
    tsk_treeseq_t ts;

    load_tree_sequence(&ts, filename);
    print_variants(&ts);
    tsk_treeseq_free(&ts);
}

static void
run_vcf(const char *filename, int verbose, int ploidy, const char *chrom)
{
    tsk_treeseq_t ts;

    load_tree_sequence(&ts, filename);
    print_vcf(&ts, (unsigned int) ploidy, chrom, verbose);
    tsk_treeseq_free(&ts);
}

static void
run_print(const char *filename, int verbose)
{
    tsk_treeseq_t ts;

    load_tree_sequence(&ts, filename);
    print_tree_sequence(&ts, verbose);
    tsk_treeseq_free(&ts);
}

static void
run_newick(const char *filename, int TSK_UNUSED(verbose))
{
    tsk_treeseq_t ts;

    load_tree_sequence(&ts, filename);
    print_newick_trees(&ts);
    tsk_treeseq_free(&ts);
}

static void
run_stats(const char *filename, int TSK_UNUSED(verbose))
{
    tsk_treeseq_t ts;

    load_tree_sequence(&ts, filename);
    print_stats(&ts);
    tsk_treeseq_free(&ts);
}

static void
run_simplify(const char *input_filename, const char *output_filename, size_t num_samples,
        bool filter_sites, int verbose)
{
    tsk_treeseq_t ts, subset;
    tsk_id_t *samples;
    int flags = 0;
    int ret;

    if (filter_sites) {
        flags |= TSK_FILTER_SITES;
    }

    load_tree_sequence(&ts, input_filename);
    if (verbose > 0) {
        printf(">>>>>>>>\nINPUT:\n>>>>>>>>\n");
        tsk_treeseq_print_state(&ts, stdout);
    }
    if (num_samples == 0) {
        num_samples = tsk_treeseq_get_num_samples(&ts);
    } else {
        num_samples = TSK_MIN(num_samples, tsk_treeseq_get_num_samples(&ts));
    }
    ret = tsk_treeseq_get_samples(&ts, &samples);
    if (ret != 0) {
        fatal_library_error(ret, "get_samples");
    }
    ret = tsk_treeseq_simplify(&ts, samples, num_samples, flags, &subset, NULL);
    if (ret != 0) {
        fatal_library_error(ret, "Subset error");
    }
    ret = tsk_treeseq_dump(&subset, output_filename, 0);
    if (ret != 0) {
        fatal_library_error(ret, "Write error");
    }
    if (verbose > 0) {
        printf(">>>>>>>>\nOUTPUT:\n>>>>>>>>\n");
        tsk_treeseq_print_state(&subset, stdout);
    }
    tsk_treeseq_free(&ts);
    tsk_treeseq_free(&subset);
}

int
main(int argc, char** argv)
{
    /* SYNTAX 1: simplify [-vi] [-s] <input-file> <output-file> */
    struct arg_rex *cmd1 = arg_rex1(NULL, NULL, "simplify", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose1 = arg_lit0("v", "verbose", NULL);
    struct arg_int *num_samples1 = arg_int0("s", "sample-size", "<sample-size>",
            "Number of samples to keep in the simplified tree sequence.");
    struct arg_lit *filter_sites1 = arg_lit0("i",
            "filter-invariant-sites", "<filter-invariant-sites>");
    struct arg_file *infiles1 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_file *outfiles1 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end1 = arg_end(20);
    void* argtable1[] = {cmd1, verbose1, filter_sites1, num_samples1,
        infiles1, outfiles1, end1};
    int nerrors1;

    /* SYNTAX 2: ld [-v] <input-file> */
    struct arg_rex *cmd2 = arg_rex1(NULL, NULL, "ld", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose2 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles2 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end2 = arg_end(20);
    void* argtable2[] = {cmd2, verbose2, infiles2, end2};
    int nerrors2;

    /* SYNTAX 3: haplotypes [-v] <input-file> */
    struct arg_rex *cmd3 = arg_rex1(NULL, NULL, "haplotypes", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose3 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles3 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end3 = arg_end(20);
    void* argtable3[] = {cmd3, verbose3, infiles3, end3};
    int nerrors3;

    /* SYNTAX 4: variants [-v] <input-file> */
    struct arg_rex *cmd4 = arg_rex1(NULL, NULL, "variants", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose4 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles4 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end4 = arg_end(20);
    void* argtable4[] = {cmd4, verbose4, infiles4, end4};
    int nerrors4;

    /* SYNTAX 5: vcf [-v] <input-file> */
    struct arg_rex *cmd5 = arg_rex1(NULL, NULL, "vcf", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose5 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles5 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_int *ploidy5 = arg_int0("p", "ploidy", "<ploidy>",
            "Ploidy level of the VCF");
    struct arg_str *chrom5 = arg_str0("c", "chrom", "<chrom>",
            "Value for the CHROM column in the VCF (default='1')");
    struct arg_end *end5 = arg_end(20);
    void* argtable5[] = {cmd5, verbose5, infiles5, ploidy5, chrom5, end5};
    int nerrors5;

    /* SYNTAX 6: print  [-v] <input-file> */
    struct arg_rex *cmd6 = arg_rex1(NULL, NULL, "print", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose6 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles6 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end6 = arg_end(20);
    void* argtable6[] = {cmd6, verbose6, infiles6, end6};
    int nerrors6;

    /* SYNTAX 7: newick [-v] <input-file> */
    struct arg_rex *cmd7 = arg_rex1(NULL, NULL, "newick", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose7 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles7 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end7 = arg_end(20);
    void* argtable7[] = {cmd7, verbose7, infiles7, end7};
    int nerrors7;

    /* SYNTAX 8: stats [-v] <input-file> */
    struct arg_rex *cmd8 = arg_rex1(NULL, NULL, "stats", NULL, REG_ICASE, NULL);
    struct arg_lit *verbose8 = arg_lit0("v", "verbose", NULL);
    struct arg_file *infiles8 = arg_file1(NULL, NULL, NULL, NULL);
    struct arg_end *end8 = arg_end(20);
    void* argtable8[] = {cmd8, verbose8, infiles8, end8};
    int nerrors8;

    int exitcode = EXIT_SUCCESS;
    const char *progname = "main";

    /* Set defaults */
    ploidy5->ival[0] = 1;
    chrom5->sval[0] = "1";
    num_samples1->ival[0] = 0;

    nerrors1 = arg_parse(argc, argv, argtable1);
    nerrors2 = arg_parse(argc, argv, argtable2);
    nerrors3 = arg_parse(argc, argv, argtable3);
    nerrors4 = arg_parse(argc, argv, argtable4);
    nerrors5 = arg_parse(argc, argv, argtable5);
    nerrors6 = arg_parse(argc, argv, argtable6);
    nerrors7 = arg_parse(argc, argv, argtable7);
    nerrors8 = arg_parse(argc, argv, argtable8);

    if (nerrors1 == 0) {
        run_simplify(infiles1->filename[0], outfiles1->filename[0],
                (size_t) num_samples1->ival[0], (bool) filter_sites1->count,
                verbose1->count);
    } else if (nerrors2 == 0) {
        run_ld(infiles2->filename[0], verbose2->count);
    } else if (nerrors3 == 0) {
        run_haplotypes(infiles3->filename[0], verbose3->count);
    } else if (nerrors4 == 0) {
        run_variants(infiles4->filename[0], verbose4->count);
    } else if (nerrors5 == 0) {
        run_vcf(infiles5->filename[0], verbose5->count, ploidy5->ival[0], chrom5->sval[0]);
    } else if (nerrors6 == 0) {
        run_print(infiles6->filename[0], verbose6->count);
    } else if (nerrors7 == 0) {
        run_newick(infiles7->filename[0], verbose7->count);
    } else if (nerrors8 == 0) {
        run_stats(infiles8->filename[0], verbose8->count);
    } else {
        /* We get here if the command line matched none of the possible syntaxes */
        if (cmd1->count > 0) {
            arg_print_errors(stdout, end1, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable1, "\n");
        } else if (cmd2->count > 0) {
            arg_print_errors(stdout, end2, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable2, "\n");
        } else if (cmd3->count > 0) {
            arg_print_errors(stdout, end3, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable3, "\n");
        } else if (cmd4->count > 0) {
            arg_print_errors(stdout, end4, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable4, "\n");
        } else if (cmd5->count > 0) {
            arg_print_errors(stdout, end5, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable5, "\n");
        } else if (cmd6->count > 0) {
            arg_print_errors(stdout, end6, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable6, "\n");
        } else if (cmd7->count > 0) {
            arg_print_errors(stdout, end7, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable7, "\n");
        } else if (cmd8->count > 0) {
            arg_print_errors(stdout, end8, progname);
            printf("usage: %s ", progname);
            arg_print_syntax(stdout, argtable8, "\n");
        } else {
            /* no correct cmd literals were given, so we cant presume which syntax was intended */
            printf("%s: missing command.\n",progname);
            printf("usage 1: %s ", progname);  arg_print_syntax(stdout, argtable1, "\n");
            printf("usage 2: %s ", progname);  arg_print_syntax(stdout, argtable2, "\n");
            printf("usage 3: %s ", progname);  arg_print_syntax(stdout, argtable3, "\n");
            printf("usage 4: %s ", progname);  arg_print_syntax(stdout, argtable4, "\n");
            printf("usage 5: %s ", progname);  arg_print_syntax(stdout, argtable5, "\n");
            printf("usage 6: %s ", progname);  arg_print_syntax(stdout, argtable6, "\n");
            printf("usage 7: %s ", progname);  arg_print_syntax(stdout, argtable7, "\n");
            printf("usage 8: %s ", progname);  arg_print_syntax(stdout, argtable8, "\n");
        }
    }

    arg_freetable(argtable1, sizeof(argtable1) / sizeof(argtable1[0]));
    arg_freetable(argtable2, sizeof(argtable2) / sizeof(argtable2[0]));
    arg_freetable(argtable3, sizeof(argtable3) / sizeof(argtable3[0]));
    arg_freetable(argtable4, sizeof(argtable4) / sizeof(argtable4[0]));
    arg_freetable(argtable5, sizeof(argtable5) / sizeof(argtable5[0]));
    arg_freetable(argtable6, sizeof(argtable6) / sizeof(argtable6[0]));
    arg_freetable(argtable7, sizeof(argtable7) / sizeof(argtable7[0]));
    arg_freetable(argtable8, sizeof(argtable8) / sizeof(argtable8[0]));

    return exitcode;
}
