#ifndef __TSK_CORE_H__
#define __TSK_CORE_H__

#include <stdbool.h>

#ifdef __GNUC__
	#define TSK_WARN_UNUSED __attribute__ ((warn_unused_result))
	#define TSK_UNUSED(x) TSK_UNUSED_ ## x __attribute__((__unused__))
#else
	#define TSK_WARN_UNUSED
	#define TSK_UNUSED(x) TSK_UNUSED_ ## x
    /* Don't bother with restrict for MSVC */
    #define restrict
#endif

/* This sets up TSK_DBL_DECIMAL_DIG, which can then be used as a
 * precision specifier when writing out doubles, if you want sufficient
 * decimal digits to be written to guarantee a lossless round-trip
 * after being read back in.  Usage:
 *
 *     printf("%.*g", TSK_DBL_DECIMAL_DIG, foo);
 *
 * See https://stackoverflow.com/a/19897395/2752221
 */
#ifdef DBL_DECIMAL_DIG
#define TSK_DBL_DECIMAL_DIG (DBL_DECIMAL_DIG)
#else
#define TSK_DBL_DECIMAL_DIG (DBL_DIG + 3)
#endif


/* Node flags */
#define TSK_NODE_IS_SAMPLE 1u

/* The root node indicator */
#define TSK_NULL_NODE (-1)
/* Indicates the that the population ID has not been set. */
#define TSK_NULL_POPULATION (-1)
/* There is no parent for a given mutation */
#define TSK_NULL_MUTATION (-1)
/* Indicates that no individual has been set */
#define TSK_NULL_INDIVIDUAL (-1)

/* Flags for simplify() */
#define TSK_FILTER_SITES                 (1 << 0)
#define TSK_REDUCE_TO_SITE_TOPOLOGY      (1 << 1)
#define TSK_FILTER_POPULATIONS           (1 << 2)
#define TSK_FILTER_INDIVIDUALS           (1 << 3)

/* Flags for check_integrity */
#define TSK_CHECK_OFFSETS                (1 << 0)
#define TSK_CHECK_EDGE_ORDERING          (1 << 1)
#define TSK_CHECK_SITE_ORDERING          (1 << 2)
#define TSK_CHECK_SITE_DUPLICATES        (1 << 3)
#define TSK_CHECK_MUTATION_ORDERING      (1 << 4)
#define TSK_CHECK_INDEXES                (1 << 5)
#define TSK_CHECK_ALL                    \
    (TSK_CHECK_OFFSETS | TSK_CHECK_EDGE_ORDERING | TSK_CHECK_SITE_ORDERING | \
     TSK_CHECK_SITE_DUPLICATES | TSK_CHECK_MUTATION_ORDERING | TSK_CHECK_INDEXES)

/* Flags for dump tables */
#define TSK_ALLOC_TABLES 1

/* Flags for load tables */
#define TSK_BUILD_INDEXES 1

/* Generic debug flag shared across all calls. Uses
 * the top bit to avoid clashes with other flags. */
#define TSK_DEBUG                       (1 << 31)

#define TSK_LOAD_EXTENDED_CHECKS  1

#define TSK_FILE_FORMAT_NAME          "tskit.trees"
#define TSK_FILE_FORMAT_NAME_LENGTH   11
#define TSK_FILE_FORMAT_VERSION_MAJOR 12
#define TSK_FILE_FORMAT_VERSION_MINOR 0



/* Error codes */
#define TSK_ERR_GENERIC                                             -1
#define TSK_ERR_NO_MEMORY                                           -2
#define TSK_ERR_IO                                                  -3
#define TSK_ERR_FILE_FORMAT                                         -4
#define TSK_ERR_REFCOUNT_NONZERO                                    -5
#define TSK_ERR_BAD_STATE                                           -6
#define TSK_ERR_BAD_PARAM_VALUE                                     -7
#define TSK_ERR_OUT_OF_BOUNDS                                       -8
#define TSK_ERR_BUFFER_OVERFLOW                                     -9
#define TSK_ERR_UNSORTED_DEMOGRAPHIC_EVENTS                         -10
#define TSK_ERR_POPULATION_OVERFLOW                                 -11
#define TSK_ERR_LINKS_OVERFLOW                                      -12
#define TSK_ERR_HDF5                                                -13
#define TSK_ERR_POPULATION_OUT_OF_BOUNDS                            -14
#define TSK_ERR_DUPLICATE_SAMPLE                                    -15
#define TSK_ERR_BAD_ORDERING                                        -16
#define TSK_ERR_BAD_MUTATION                                        -17
#define TSK_ERR_UNSUPPORTED_OPERATION                               -18
#define TSK_ERR_BAD_POPULATION_CONFIGURATION                        -19
#define TSK_ERR_BAD_MIGRATION_MATRIX                                -20
#define TSK_ERR_BAD_MIGRATION_MATRIX_INDEX                          -21
#define TSK_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX                     -22
#define TSK_ERR_INFINITE_WAITING_TIME                               -23
#define TSK_ERR_ASSERTION_FAILED                                    -24
#define TSK_ERR_SOURCE_DEST_EQUAL                                   -25
#define TSK_ERR_BAD_RECOMBINATION_MAP                               -26
#define TSK_ERR_BAD_POPULATION_SIZE                                 -27
#define TSK_ERR_BAD_SAMPLES                                         -28
#define TSK_ERR_BAD_TABLE_POSITION                                  -29
#define TSK_ERR_FILE_VERSION_TOO_OLD                                -30
#define TSK_ERR_FILE_VERSION_TOO_NEW                                -31
#define TSK_ERR_CANNOT_SIMPLIFY                                     -32
#define TSK_ERR_BAD_MODEL                                           -33
#define TSK_ERR_NULL_PARENT                                         -34
#define TSK_ERR_NULL_CHILD                                          -35
#define TSK_ERR_EDGES_NOT_SORTED_PARENT_TIME                        -36
#define TSK_ERR_EDGES_NONCONTIGUOUS_PARENTS                         -37
#define TSK_ERR_EDGES_NOT_SORTED_CHILD                              -38
#define TSK_ERR_EDGES_NOT_SORTED_LEFT                               -39
#define TSK_ERR_BAD_NODE_TIME_ORDERING                              -40
#define TSK_ERR_BAD_EDGE_INTERVAL                                   -41
#define TSK_ERR_DUPLICATE_EDGES                                     -42
#define TSK_ERR_NOT_INITIALISED                                     -43
#define TSK_ERR_BAD_OFFSET                                          -44
#define TSK_ERR_TOO_MANY_ALLELES                                    -45
#define TSK_ERR_DUPLICATE_MUTATION_NODES                            -46
#define TSK_ERR_NONBINARY_MUTATIONS_UNSUPPORTED                     -47
#define TSK_ERR_INCONSISTENT_MUTATIONS                              -48
#define TSK_ERR_INSUFFICIENT_SAMPLES                                -49
#define TSK_ERR_ZERO_RECORDS                                        -50
#define TSK_ERR_COORDINATE_NOT_FOUND                                -51
#define TSK_ERR_BAD_NODES_ARRAY                                     -52
#define TSK_ERR_BAD_CHILDREN_ARRAY                                  -53
#define TSK_ERR_SITE_OUT_OF_BOUNDS                                  -54
#define TSK_ERR_NODE_OUT_OF_BOUNDS                                  -55
#define TSK_ERR_LENGTH_MISMATCH                                     -56
#define TSK_ERR_DUPLICATE_SITE_POSITION                             -57
#define TSK_ERR_NON_SINGLE_CHAR_MUTATION                            -58
#define TSK_ERR_UNSORTED_SITES                                      -59
#define TSK_ERR_BAD_SITE_POSITION                                   -60
#define TSK_ERR_UNSORTED_MUTATIONS                                  -61
#define TSK_ERR_UNDEFINED_MULTIPLE_MERGER_COALESCENT                -62
#define TSK_ERR_EDGESETS_FOR_PARENT_NOT_ADJACENT                    -63
#define TSK_ERR_BAD_EDGESET_CONTRADICTORY_CHILDREN                  -64
#define TSK_ERR_BAD_EDGESET_OVERLAPPING_PARENT                      -65
#define TSK_ERR_BAD_SEQUENCE_LENGTH                                 -66
#define TSK_ERR_RIGHT_GREATER_SEQ_LENGTH                            -67
#define TSK_ERR_MUTATION_OUT_OF_BOUNDS                              -68
#define TSK_ERR_MUTATION_PARENT_DIFFERENT_SITE                      -69
#define TSK_ERR_MUTATION_PARENT_EQUAL                               -70
#define TSK_ERR_MUTATION_PARENT_AFTER_CHILD                         -71
#define TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS                            -72
#define TSK_ERR_GENERATE_UUID                                       -73
#define TSK_ERR_BAD_EDGE_INDEX                                      -74
#define TSK_ERR_LEFT_LESS_ZERO                                      -75
#define TSK_ERR_TABLES_NOT_INDEXED                                  -76
#define TSK_ERR_SIMPLIFY_MIGRATIONS_NOT_SUPPORTED                   -77
#define TSK_ERR_INCOMPATIBLE_FROM_TS                                -78
#define TSK_ERR_BAD_START_TIME_FROM_TS                              -79
#define TSK_ERR_BAD_START_TIME                                      -80
#define TSK_ERR_BAD_DEMOGRAPHIC_EVENT_TIME                          -81
#define TSK_ERR_RECOMB_MAP_TOO_COARSE                               -82
#define TSK_ERR_TIME_TRAVEL                                         -83
#define TSK_ERR_ONLY_INFINITE_SITES                                 -84


#define tsk_safe_free(pointer) __tsk_safe_free((void **) &(pointer))
#define TSK_MAX(a,b) ((a) > (b) ? (a) : (b))
#define TSK_MIN(a,b) ((a) < (b) ? (a) : (b))


/* This bit is 0 for any errors originating from kastore */
#define TSK_KAS_ERR_BIT 14

int tsk_set_kas_error(int err);
bool tsk_is_kas_error(int err);
const char * tsk_strerror(int err);
void __tsk_safe_free(void **ptr);


/* This is a simple allocator that is optimised to efficiently allocate a
 * large number of small objects without large numbers of calls to malloc.
 * The allocator mallocs memory in chunks of a configurable size. When
 * responding to calls to get(), it will return a chunk of this memory.
 * This memory cannot be subsequently handed back to the allocator. However,
 * all memory allocated by the allocator can be returned at once by calling
 * reset.
 */

typedef struct {
    size_t chunk_size;        /* number of bytes per chunk */
    size_t top;               /* the offset of the next available byte in the current chunk */
    size_t current_chunk;     /* the index of the chunk currently being used */
    size_t total_size;        /* the total number of bytes allocated + overhead. */
    size_t total_allocated;   /* the total number of bytes allocated. */
    size_t num_chunks;        /* the number of memory chunks. */
    char **mem_chunks;        /* the memory chunks */
} tsk_blkalloc_t;

extern void tsk_blkalloc_print_state(tsk_blkalloc_t *self, FILE *out);
extern int tsk_blkalloc_reset(tsk_blkalloc_t *self);
extern int tsk_blkalloc_alloc(tsk_blkalloc_t *self, size_t chunk_size);
extern void * tsk_blkalloc_get(tsk_blkalloc_t *self, size_t size);
extern void tsk_blkalloc_free(tsk_blkalloc_t *self);

size_t tsk_search_sorted(const double *array, size_t size, double value);

#define TSK_UUID_SIZE 36
int tsk_generate_uuid(char *dest, int flags);

#endif
