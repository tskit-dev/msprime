#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <kastore.h>
#include "tsk_core.h"

#define UUID_NUM_BYTES 16

#if defined(_WIN32)

#include <windows.h>
#include <wincrypt.h>

static int TSK_WARN_UNUSED
get_random_bytes(uint8_t *buf)
{
    /* Based on CPython's code in bootstrap_hash.c */
    int ret = TSK_ERR_GENERATE_UUID;
    HCRYPTPROV hCryptProv = NULL;

    if (!CryptAcquireContext(&hCryptProv, NULL, NULL, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT)) {
        goto out;
    }
    if (!CryptGenRandom(hCryptProv, (DWORD) UUID_NUM_BYTES, buf)) {
        goto out;
    }
    if (!CryptReleaseContext(hCryptProv, 0)) {
        hCryptProv = NULL;
        goto out;
    }
    hCryptProv = NULL;
    ret = 0;
out:
    if (hCryptProv != NULL) {
        CryptReleaseContext(hCryptProv, 0);
    }
    return ret;
}

#else

/* Assuming the existance of /dev/urandom on Unix platforms */
static int TSK_WARN_UNUSED
get_random_bytes(uint8_t *buf)
{
    int ret = TSK_ERR_GENERATE_UUID;
    FILE *f = fopen("/dev/urandom", "r");

    if (f == NULL) {
        goto out;
    }
    if (fread(buf, UUID_NUM_BYTES, 1, f) != 1) {
        goto out;
    }
    if (fclose(f) != 0) {
        goto out;
    }
    ret = 0;
out:
    return ret;
}

#endif

/* Generate a new UUID4 using a system-generated source of randomness.
 * Note that this function writes a NULL terminator to the end of this
 * string, so that the total length of the buffer must be 37 bytes.
 */
int
tsk_generate_uuid(char *dest, int TSK_UNUSED(flags))
{
    int ret = 0;
    uint8_t buf[UUID_NUM_BYTES];
    const char *pattern =
        "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x";

    ret = get_random_bytes(buf);
    if (ret != 0) {
        goto out;
    }
    if (snprintf(dest, TSK_UUID_SIZE + 1, pattern,
            buf[0], buf[1], buf[2], buf[3],
            buf[4], buf[5], buf[6], buf[7],
            buf[8], buf[9], buf[10], buf[11],
            buf[12], buf[13], buf[14], buf[15]) < 0) {
        ret = TSK_ERR_GENERATE_UUID;
        goto out;
    }
out:
    return ret;
}
static const char *
tsk_strerror_internal(int err)
{
    const char *ret = "Unknown error";

    switch (err) {
        case 0:
            ret = "Normal exit condition. This is not an error!";
            goto out;
        case TSK_ERR_NO_MEMORY:
            ret = "Out of memory.";
            break;
        case TSK_ERR_GENERIC:
            ret = "Generic error; please file a bug report";
            break;
        case TSK_ERR_FILE_FORMAT:
            ret = "File format error";
            break;
        case TSK_ERR_BAD_STATE:
            ret = "Bad simulator state. Initialise or reset must be called.";
            break;
        case TSK_ERR_BUFFER_OVERFLOW:
            ret = "Supplied buffer if too small";
            break;
        case TSK_ERR_UNSORTED_DEMOGRAPHIC_EVENTS:
            ret = "Demographic events must be time sorted.";
            break;
        case TSK_ERR_POPULATION_OVERFLOW:
            ret = "Population Overflow occurred.";
            break;
        case TSK_ERR_LINKS_OVERFLOW:
            ret = "Links Overflow occurred.";
            break;
        case TSK_ERR_OUT_OF_BOUNDS:
            ret = "Object reference out of bounds";
            break;
        case TSK_ERR_BAD_ORDERING:
            ret = "Bad record ordering requested";
            break;
        case TSK_ERR_BAD_MUTATION:
            ret = "Bad mutation provided";
            break;
        case TSK_ERR_BAD_PARAM_VALUE:
            ret = "Bad parameter value provided";
            break;
        case TSK_ERR_UNSUPPORTED_OPERATION:
            ret = "Operation cannot be performed in current configuration";
            break;
        case TSK_ERR_BAD_POPULATION_CONFIGURATION:
            ret = "Bad population configuration provided.";
            break;
        case TSK_ERR_BAD_POPULATION_SIZE:
            ret = "Bad population size provided. Must be > 0.";
            break;
        case TSK_ERR_POPULATION_OUT_OF_BOUNDS:
            ret = "Population ID out of bounds.";
            break;
        case TSK_ERR_BAD_MIGRATION_MATRIX:
            ret = "Bad migration matrix provided.";
            break;
        case TSK_ERR_BAD_MIGRATION_MATRIX_INDEX:
            ret = "Bad migration matrix index provided.";
            break;
        case TSK_ERR_DIAGONAL_MIGRATION_MATRIX_INDEX:
            ret = "Cannot set diagonal migration matrix elements.";
            break;
        case TSK_ERR_INFINITE_WAITING_TIME:
            ret = "Infinite waiting time until next simulation event.";
            break;
        case TSK_ERR_ASSERTION_FAILED:
            ret = "Internal error; please file a bug report.";
            break;
        case TSK_ERR_SOURCE_DEST_EQUAL:
            ret = "Source and destination populations equal.";
            break;
        case TSK_ERR_BAD_RECOMBINATION_MAP:
            ret = "Bad recombination map provided.";
            break;
        case TSK_ERR_INSUFFICIENT_SAMPLES:
            ret = "At least two samples needed.";
            break;
        case TSK_ERR_ZERO_RECORDS:
            ret = "At least one record must be supplied";
            break;
        case TSK_ERR_EDGES_NOT_SORTED_PARENT_TIME:
            ret = "Edges must be listed in (time[parent], child, left) order;"
                " time[parent] order violated";
            break;
        case TSK_ERR_EDGES_NONCONTIGUOUS_PARENTS:
            ret = "All edges for a given parent must be contiguous";
            break;
        case TSK_ERR_EDGES_NOT_SORTED_CHILD:
            ret = "Edges must be listed in (time[parent], child, left) order;"
                " child order violated";
            break;
        case TSK_ERR_EDGES_NOT_SORTED_LEFT:
            ret = "Edges must be listed in (time[parent], child, left) order;"
                " left order violated";
            break;
        case TSK_ERR_NULL_PARENT:
            ret = "Edge in parent is null.";
            break;
        case TSK_ERR_NULL_CHILD:
            ret = "Edge in parent is null.";
            break;
        case TSK_ERR_BAD_NODE_TIME_ORDERING:
            ret = "time[parent] must be greater than time[child]";
            break;
        case TSK_ERR_BAD_EDGE_INTERVAL:
            ret = "Bad edge interval where right <= left";
            break;
        case TSK_ERR_DUPLICATE_EDGES:
            ret = "Duplicate edges provided.";
            break;
        case TSK_ERR_CANNOT_SIMPLIFY:
            ret = "Cannot simplify the tree sequence; no output records.";
            break;
        case TSK_ERR_BAD_SAMPLES:
            ret = "Bad sample configuration provided.";
            break;
        case TSK_ERR_FILE_VERSION_TOO_OLD:
            ret = "tskit file version too old. Please upgrade using the "
                "'msp upgrade' command";
            break;
        case TSK_ERR_FILE_VERSION_TOO_NEW:
            ret = "tskit file version is too new for this instance. "
                "Please upgrade msprime to the latest version.";
            break;
        case TSK_ERR_DUPLICATE_SAMPLE:
            ret = "Duplicate value provided in tracked leaf list.";
            break;
        case TSK_ERR_REFCOUNT_NONZERO:
            ret = "Cannot change the state of the tree sequence when "
                "other objects reference it. Make sure all trees are freed first.";
            break;
        case TSK_ERR_BAD_MODEL:
            ret = "Model error. Either a bad model, or the requested operation "
                "is not supported for the current model";
            break;
        case TSK_ERR_NOT_INITIALISED:
            ret = "object not initialised. Please file a bug report.";
            break;
        case TSK_ERR_DUPLICATE_MUTATION_NODES:
            ret = "Cannot have more than one mutation at a node for a given site.";
            break;
        case TSK_ERR_NONBINARY_MUTATIONS_UNSUPPORTED:
            ret = "Only binary mutations are supported for this operation.";
            break;
        case TSK_ERR_INCONSISTENT_MUTATIONS:
            ret = "Inconsistent mutations: state already equal to derived state.";
            break;
        case TSK_ERR_COORDINATE_NOT_FOUND:
            ret = "Coordinate not found.";
            break;
        case TSK_ERR_BAD_NODES_ARRAY:
            ret = "Malformed nodes array.";
            break;
        case TSK_ERR_BAD_CHILDREN_ARRAY:
            ret = "Malformed array of children.";
            break;
        case TSK_ERR_SITE_OUT_OF_BOUNDS:
            ret = "Site out of bounds";
            break;
        case TSK_ERR_NODE_OUT_OF_BOUNDS:
            ret = "Node out of bounds";
            break;
        case TSK_ERR_LENGTH_MISMATCH:
            ret = "Mismatch in stored total column length and sum of row lengths";
            break;
        case TSK_ERR_NON_SINGLE_CHAR_MUTATION:
            ret = "Only single char mutations supported.";
            break;
        case TSK_ERR_UNSORTED_SITES:
            ret = "Sites must be provided in strictly increasing position order.";
            break;
        case TSK_ERR_BAD_SITE_POSITION:
            ret = "Sites positions must be between 0 and sequence_length";
            break;
        case TSK_ERR_UNSORTED_MUTATIONS:
            ret = "Mutations must be provided in non-decreasing site order";
            break;
        case TSK_ERR_EDGESETS_FOR_PARENT_NOT_ADJACENT:
            ret = "All edges for a given parent must be adjacent.";
            break;
        case TSK_ERR_BAD_EDGESET_CONTRADICTORY_CHILDREN:
            ret = "Bad edges: contradictory children for a given parent over "
                "an interval.";
            break;
        case TSK_ERR_BAD_EDGESET_OVERLAPPING_PARENT:
            ret = "Bad edges: multiple definitions of a given parent over an interval";
            break;
        case TSK_ERR_BAD_SEQUENCE_LENGTH:
            ret = "Sequence length must be > 0.";
            break;
        case TSK_ERR_LEFT_LESS_ZERO:
            ret = "Left coordinate must be >= 0";
            break;
        case TSK_ERR_RIGHT_GREATER_SEQ_LENGTH:
            ret = "Right coordinate > sequence length.";
            break;
        case TSK_ERR_MUTATION_OUT_OF_BOUNDS:
            ret = "mutation ID out of bounds";
            break;
        case TSK_ERR_MUTATION_PARENT_DIFFERENT_SITE:
            ret = "Specified parent mutation is at a different site.";
            break;
        case TSK_ERR_MUTATION_PARENT_EQUAL:
            ret = "Parent mutation refers to itself.";
            break;
        case TSK_ERR_MUTATION_PARENT_AFTER_CHILD:
            ret = "Parent mutation ID must be < current ID.";
            break;
        case TSK_ERR_BAD_OFFSET:
            ret = "Bad offset provided in input array.";
            break;
        case TSK_ERR_TOO_MANY_ALLELES:
            ret = "Cannot have more than 255 alleles.";
            break;
        case TSK_ERR_INDIVIDUAL_OUT_OF_BOUNDS:
            ret = "Individual ID out of bounds";
            break;
        case TSK_ERR_GENERATE_UUID:
            ret = "Error generating UUID";
            break;
        case TSK_ERR_DUPLICATE_SITE_POSITION:
            ret = "Duplicate site positions.";
            break;
        case TSK_ERR_BAD_TABLE_POSITION:
            ret = "Bad table position provided to truncate/reset.";
            break;
        case TSK_ERR_BAD_EDGE_INDEX:
            ret = "Invalid edge index value.";
            break;
        case TSK_ERR_TABLES_NOT_INDEXED:
            ret = "Table collection must be indexed.";
            break;
        case TSK_ERR_SIMPLIFY_MIGRATIONS_NOT_SUPPORTED:
            ret = "Migrations currently not supported in simplify. Please open an "
                "issue on GitHub if this operation is important to you.";
            break;
        case TSK_ERR_INCOMPATIBLE_FROM_TS:
            ret = "The specified tree sequence is not a compatible starting point "
                "for the current simulation";
            break;
        case TSK_ERR_BAD_START_TIME_FROM_TS:
            ret = "The specified start_time and from_ts are not compatible. All "
                "node times in the tree sequence must be <= start_time.";
            break;
        case TSK_ERR_BAD_START_TIME:
            ret = "start_time must be >= 0.";
            break;
        case TSK_ERR_BAD_DEMOGRAPHIC_EVENT_TIME:
            ret = "demographic event time must be >= start_time.";
            break;
        case TSK_ERR_RECOMB_MAP_TOO_COARSE:
            ret = "The specified recombination map is cannot translate the coordinates"
                "for the specified tree sequence. It is either too coarse (num_loci "
                "is too small) or contains zero recombination rates. Please either "
                "increase the number of loci or recombination rate";
            break;
        case TSK_ERR_TIME_TRAVEL:
            ret = "The simulation model supplied resulted in a parent node having "
                "a time value <= to its child. This can occur either as a result "
                "of multiple bottlenecks happening at the same time or because of "
                "numerical imprecision with very small population sizes.";
            break;
        case TSK_ERR_ONLY_INFINITE_SITES:
            ret = "Only infinite sites mutations are supported for this operation.";
            break;

        case TSK_ERR_IO:
            if (errno != 0) {
                ret = strerror(errno);
            } else {
                ret = "Unspecified IO error";
            }
            break;
        default:
            ret = "Error occurred generating error string. Please file a bug "
                "report!";
            break;
    }
out:
    return ret;
}

int
tsk_set_kas_error(int err)
{
    /* Flip this bit. As the error is negative, this sets the bit to 0 */
    return err ^ (1 << TSK_KAS_ERR_BIT);
}

bool
tsk_is_kas_error(int err)
{
    return !(err & (1 << TSK_KAS_ERR_BIT));
}

const char *
tsk_strerror(int err)
{
    if (tsk_is_kas_error(err)) {
        err ^= (1 << TSK_KAS_ERR_BIT);
        return kas_strerror(err);
    } else {
        return tsk_strerror_internal(err);
    }
}

void
__tsk_safe_free(void **ptr) {
    if (ptr != NULL) {
        if (*ptr != NULL) {
            free(*ptr);
            *ptr = NULL;
        }
    }
}


/* Block allocator. Simple allocator when we lots of chunks of memory
 * and don't need to free them individually.
 */

void
tsk_blkalloc_print_state(tsk_blkalloc_t *self, FILE *out)
{
    fprintf(out, "Block allocator%p::\n", (void *) self);
    fprintf(out, "\ttop = %d\n", (int) self->top);
    fprintf(out, "\tchunk_size = %d\n", (int) self->chunk_size);
    fprintf(out, "\tnum_chunks = %d\n", (int) self->num_chunks);
    fprintf(out, "\ttotal_allocated = %d\n", (int) self->total_allocated);
    fprintf(out, "\ttotal_size = %d\n", (int) self->total_size);
}

int TSK_WARN_UNUSED
tsk_blkalloc_reset(tsk_blkalloc_t *self)
{
    int ret = 0;

    self->top = 0;
    self->current_chunk = 0;
    self->total_allocated = 0;
    return ret;
}

int TSK_WARN_UNUSED
tsk_blkalloc_alloc(tsk_blkalloc_t *self, size_t chunk_size)
{
    int ret = 0;

    assert(chunk_size > 0);
    memset(self, 0, sizeof(tsk_blkalloc_t));
    self->chunk_size = chunk_size;
    self->top = 0;
    self->current_chunk = 0;
    self->total_allocated = 0;
    self->total_size = 0;
    self->num_chunks = 0;
    self->mem_chunks = malloc(sizeof(char *));
    if (self->mem_chunks == NULL) {
        ret = TSK_ERR_NO_MEMORY;
        goto out;
    }
    self->mem_chunks[0] = malloc(chunk_size);
    if (self->mem_chunks[0] == NULL) {
        ret = TSK_ERR_NO_MEMORY;
        goto out;
    }
    self->num_chunks = 1;
    self->total_size = chunk_size + sizeof(void *);
out:
    return ret;
}

void * TSK_WARN_UNUSED
tsk_blkalloc_get(tsk_blkalloc_t *self, size_t size)
{
    void *ret = NULL;
    void *p;

    assert(size < self->chunk_size);
    if ((self->top + size) > self->chunk_size) {
        if (self->current_chunk == (self->num_chunks - 1)) {
            p = realloc(self->mem_chunks, (self->num_chunks + 1) * sizeof(void *));
            if (p == NULL) {
                goto out;
            }
            self->mem_chunks = p;
            p = malloc(self->chunk_size);
            if (p == NULL) {
                goto out;
            }
            self->mem_chunks[self->num_chunks] = p;
            self->num_chunks++;
            self->total_size += self->chunk_size + sizeof(void *);
        }
        self->current_chunk++;
        self->top = 0;
    }
    ret = self->mem_chunks[self->current_chunk] + self->top;
    self->top += size;
    self->total_allocated += size;
out:
    return ret;
}

void
tsk_blkalloc_free(tsk_blkalloc_t *self)
{
    size_t j;

    for (j = 0; j < self->num_chunks; j++) {
        if (self->mem_chunks[j] != NULL) {
            free(self->mem_chunks[j]);
        }
    }
    if (self->mem_chunks != NULL) {
        free(self->mem_chunks);
    }
}

/* Mirrors the semantics of numpy's searchsorted function. Uses binary
 * search to find the index of the closest value in the array. */
size_t
tsk_search_sorted(const double *restrict array, size_t size, double value)
{
    int64_t upper = (int64_t) size;
    int64_t lower = 0;
    int64_t offset = 0;
    int64_t mid;

    if (upper == 0) {
        return 0;
    }

    while (upper - lower > 1) {
        mid = (upper + lower) / 2;
        if (value >= array[mid]) {
            lower = mid;
        } else {
            upper = mid;
        }
    }
    offset = (int64_t) (array[lower] < value);
    return (size_t) (lower + offset);
}
