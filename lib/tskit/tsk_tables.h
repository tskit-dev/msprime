#ifndef TSK_TABLES_H
#define TSK_TABLES_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

#include <kastore.h>

#include "tsk_core.h"

typedef int32_t tsk_id_t;
typedef uint32_t tsk_tbl_size_t;

/****************************************************************************/
/* Definitions for the basic objects */
/****************************************************************************/

typedef struct {
    tsk_id_t id;
    uint32_t flags;
    double *location;
    tsk_tbl_size_t location_length;
    const char *metadata;
    tsk_tbl_size_t metadata_length;
    tsk_id_t *nodes;
    tsk_tbl_size_t nodes_length;
} tsk_individual_t;

typedef struct {
    tsk_id_t id;
    uint32_t flags;
    double time;
    tsk_id_t population;
    tsk_id_t individual;
    const char *metadata;
    tsk_tbl_size_t metadata_length;
} tsk_node_t;

typedef struct {
    tsk_id_t id;
    tsk_id_t parent;
    tsk_id_t child;
    double left;
    double right;
} tsk_edge_t;

typedef struct {
    tsk_id_t id;
    tsk_id_t site;
    tsk_id_t node;
    tsk_id_t parent;
    const char *derived_state;
    tsk_tbl_size_t derived_state_length;
    const char *metadata;
    tsk_tbl_size_t metadata_length;
} tsk_mutation_t;

typedef struct {
    tsk_id_t id;
    double position;
    const char *ancestral_state;
    tsk_tbl_size_t ancestral_state_length;
    const char *metadata;
    tsk_tbl_size_t metadata_length;
    tsk_mutation_t *mutations;
    tsk_tbl_size_t mutations_length;
} tsk_site_t;

typedef struct {
    tsk_id_t id;
    tsk_id_t source;
    tsk_id_t dest;
    tsk_id_t node;
    double left;
    double right;
    double time;
} tsk_migration_t;

typedef struct {
    tsk_id_t id;
    const char *metadata;
    tsk_tbl_size_t metadata_length;
} tsk_population_t;

typedef struct {
    tsk_id_t id;
    const char *timestamp;
    tsk_tbl_size_t timestamp_length;
    const char *record;
    tsk_tbl_size_t record_length;
} tsk_provenance_t;

/****************************************************************************/
/* Table definitions */
/****************************************************************************/

typedef struct {
    tsk_tbl_size_t num_rows;
    tsk_tbl_size_t max_rows;
    tsk_tbl_size_t max_rows_increment;
    tsk_tbl_size_t location_length;
    tsk_tbl_size_t max_location_length;
    tsk_tbl_size_t max_location_length_increment;
    tsk_tbl_size_t metadata_length;
    tsk_tbl_size_t max_metadata_length;
    tsk_tbl_size_t max_metadata_length_increment;
    uint32_t *flags;
    double *location;
    tsk_tbl_size_t *location_offset;
    char *metadata;
    tsk_tbl_size_t *metadata_offset;
} tsk_individual_tbl_t;

typedef struct {
    tsk_tbl_size_t num_rows;
    tsk_tbl_size_t max_rows;
    tsk_tbl_size_t max_rows_increment;
    tsk_tbl_size_t metadata_length;
    tsk_tbl_size_t max_metadata_length;
    tsk_tbl_size_t max_metadata_length_increment;
    uint32_t *flags;
    double *time;
    tsk_id_t *population;
    tsk_id_t *individual;
    char *metadata;
    tsk_tbl_size_t *metadata_offset;
} tsk_node_tbl_t;

typedef struct {
    tsk_tbl_size_t num_rows;
    tsk_tbl_size_t max_rows;
    tsk_tbl_size_t max_rows_increment;
    tsk_tbl_size_t ancestral_state_length;
    tsk_tbl_size_t max_ancestral_state_length;
    tsk_tbl_size_t max_ancestral_state_length_increment;
    tsk_tbl_size_t metadata_length;
    tsk_tbl_size_t max_metadata_length;
    tsk_tbl_size_t max_metadata_length_increment;
    double *position;
    char *ancestral_state;
    tsk_tbl_size_t *ancestral_state_offset;
    char *metadata;
    tsk_tbl_size_t *metadata_offset;
} tsk_site_tbl_t;

typedef struct {
    tsk_tbl_size_t num_rows;
    tsk_tbl_size_t max_rows;
    tsk_tbl_size_t max_rows_increment;
    tsk_tbl_size_t derived_state_length;
    tsk_tbl_size_t max_derived_state_length;
    tsk_tbl_size_t max_derived_state_length_increment;
    tsk_tbl_size_t metadata_length;
    tsk_tbl_size_t max_metadata_length;
    tsk_tbl_size_t max_metadata_length_increment;
    tsk_id_t *node;
    tsk_id_t *site;
    tsk_id_t *parent;
    char *derived_state;
    tsk_tbl_size_t *derived_state_offset;
    char *metadata;
    tsk_tbl_size_t *metadata_offset;
} tsk_mutation_tbl_t;

typedef struct {
    tsk_tbl_size_t num_rows;
    tsk_tbl_size_t max_rows;
    tsk_tbl_size_t max_rows_increment;
    double *left;
    double *right;
    tsk_id_t *parent;
    tsk_id_t *child;
} tsk_edge_tbl_t;

typedef struct {
    tsk_tbl_size_t num_rows;
    tsk_tbl_size_t max_rows;
    tsk_tbl_size_t max_rows_increment;
    tsk_id_t *source;
    tsk_id_t *dest;
    tsk_id_t *node;
    double *left;
    double *right;
    double *time;
} tsk_migration_tbl_t;

typedef struct {
    tsk_tbl_size_t num_rows;
    tsk_tbl_size_t max_rows;
    tsk_tbl_size_t max_rows_increment;
    tsk_tbl_size_t metadata_length;
    tsk_tbl_size_t max_metadata_length;
    tsk_tbl_size_t max_metadata_length_increment;
    char *metadata;
    tsk_tbl_size_t *metadata_offset;
} tsk_population_tbl_t;

typedef struct {
    tsk_tbl_size_t num_rows;
    tsk_tbl_size_t max_rows;
    tsk_tbl_size_t max_rows_increment;
    tsk_tbl_size_t timestamp_length;
    tsk_tbl_size_t max_timestamp_length;
    tsk_tbl_size_t max_timestamp_length_increment;
    tsk_tbl_size_t record_length;
    tsk_tbl_size_t max_record_length;
    tsk_tbl_size_t max_record_length_increment;
    char *timestamp;
    tsk_tbl_size_t *timestamp_offset;
    char *record;
    tsk_tbl_size_t *record_offset;
} tsk_provenance_tbl_t;

typedef struct {
    double sequence_length;
    char *file_uuid;
    tsk_individual_tbl_t *individuals;
    tsk_node_tbl_t *nodes;
    tsk_edge_tbl_t *edges;
    tsk_migration_tbl_t *migrations;
    tsk_site_tbl_t *sites;
    tsk_mutation_tbl_t *mutations;
    tsk_population_tbl_t *populations;
    tsk_provenance_tbl_t *provenances;
    struct {
        tsk_id_t *edge_insertion_order;
        tsk_id_t *edge_removal_order;
        bool malloced_locally;
    } indexes;
    kastore_t *store;
    /* TODO Add in reserved space for future tables. */
} tsk_tbl_collection_t;

typedef struct {
    tsk_tbl_size_t individuals;
    tsk_tbl_size_t nodes;
    tsk_tbl_size_t edges;
    tsk_tbl_size_t migrations;
    tsk_tbl_size_t sites;
    tsk_tbl_size_t mutations;
    tsk_tbl_size_t populations;
    tsk_tbl_size_t provenances;
    /* TODO add reserved space for future tables. */
} tsk_tbl_collection_position_t;


/****************************************************************************/
/* Function signatures */
/****************************************************************************/

int tsk_individual_tbl_alloc(tsk_individual_tbl_t *self, int flags);
int tsk_individual_tbl_set_max_rows_increment(tsk_individual_tbl_t *self, size_t max_rows_increment);
int tsk_individual_tbl_set_max_metadata_length_increment(tsk_individual_tbl_t *self,
        size_t max_metadata_length_increment);
int tsk_individual_tbl_set_max_location_length_increment(tsk_individual_tbl_t *self,
        size_t max_location_length_increment);
tsk_id_t tsk_individual_tbl_add_row(tsk_individual_tbl_t *self, uint32_t flags,
        double *location, size_t location_length,
        const char *metadata, size_t metadata_length);
int tsk_individual_tbl_set_columns(tsk_individual_tbl_t *self, size_t num_rows, uint32_t *flags,
        double *location, tsk_tbl_size_t *location_length,
        const char *metadata, tsk_tbl_size_t *metadata_length);
int tsk_individual_tbl_append_columns(tsk_individual_tbl_t *self, size_t num_rows, uint32_t *flags,
        double *location, tsk_tbl_size_t *location_length,
        const char *metadata, tsk_tbl_size_t *metadata_length);
int tsk_individual_tbl_clear(tsk_individual_tbl_t *self);
int tsk_individual_tbl_truncate(tsk_individual_tbl_t *self, size_t num_rows);
int tsk_individual_tbl_free(tsk_individual_tbl_t *self);
int tsk_individual_tbl_dump_text(tsk_individual_tbl_t *self, FILE *out);
int tsk_individual_tbl_copy(tsk_individual_tbl_t *self, tsk_individual_tbl_t *dest);
void tsk_individual_tbl_print_state(tsk_individual_tbl_t *self, FILE *out);
bool tsk_individual_tbl_equals(tsk_individual_tbl_t *self, tsk_individual_tbl_t *other);
int tsk_individual_tbl_get_row(tsk_individual_tbl_t *self, size_t index,
        tsk_individual_t *row);

int tsk_node_tbl_alloc(tsk_node_tbl_t *self, int flags);
int tsk_node_tbl_set_max_rows_increment(tsk_node_tbl_t *self, size_t max_rows_increment);
int tsk_node_tbl_set_max_metadata_length_increment(tsk_node_tbl_t *self,
        size_t max_metadata_length_increment);
tsk_id_t tsk_node_tbl_add_row(tsk_node_tbl_t *self, uint32_t flags, double time,
        tsk_id_t population, tsk_id_t individual,
        const char *metadata, size_t metadata_length);
int tsk_node_tbl_set_columns(tsk_node_tbl_t *self, size_t num_rows,
        uint32_t *flags, double *time,
        tsk_id_t *population, tsk_id_t *individual,
        const char *metadata, tsk_tbl_size_t *metadata_length);
int tsk_node_tbl_append_columns(tsk_node_tbl_t *self, size_t num_rows,
        uint32_t *flags, double *time,
        tsk_id_t *population, tsk_id_t *individual,
        const char *metadata, tsk_tbl_size_t *metadata_length);
int tsk_node_tbl_clear(tsk_node_tbl_t *self);
int tsk_node_tbl_truncate(tsk_node_tbl_t *self, size_t num_rows);
int tsk_node_tbl_free(tsk_node_tbl_t *self);
int tsk_node_tbl_dump_text(tsk_node_tbl_t *self, FILE *out);
int tsk_node_tbl_copy(tsk_node_tbl_t *self, tsk_node_tbl_t *dest);
void tsk_node_tbl_print_state(tsk_node_tbl_t *self, FILE *out);
bool tsk_node_tbl_equals(tsk_node_tbl_t *self, tsk_node_tbl_t *other);
int tsk_node_tbl_get_row(tsk_node_tbl_t *self, size_t index, tsk_node_t *row);

int tsk_edge_tbl_alloc(tsk_edge_tbl_t *self, int flags);
int tsk_edge_tbl_set_max_rows_increment(tsk_edge_tbl_t *self, size_t max_rows_increment);
tsk_id_t tsk_edge_tbl_add_row(tsk_edge_tbl_t *self, double left, double right, tsk_id_t parent,
        tsk_id_t child);
int tsk_edge_tbl_set_columns(tsk_edge_tbl_t *self, size_t num_rows, double *left,
        double *right, tsk_id_t *parent, tsk_id_t *child);
int tsk_edge_tbl_append_columns(tsk_edge_tbl_t *self, size_t num_rows, double *left,
        double *right, tsk_id_t *parent, tsk_id_t *child);
int tsk_edge_tbl_clear(tsk_edge_tbl_t *self);
int tsk_edge_tbl_truncate(tsk_edge_tbl_t *self, size_t num_rows);
int tsk_edge_tbl_free(tsk_edge_tbl_t *self);
int tsk_edge_tbl_dump_text(tsk_edge_tbl_t *self, FILE *out);
int tsk_edge_tbl_copy(tsk_edge_tbl_t *self, tsk_edge_tbl_t *dest);
void tsk_edge_tbl_print_state(tsk_edge_tbl_t *self, FILE *out);
bool tsk_edge_tbl_equals(tsk_edge_tbl_t *self, tsk_edge_tbl_t *other);
int tsk_edge_tbl_get_row(tsk_edge_tbl_t *self, size_t index, tsk_edge_t *row);

int tsk_site_tbl_alloc(tsk_site_tbl_t *self, int flags);
int tsk_site_tbl_set_max_rows_increment(tsk_site_tbl_t *self, size_t max_rows_increment);
int tsk_site_tbl_set_max_metadata_length_increment(tsk_site_tbl_t *self,
        size_t max_metadata_length_increment);
int tsk_site_tbl_set_max_ancestral_state_length_increment(tsk_site_tbl_t *self,
        size_t max_ancestral_state_length_increment);
tsk_id_t tsk_site_tbl_add_row(tsk_site_tbl_t *self, double position,
        const char *ancestral_state, tsk_tbl_size_t ancestral_state_length,
        const char *metadata, tsk_tbl_size_t metadata_length);
int tsk_site_tbl_set_columns(tsk_site_tbl_t *self, size_t num_rows, double *position,
        const char *ancestral_state, tsk_tbl_size_t *ancestral_state_length,
        const char *metadata, tsk_tbl_size_t *metadata_length);
int tsk_site_tbl_append_columns(tsk_site_tbl_t *self, size_t num_rows, double *position,
        const char *ancestral_state, tsk_tbl_size_t *ancestral_state_length,
        const char *metadata, tsk_tbl_size_t *metadata_length);
bool tsk_site_tbl_equals(tsk_site_tbl_t *self, tsk_site_tbl_t *other);
int tsk_site_tbl_clear(tsk_site_tbl_t *self);
int tsk_site_tbl_truncate(tsk_site_tbl_t *self, size_t num_rows);
int tsk_site_tbl_copy(tsk_site_tbl_t *self, tsk_site_tbl_t *dest);
int tsk_site_tbl_free(tsk_site_tbl_t *self);
int tsk_site_tbl_dump_text(tsk_site_tbl_t *self, FILE *out);
void tsk_site_tbl_print_state(tsk_site_tbl_t *self, FILE *out);
int tsk_site_tbl_get_row(tsk_site_tbl_t *self, size_t index, tsk_site_t *row);

void tsk_mutation_tbl_print_state(tsk_mutation_tbl_t *self, FILE *out);
int tsk_mutation_tbl_alloc(tsk_mutation_tbl_t *self, int flags);
int tsk_mutation_tbl_set_max_rows_increment(tsk_mutation_tbl_t *self, size_t max_rows_increment);
int tsk_mutation_tbl_set_max_metadata_length_increment(tsk_mutation_tbl_t *self,
        size_t max_metadata_length_increment);
int tsk_mutation_tbl_set_max_derived_state_length_increment(tsk_mutation_tbl_t *self,
        size_t max_derived_state_length_increment);
tsk_id_t tsk_mutation_tbl_add_row(tsk_mutation_tbl_t *self, tsk_id_t site,
        tsk_id_t node, tsk_id_t parent,
        const char *derived_state, tsk_tbl_size_t derived_state_length,
        const char *metadata, tsk_tbl_size_t metadata_length);
int tsk_mutation_tbl_set_columns(tsk_mutation_tbl_t *self, size_t num_rows,
        tsk_id_t *site, tsk_id_t *node, tsk_id_t *parent,
        const char *derived_state, tsk_tbl_size_t *derived_state_length,
        const char *metadata, tsk_tbl_size_t *metadata_length);
int tsk_mutation_tbl_append_columns(tsk_mutation_tbl_t *self, size_t num_rows,
        tsk_id_t *site, tsk_id_t *node, tsk_id_t *parent,
        const char *derived_state, tsk_tbl_size_t *derived_state_length,
        const char *metadata, tsk_tbl_size_t *metadata_length);
bool tsk_mutation_tbl_equals(tsk_mutation_tbl_t *self, tsk_mutation_tbl_t *other);
int tsk_mutation_tbl_clear(tsk_mutation_tbl_t *self);
int tsk_mutation_tbl_truncate(tsk_mutation_tbl_t *self, size_t num_rows);
int tsk_mutation_tbl_copy(tsk_mutation_tbl_t *self, tsk_mutation_tbl_t *dest);
int tsk_mutation_tbl_free(tsk_mutation_tbl_t *self);
int tsk_mutation_tbl_dump_text(tsk_mutation_tbl_t *self, FILE *out);
void tsk_mutation_tbl_print_state(tsk_mutation_tbl_t *self, FILE *out);
int tsk_mutation_tbl_get_row(tsk_mutation_tbl_t *self, size_t index, tsk_mutation_t *row);

int tsk_migration_tbl_alloc(tsk_migration_tbl_t *self, int flags);
int tsk_migration_tbl_set_max_rows_increment(tsk_migration_tbl_t *self, size_t max_rows_increment);
tsk_id_t tsk_migration_tbl_add_row(tsk_migration_tbl_t *self, double left,
        double right, tsk_id_t node, tsk_id_t source,
        tsk_id_t dest, double time);
int tsk_migration_tbl_set_columns(tsk_migration_tbl_t *self, size_t num_rows,
        double *left, double *right, tsk_id_t *node, tsk_id_t *source,
        tsk_id_t *dest, double *time);
int tsk_migration_tbl_append_columns(tsk_migration_tbl_t *self, size_t num_rows,
        double *left, double *right, tsk_id_t *node, tsk_id_t *source,
        tsk_id_t *dest, double *time);
int tsk_migration_tbl_clear(tsk_migration_tbl_t *self);
int tsk_migration_tbl_truncate(tsk_migration_tbl_t *self, size_t num_rows);
int tsk_migration_tbl_free(tsk_migration_tbl_t *self);
int tsk_migration_tbl_copy(tsk_migration_tbl_t *self, tsk_migration_tbl_t *dest);
int tsk_migration_tbl_dump_text(tsk_migration_tbl_t *self, FILE *out);
void tsk_migration_tbl_print_state(tsk_migration_tbl_t *self, FILE *out);
bool tsk_migration_tbl_equals(tsk_migration_tbl_t *self, tsk_migration_tbl_t *other);
int tsk_migration_tbl_get_row(tsk_migration_tbl_t *self, size_t index, tsk_migration_t *row);

int tsk_population_tbl_alloc(tsk_population_tbl_t *self, int flags);
int tsk_population_tbl_set_max_rows_increment(tsk_population_tbl_t *self, size_t max_rows_increment);
int tsk_population_tbl_set_max_metadata_length_increment(tsk_population_tbl_t *self,
        size_t max_metadata_length_increment);
tsk_id_t tsk_population_tbl_add_row(tsk_population_tbl_t *self,
        const char *metadata, size_t metadata_length);
int tsk_population_tbl_set_columns(tsk_population_tbl_t *self, size_t num_rows,
        const char *metadata, tsk_tbl_size_t *metadata_offset);
int tsk_population_tbl_append_columns(tsk_population_tbl_t *self, size_t num_rows,
        const char *metadata, tsk_tbl_size_t *metadata_offset);
int tsk_population_tbl_clear(tsk_population_tbl_t *self);
int tsk_population_tbl_truncate(tsk_population_tbl_t *self, size_t num_rows);
int tsk_population_tbl_copy(tsk_population_tbl_t *self, tsk_population_tbl_t *dest);
int tsk_population_tbl_free(tsk_population_tbl_t *self);
void tsk_population_tbl_print_state(tsk_population_tbl_t *self, FILE *out);
int tsk_population_tbl_dump_text(tsk_population_tbl_t *self, FILE *out);
bool tsk_population_tbl_equals(tsk_population_tbl_t *self, tsk_population_tbl_t *other);
int tsk_population_tbl_get_row(tsk_population_tbl_t *self, size_t index, tsk_population_t *row);

int tsk_provenance_tbl_alloc(tsk_provenance_tbl_t *self, int flags);
int tsk_provenance_tbl_set_max_rows_increment(tsk_provenance_tbl_t *self, size_t max_rows_increment);
int tsk_provenance_tbl_set_max_timestamp_length_increment(tsk_provenance_tbl_t *self,
        size_t max_timestamp_length_increment);
int tsk_provenance_tbl_set_max_record_length_increment(tsk_provenance_tbl_t *self,
        size_t max_record_length_increment);
tsk_id_t tsk_provenance_tbl_add_row(tsk_provenance_tbl_t *self,
        const char *timestamp, size_t timestamp_length,
        const char *record, size_t record_length);
int tsk_provenance_tbl_set_columns(tsk_provenance_tbl_t *self, size_t num_rows,
       char *timestamp, tsk_tbl_size_t *timestamp_offset,
       char *record, tsk_tbl_size_t *record_offset);
int tsk_provenance_tbl_append_columns(tsk_provenance_tbl_t *self, size_t num_rows,
        char *timestamp, tsk_tbl_size_t *timestamp_offset,
        char *record, tsk_tbl_size_t *record_offset);
int tsk_provenance_tbl_clear(tsk_provenance_tbl_t *self);
int tsk_provenance_tbl_truncate(tsk_provenance_tbl_t *self, size_t num_rows);
int tsk_provenance_tbl_copy(tsk_provenance_tbl_t *self, tsk_provenance_tbl_t *dest);
int tsk_provenance_tbl_free(tsk_provenance_tbl_t *self);
int tsk_provenance_tbl_dump_text(tsk_provenance_tbl_t *self, FILE *out);
void tsk_provenance_tbl_print_state(tsk_provenance_tbl_t *self, FILE *out);
bool tsk_provenance_tbl_equals(tsk_provenance_tbl_t *self, tsk_provenance_tbl_t *other);
int tsk_provenance_tbl_get_row(tsk_provenance_tbl_t *self, size_t index, tsk_provenance_t *row);

/****************************************************************************/
/* Table collection .*/
/****************************************************************************/

int tsk_tbl_collection_alloc(tsk_tbl_collection_t *self, int flags);
int tsk_tbl_collection_load(tsk_tbl_collection_t *self, const char *filename, int flags);
int tsk_tbl_collection_dump(tsk_tbl_collection_t *tables, const char *filename, int flags);
int tsk_tbl_collection_copy(tsk_tbl_collection_t *self, tsk_tbl_collection_t *dest);
int tsk_tbl_collection_print_state(tsk_tbl_collection_t *self, FILE *out);
int tsk_tbl_collection_free(tsk_tbl_collection_t *self);

bool tsk_tbl_collection_is_indexed(tsk_tbl_collection_t *self);
int tsk_tbl_collection_drop_indexes(tsk_tbl_collection_t *self);
int tsk_tbl_collection_build_indexes(tsk_tbl_collection_t *self, int flags);
int tsk_tbl_collection_simplify(tsk_tbl_collection_t *self,
        tsk_id_t *samples, size_t num_samples, int flags, tsk_id_t *node_map);
int tsk_tbl_collection_sort(tsk_tbl_collection_t *self, size_t edge_start, int flags);
int tsk_tbl_collection_deduplicate_sites(tsk_tbl_collection_t *tables, int flags);
int tsk_tbl_collection_compute_mutation_parents(tsk_tbl_collection_t *self, int flags);
bool tsk_tbl_collection_equals(tsk_tbl_collection_t *self, tsk_tbl_collection_t *other);
int tsk_tbl_collection_record_position(tsk_tbl_collection_t *self,
        tsk_tbl_collection_position_t *position);
int tsk_tbl_collection_reset_position(tsk_tbl_collection_t *self,
        tsk_tbl_collection_position_t *position);
int tsk_tbl_collection_clear(tsk_tbl_collection_t *self);
int tsk_tbl_collection_check_integrity(tsk_tbl_collection_t *self, int flags);

int tsk_squash_edges(tsk_edge_t *edges, size_t num_edges, size_t *num_output_edges);

#ifdef __cplusplus
}
#endif
#endif
