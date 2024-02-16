/*
 * MIT License
 *
 * Copyright (c) 2019-2022 Tskit Developers
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, mergetest, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "testlib.h"
#include <tskit/tables.h>

typedef struct {
    const char *name;
    void *array;
    tsk_size_t len;
    int type;
} write_table_col_t;

static void
write_table_cols(kastore_t *store, write_table_col_t *write_cols, size_t num_cols)
{
    size_t j;
    int ret;

    for (j = 0; j < num_cols; j++) {
        ret = kastore_puts(store, write_cols[j].name, write_cols[j].array,
            (size_t) write_cols[j].len, write_cols[j].type, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
copy_store_drop_columns(
    tsk_treeseq_t *ts, size_t num_drop_cols, const char **drop_cols, const char *outfile)
{
    int ret = 0;
    char tmpfile[] = "/tmp/tsk_c_test_copy_XXXXXX";
    int fd;
    kastore_t read_store, write_store;
    kaitem_t *item;
    size_t j, k;
    bool keep;

    fd = mkstemp(tmpfile);
    CU_ASSERT_FATAL(fd != -1);
    close(fd);

    ret = tsk_treeseq_dump(ts, tmpfile, 0);
    if (ret != 0) {
        unlink(tmpfile);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }

    ret = kastore_open(&read_store, tmpfile, "r", KAS_READ_ALL);
    /* We can now unlink the file as either kastore has read it all, or failed */
    unlink(tmpfile);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = kastore_open(&write_store, outfile, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Note: this API is not a documented part of kastore, so may be subject to
     * change. */
    for (j = 0; j < read_store.num_items; j++) {
        item = &read_store.items[j];
        keep = true;
        for (k = 0; k < num_drop_cols; k++) {
            if (strlen(drop_cols[k]) == item->key_len
                && strncmp(drop_cols[k], item->key, item->key_len) == 0) {
                keep = false;
                break;
            }
        }
        if (keep) {
            ret = kastore_put(&write_store, item->key, item->key_len, item->array,
                item->array_len, item->type, 0);
            CU_ASSERT_EQUAL_FATAL(ret, 0);
        }
    }
    kastore_close(&read_store);
    ret = kastore_close(&write_store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
}

static void
test_format_data_load_errors(void)
{
    size_t uuid_size = 36;
    char uuid[uuid_size];
    char format_name[TSK_FILE_FORMAT_NAME_LENGTH];
    double L[2];
    uint32_t version[2]
        = { TSK_FILE_FORMAT_VERSION_MAJOR, TSK_FILE_FORMAT_VERSION_MINOR };
    write_table_col_t write_cols[] = {
        { "format/name", (void *) format_name, sizeof(format_name), KAS_INT8 },
        { "format/version", (void *) version, 2, KAS_UINT32 },
        { "sequence_length", (void *) L, 1, KAS_FLOAT64 },
        { "uuid", (void *) uuid, (tsk_size_t) uuid_size, KAS_INT8 },
    };
    tsk_table_collection_t tables;
    kastore_t store;
    size_t j;
    int ret;

    L[0] = 1;
    L[1] = 0;
    tsk_memcpy(format_name, TSK_FILE_FORMAT_NAME, sizeof(format_name));
    /* Note: this will fail if we ever start parsing the form of the UUID */
    tsk_memset(uuid, 0, uuid_size);

    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    /* We've only defined the format headers, so we should fail immediately
     * after with required columns not found */
    CU_ASSERT_FALSE(tsk_is_kas_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_REQUIRED_COL_NOT_FOUND);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Version too old */
    version[0] = TSK_FILE_FORMAT_VERSION_MAJOR - 1;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_VERSION_TOO_OLD);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Version too new */
    version[0] = TSK_FILE_FORMAT_VERSION_MAJOR + 1;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_VERSION_TOO_NEW);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    version[0] = TSK_FILE_FORMAT_VERSION_MAJOR;

    /* Bad version length */
    write_cols[1].len = 0;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[1].len = 2;

    /* Bad format name length */
    write_cols[0].len = 0;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[0].len = TSK_FILE_FORMAT_NAME_LENGTH;

    /* Bad format name */
    format_name[0] = 'X';
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    format_name[0] = 't';

    /* Bad type for sequence length. */
    write_cols[2].type = KAS_FLOAT32;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_TRUE(tsk_is_kas_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret ^ (1 << TSK_KAS_ERR_BIT), KAS_ERR_TYPE_MISMATCH);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[2].type = KAS_FLOAT64;

    /* Bad length for sequence length. */
    write_cols[2].len = 2;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[2].len = 1;

    /* Bad value for sequence length. */
    L[0] = -1;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_SEQUENCE_LENGTH);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    L[0] = 1;

    /* Wrong length for uuid */
    write_cols[3].len = 1;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[3].len = (tsk_size_t) uuid_size;

    /* Missing keys */
    for (j = 0; j < sizeof(write_cols) / sizeof(*write_cols) - 1; j++) {
        ret = kastore_open(&store, _tmp_file_name, "w", 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        write_table_cols(&store, write_cols, j);
        ret = kastore_close(&store);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_REQUIRED_COL_NOT_FOUND);
        ret = tsk_table_collection_free(&tables);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
    }
}

static void
test_missing_optional_column_pairs(void)
{
    int ret;
    size_t j;
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    tsk_table_collection_t t1, t2;
    const char *required_cols[][2] = { { "edges/metadata", "edges/metadata_offset" },
        { "migrations/metadata", "migrations/metadata_offset" },
        { "individuals/parents", "individuals/parents_offset" } };
    const char *drop_cols[2];

    ret = tsk_treeseq_copy_tables(ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < sizeof(required_cols) / sizeof(*required_cols); j++) {
        drop_cols[0] = required_cols[j][0];
        copy_store_drop_columns(ts, 1, drop_cols, _tmp_file_name);
        ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BOTH_COLUMNS_REQUIRED);
        tsk_table_collection_free(&t2);

        drop_cols[0] = required_cols[j][1];
        copy_store_drop_columns(ts, 1, drop_cols, _tmp_file_name);
        ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BOTH_COLUMNS_REQUIRED);
        tsk_table_collection_free(&t2);

        drop_cols[0] = required_cols[j][0];
        drop_cols[1] = required_cols[j][1];
        copy_store_drop_columns(ts, 2, drop_cols, _tmp_file_name);
        ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_FALSE(tsk_table_collection_equals(&t1, &t2, 0));
        tsk_table_collection_free(&t2);
    }

    tsk_table_collection_free(&t1);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_missing_required_column_pairs(void)
{
    int ret;
    size_t j;
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    tsk_table_collection_t t;
    const char *required_cols[][2] = {
        { "individuals/location", "individuals/location_offset" },
        { "individuals/metadata", "individuals/metadata_offset" },
        { "mutations/derived_state", "mutations/derived_state_offset" },
        { "mutations/metadata", "mutations/metadata_offset" },
        { "nodes/metadata", "nodes/metadata_offset" },
        { "populations/metadata", "populations/metadata_offset" },
        { "provenances/record", "provenances/record_offset" },
        { "provenances/timestamp", "provenances/timestamp_offset" },
        { "sites/ancestral_state", "sites/ancestral_state_offset" },
        { "sites/metadata", "sites/metadata_offset" },
    };
    const char *drop_cols[2];

    for (j = 0; j < sizeof(required_cols) / sizeof(*required_cols); j++) {
        drop_cols[0] = required_cols[j][0];
        copy_store_drop_columns(ts, 1, drop_cols, _tmp_file_name);
        ret = tsk_table_collection_load(&t, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_REQUIRED_COL_NOT_FOUND);
        tsk_table_collection_free(&t);

        drop_cols[0] = required_cols[j][1];
        copy_store_drop_columns(ts, 1, drop_cols, _tmp_file_name);
        ret = tsk_table_collection_load(&t, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BOTH_COLUMNS_REQUIRED);
        tsk_table_collection_free(&t);

        copy_store_drop_columns(ts, 2, required_cols[j], _tmp_file_name);
        ret = tsk_table_collection_load(&t, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_REQUIRED_COL_NOT_FOUND);
        tsk_table_collection_free(&t);
    }

    tsk_treeseq_free(ts);
    free(ts);
}

static void
verify_bad_offset_columns(tsk_treeseq_t *ts, const char *offset_col)
{
    int ret = 0;
    kastore_t store;
    tsk_table_collection_t tables;
    uint32_t *offset_array, *offset_copy;
    size_t offset_len;
    int type;
    uint32_t data_len;

    ret = tsk_treeseq_dump(ts, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_open(&store, _tmp_file_name, "r", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = kastore_gets(&store, offset_col, (void **) &offset_array, &offset_len, &type);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(type, KAS_UINT32);
    offset_copy = malloc(offset_len * sizeof(*offset_array));
    CU_ASSERT_FATAL(offset_copy != NULL);
    tsk_memcpy(offset_copy, offset_array, offset_len * sizeof(*offset_array));
    data_len = offset_array[offset_len - 1];
    CU_ASSERT_TRUE(data_len > 0);
    kastore_close(&store);

    offset_copy[0] = UINT32_MAX;
    copy_store_drop_columns(ts, 1, &offset_col, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, offset_col, offset_copy, offset_len, KAS_UINT32, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_OFFSET);
    tsk_table_collection_free(&tables);

    offset_copy[0] = 0;
    offset_copy[offset_len - 1] = 0;
    copy_store_drop_columns(ts, 1, &offset_col, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, offset_col, offset_copy, offset_len, KAS_UINT32, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_OFFSET);
    tsk_table_collection_free(&tables);

    offset_copy[offset_len - 1] = data_len + 1;
    copy_store_drop_columns(ts, 1, &offset_col, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, offset_col, offset_copy, offset_len, KAS_UINT32, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_OFFSET);
    tsk_table_collection_free(&tables);

    copy_store_drop_columns(ts, 1, &offset_col, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, offset_col, NULL, 0, KAS_UINT32, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    tsk_table_collection_free(&tables);

    copy_store_drop_columns(ts, 1, &offset_col, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, offset_col, offset_copy, offset_len, KAS_FLOAT32, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_COLUMN_TYPE);
    tsk_table_collection_free(&tables);

    free(offset_copy);
}

static void
test_bad_offset_columns(void)
{
    size_t j;
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    /* We exclude "provenances/timestamp_offset" here because there are no
     * non-ragged columns in the provenances table, so this doesn't quite
     * fit into the same pattern as the other tables */
    const char *cols[] = {
        "edges/metadata_offset",
        "migrations/metadata_offset",
        "individuals/location_offset",
        "individuals/parents_offset",
        "individuals/metadata_offset",
        "mutations/derived_state_offset",
        "mutations/metadata_offset",
        "nodes/metadata_offset",
        "populations/metadata_offset",
        "provenances/record_offset",
        "sites/ancestral_state_offset",
        "sites/metadata_offset",
    };

    for (j = 0; j < sizeof(cols) / sizeof(*cols); j++) {
        verify_bad_offset_columns(ts, cols[j]);
    }
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_force_offset_64(void)
{
    int ret;
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    tsk_table_collection_t t1;
    tsk_table_collection_t t2;
    kastore_t store;
    kaitem_t *item;
    const char *suffix;
    const char *offset_str = "_offset";
    int num_found = 0;
    size_t j;

    ret = tsk_treeseq_dump(ts, _tmp_file_name, TSK_DUMP_FORCE_OFFSET_64);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = kastore_open(&store, _tmp_file_name, "r", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    for (j = 0; j < store.num_items; j++) {
        item = &store.items[j];
        /* Does the key end in "_offset"? */
        if (item->key_len > strlen(offset_str)) {
            suffix = item->key + (item->key_len - strlen(offset_str));
            if (strncmp(suffix, offset_str, strlen(offset_str)) == 0) {
                CU_ASSERT_EQUAL(item->type, KAS_UINT64);
                num_found++;
            }
        }
    }
    CU_ASSERT_TRUE(num_found > 0);
    kastore_close(&store);

    ret = tsk_table_collection_load(&t1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    ret = tsk_treeseq_copy_tables(ts, &t2, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));

    tsk_table_collection_free(&t1);
    tsk_table_collection_free(&t2);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_missing_indexes(void)
{
    int ret;
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    tsk_table_collection_t t1, t2;
    const char *cols[]
        = { "indexes/edge_insertion_order", "indexes/edge_removal_order" };
    const char *drop_cols[2];

    ret = tsk_treeseq_copy_tables(ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    drop_cols[0] = cols[0];
    copy_store_drop_columns(ts, 1, drop_cols, _tmp_file_name);
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BOTH_COLUMNS_REQUIRED);
    tsk_table_collection_free(&t2);

    drop_cols[0] = cols[1];
    copy_store_drop_columns(ts, 1, drop_cols, _tmp_file_name);
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BOTH_COLUMNS_REQUIRED);
    tsk_table_collection_free(&t2);

    copy_store_drop_columns(ts, 2, cols, _tmp_file_name);
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    CU_ASSERT_FALSE(tsk_table_collection_has_index(&t2, 0));
    tsk_table_collection_free(&t2);

    tsk_table_collection_free(&t1);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_malformed_indexes(void)
{
    int ret;
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    tsk_table_collection_t tables;
    tsk_treeseq_t ts2;
    tsk_size_t num_edges = tsk_treeseq_get_num_edges(ts);
    tsk_id_t *bad_index = tsk_calloc(num_edges, sizeof(tsk_id_t));
    tsk_id_t *good_index = tsk_calloc(num_edges, sizeof(tsk_id_t));
    kastore_t store;
    const char *cols[]
        = { "indexes/edge_insertion_order", "indexes/edge_removal_order" };

    CU_ASSERT_FATAL(bad_index != NULL);
    CU_ASSERT_FATAL(good_index != NULL);

    /* If both columns are not the same length as the number of edges we
     * should raise an error */
    copy_store_drop_columns(ts, 2, cols, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, cols[0], NULL, 0, TSK_ID_STORAGE_TYPE, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, cols[1], NULL, 0, TSK_ID_STORAGE_TYPE, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    tsk_table_collection_free(&tables);

    bad_index[0] = -1;

    copy_store_drop_columns(ts, 2, cols, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(
        &store, cols[0], good_index, (size_t) num_edges, TSK_ID_STORAGE_TYPE, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(
        &store, cols[1], bad_index, (size_t) num_edges, TSK_ID_STORAGE_TYPE, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_load(&ts2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts2);

    copy_store_drop_columns(ts, 2, cols, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(
        &store, cols[0], bad_index, (size_t) num_edges, TSK_ID_STORAGE_TYPE, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(
        &store, cols[1], good_index, (size_t) num_edges, TSK_ID_STORAGE_TYPE, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_load(&ts2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EDGE_OUT_OF_BOUNDS);
    tsk_treeseq_free(&ts2);

    copy_store_drop_columns(ts, 1, cols, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, cols[0], bad_index, (size_t) num_edges, KAS_FLOAT32, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_load(&ts2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_COLUMN_TYPE);
    tsk_treeseq_free(&ts2);

    free(good_index);
    free(bad_index);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_missing_reference_sequence(void)
{
    int ret;
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    tsk_table_collection_t t1, t2;
    const char *cols[] = { "reference_sequence/data", "reference_sequence/url",
        "reference_sequence/metadata_schema", "reference_sequence/metadata" };

    CU_ASSERT_TRUE(tsk_treeseq_has_reference_sequence(ts));

    ret = tsk_treeseq_copy_tables(ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    copy_store_drop_columns(ts, 1, cols, _tmp_file_name);
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_TRUE(tsk_table_collection_has_reference_sequence(&t2));
    tsk_table_collection_free(&t2);

    copy_store_drop_columns(ts, 2, cols, _tmp_file_name);
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_TRUE(tsk_table_collection_has_reference_sequence(&t2));
    tsk_table_collection_free(&t2);

    copy_store_drop_columns(ts, 3, cols, _tmp_file_name);
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_TRUE(tsk_table_collection_has_reference_sequence(&t2));
    tsk_table_collection_free(&t2);

    /* Dropping all the columns gives us a NULL reference_sequence, though */
    copy_store_drop_columns(ts, 4, cols, _tmp_file_name);
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_FALSE(tsk_table_collection_has_reference_sequence(&t2));
    tsk_table_collection_free(&t2);

    tsk_table_collection_free(&t1);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_bad_column_types(void)
{
    int ret;
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    tsk_table_collection_t tables;
    tsk_size_t num_edges = tsk_treeseq_get_num_edges(ts);
    /* make sure we have enough memory in all cases */
    tsk_id_t *col_memory = tsk_calloc(num_edges + 1, sizeof(double));
    kastore_t store;
    const char *cols[1];

    CU_ASSERT_FATAL(col_memory != NULL);

    cols[0] = "edges/left";
    copy_store_drop_columns(ts, 1, cols, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, cols[0], col_memory, (size_t) num_edges, KAS_FLOAT32, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_COLUMN_TYPE);
    tsk_table_collection_free(&tables);

    cols[0] = "edges/metadata_offset";
    copy_store_drop_columns(ts, 1, cols, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(
        &store, cols[0], col_memory, (size_t) num_edges + 1, KAS_FLOAT32, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_COLUMN_TYPE);
    tsk_table_collection_free(&tables);

    cols[0] = "edges/metadata";
    copy_store_drop_columns(ts, 1, cols, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, cols[0], NULL, 0, KAS_FLOAT32, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_COLUMN_TYPE);
    tsk_table_collection_free(&tables);

    cols[0] = "edges/metadata_schema";
    copy_store_drop_columns(ts, 1, cols, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, cols[0], NULL, 0, KAS_FLOAT32, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_COLUMN_TYPE);
    tsk_table_collection_free(&tables);

    cols[0] = "reference_sequence/metadata";
    copy_store_drop_columns(ts, 1, cols, _tmp_file_name);
    ret = kastore_open(&store, _tmp_file_name, "a", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_puts(&store, cols[0], NULL, 0, KAS_FLOAT32, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_COLUMN_TYPE);
    tsk_table_collection_free(&tables);

    free(col_memory);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_missing_required_columns(void)
{
    int ret;
    size_t j;
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    tsk_table_collection_t t;
    const char *required_cols[] = {
        "edges/child",
        "edges/left",
        "edges/parent",
        "edges/right",
        "format/name",
        "format/version",
        "individuals/flags",
        "migrations/dest",
        "migrations/left",
        "migrations/node",
        "migrations/right",
        "migrations/source",
        "migrations/time",
        "mutations/node",
        "mutations/parent",
        "mutations/site",
        "nodes/flags",
        "nodes/individual",
        "nodes/population",
        "nodes/time",
        "sequence_length",
        "sites/position",
        "uuid",
    };
    const char *drop_cols[1];

    for (j = 0; j < sizeof(required_cols) / sizeof(*required_cols); j++) {
        drop_cols[0] = required_cols[j];
        copy_store_drop_columns(ts, 1, drop_cols, _tmp_file_name);
        ret = tsk_table_collection_load(&t, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_REQUIRED_COL_NOT_FOUND);
        tsk_table_collection_free(&t);
    }

    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_metadata_schemas_optional(void)
{
    int ret;
    size_t j;
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    tsk_table_collection_t t1, t2;
    const char *cols[] = {
        "metadata",
        "metadata_schema",
        "reference_sequence/metadata",
        "reference_sequence/metadata_schema",
        "individuals/metadata_schema",
        "populations/metadata_schema",
        "nodes/metadata_schema",
        "edges/metadata_schema",
        "sites/metadata_schema",
        "mutations/metadata_schema",
        "migrations/metadata_schema",
    };
    const char *drop_cols[1];

    ret = tsk_treeseq_copy_tables(ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    for (j = 0; j < sizeof(cols) / sizeof(*cols); j++) {
        drop_cols[0] = cols[j];
        copy_store_drop_columns(ts, 1, drop_cols, _tmp_file_name);
        ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        /* metadata schemas are included in data comparisons */
        CU_ASSERT_FALSE(tsk_table_collection_equals(&t1, &t2, 0));
        tsk_table_collection_free(&t2);
    }

    tsk_table_collection_free(&t1);
    tsk_treeseq_free(ts);
    free(ts);
}

/* This test is problematic on windows because of the different off_t
 * types. Doesn't seem worth the trouble of getting it working.
 */
static void
test_load_bad_file_formats(void)
{
#if !defined(_WIN32)
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    int ret, ret2;
    off_t offset;
    FILE *f;

    /* A zero byte file is TSK_ERR_EOF */
    f = fopen(_tmp_file_name, "w+");
    ret = tsk_table_collection_loadf(&tables, f, 0);
    ret2 = tsk_treeseq_loadf(&ts, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, ret2);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EOF);
    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
    fclose(f);

    for (offset = 1; offset < 100; offset++) {
        ret = tsk_table_collection_init(&tables, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tables.sequence_length = 1.0;
        ret = tsk_table_collection_dump(&tables, _tmp_file_name, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);

        truncate(_tmp_file_name, offset);
        ret = tsk_table_collection_load(&tables, _tmp_file_name, TSK_NO_INIT);
        CU_ASSERT_EQUAL_FATAL(ret ^ (1 << TSK_KAS_ERR_BIT), KAS_ERR_BAD_FILE_FORMAT);
        tsk_table_collection_free(&tables);
    }
#endif
}

static void
test_load_errors(void)
{
    tsk_table_collection_t tables;
    tsk_treeseq_t ts;
    int ret, ret2;
    const char *str;
    FILE *f;

    ret = tsk_table_collection_load(&tables, "/", 0);
    ret2 = tsk_treeseq_load(&ts, "/", 0);
    CU_ASSERT_EQUAL_FATAL(ret, ret2);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_IO);
    str = tsk_strerror(ret);
    CU_ASSERT_TRUE(strlen(str) > 0);
    CU_ASSERT_STRING_EQUAL(str, strerror(EISDIR));
    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);

    ret = tsk_table_collection_load(&tables, "/bin/theres_no_way_this_file_exists", 0);
    ret2 = tsk_treeseq_load(&ts, "/bin/theres_no_way_this_file_exists", 0);
    CU_ASSERT_EQUAL_FATAL(ret, ret2);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_IO);
    str = tsk_strerror(ret);
    CU_ASSERT_TRUE(strlen(str) > 0);
    CU_ASSERT_STRING_EQUAL(str, strerror(ENOENT));
    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);

    ret = tsk_table_collection_load(&tables, "/bin/sh", 0);
    ret2 = tsk_treeseq_load(&ts, "/bin/sh", 0);
    CU_ASSERT_EQUAL_FATAL(ret, ret2);
    CU_ASSERT_TRUE(tsk_is_kas_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret ^ (1 << TSK_KAS_ERR_BIT), KAS_ERR_BAD_FILE_FORMAT);
    str = tsk_strerror(ret);
    CU_ASSERT_TRUE(strlen(str) > 0);
    tsk_table_collection_free(&tables);

    /* open a file in the wrong mode */
    f = fopen(_tmp_file_name, "w");
    ret = tsk_table_collection_loadf(&tables, f, 0);
    ret2 = tsk_treeseq_loadf(&ts, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, ret2);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_IO);
    str = tsk_strerror(ret);
    CU_ASSERT_TRUE(strlen(str) > 0);
    CU_ASSERT_STRING_EQUAL(str, strerror(EBADF));
    tsk_table_collection_free(&tables);
    tsk_treeseq_free(&ts);
    fclose(f);
}

static void
test_load_eof(void)
{
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    tsk_table_collection_t tables;
    int ret;
    FILE *f;

    f = fopen(_tmp_file_name, "w+");
    CU_ASSERT_NOT_EQUAL(f, NULL);
    ret = tsk_table_collection_loadf(&tables, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EOF);
    fclose(f);
    tsk_table_collection_free(&tables);

    /* Reading an empty file also returns EOF */
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EOF);
    tsk_table_collection_free(&tables);

    f = fopen(_tmp_file_name, "w+");
    CU_ASSERT_NOT_EQUAL(f, NULL);
    ret = tsk_treeseq_dumpf(ts, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Reading from the end of the stream gives EOF */
    ret = tsk_table_collection_loadf(&tables, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EOF);
    tsk_table_collection_free(&tables);

    /* Reading the start of the stream is fine */
    fseek(f, 0, SEEK_SET);
    ret = tsk_table_collection_loadf(&tables, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_table_collection_free(&tables);

    /* And we should be back to the end of the stream */
    ret = tsk_table_collection_loadf(&tables, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EOF);
    tsk_table_collection_free(&tables);

    /* Trying to read the same end stream should give the same
     * result. */
    ret = tsk_table_collection_loadf(&tables, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EOF);
    tsk_table_collection_free(&tables);

    /* A previously init'd tables should be good too */
    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_loadf(&tables, f, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_EOF);
    tsk_table_collection_free(&tables);

    fclose(f);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_dump_errors(void)
{
    tsk_table_collection_t tables;
    int ret;
    FILE *f;
    const char *str;

    ret = tsk_table_collection_init(&tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tables.sequence_length = 1.0;

    ret = tsk_table_collection_dump(&tables, "/", 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_IO);
    str = tsk_strerror(ret);
    CU_ASSERT_TRUE(strlen(str) > 0);
    CU_ASSERT_STRING_EQUAL(str, strerror(EISDIR));

    /* We're assuming that we don't have write access to /bin, so don't run this
     * as root! */
    ret = tsk_table_collection_dump(&tables, "/bin/theres_no_way_this_file_exists", 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_IO);
    str = tsk_strerror(ret);
    CU_ASSERT_TRUE(strlen(str) > 0);
    CU_ASSERT_TRUE(
        (strcmp(str, strerror(EACCES)) == 0) || (strcmp(str, strerror(EPERM)) == 0));

    /* open a file in the wrong mode */
    f = fopen(_tmp_file_name, "r");
    ret = tsk_table_collection_dumpf(&tables, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_IO);
    str = tsk_strerror(ret);
    CU_ASSERT_TRUE(strlen(str) > 0);
    CU_ASSERT_STRING_EQUAL(str, strerror(EBADF));
    fclose(f);

    /* We'd like to catch close errors also, but it's hard to provoke them
     * without intercepting calls to fclose() */

    tsk_table_collection_free(&tables);
}

/* FIXME these are good tests, but we want to make them more general so that
 * they can be applied to other tables.*/
static void
test_load_node_table_errors(void)
{
    char format_name[TSK_FILE_FORMAT_NAME_LENGTH];
    size_t uuid_size = 36;
    char uuid[uuid_size];
    double L = 1;
    double time = 0;
    double flags = 0;
    tsk_id_t population = 0;
    tsk_id_t individual = 0;
    int8_t metadata = 0;
    uint32_t metadata_offset[] = { 0, 1 };
    uint32_t version[2]
        = { TSK_FILE_FORMAT_VERSION_MAJOR, TSK_FILE_FORMAT_VERSION_MINOR };
    write_table_col_t write_cols[] = {
        { "nodes/time", (void *) &time, 1, KAS_FLOAT64 },
        { "nodes/flags", (void *) &flags, 1, TSK_FLAGS_STORAGE_TYPE },
        { "nodes/population", (void *) &population, 1, TSK_ID_STORAGE_TYPE },
        { "nodes/individual", (void *) &individual, 1, TSK_ID_STORAGE_TYPE },
        { "nodes/metadata", (void *) &metadata, 1, KAS_UINT8 },
        { "nodes/metadata_offset", (void *) metadata_offset, 2, KAS_UINT32 },
        { "format/name", (void *) format_name, sizeof(format_name), KAS_INT8 },
        { "format/version", (void *) version, 2, KAS_UINT32 },
        { "uuid", (void *) uuid, uuid_size, KAS_INT8 },
        { "sequence_length", (void *) &L, 1, KAS_FLOAT64 },
    };
    tsk_table_collection_t tables;
    kastore_t store;
    int ret;

    tsk_memcpy(format_name, TSK_FILE_FORMAT_NAME, sizeof(format_name));
    /* Note: this will fail if we ever start parsing the form of the UUID */
    tsk_memset(uuid, 0, uuid_size);

    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    /* We've only defined the format headers and nodes, so we should fail immediately
     * after with key not found */
    CU_ASSERT_FALSE(tsk_is_kas_error(ret));
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_REQUIRED_COL_NOT_FOUND);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    /* Wrong type for time */
    write_cols[0].type = KAS_INT64;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_BAD_COLUMN_TYPE);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[0].type = KAS_FLOAT64;

    /* Wrong length for flags */
    write_cols[1].len = 0;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[1].len = 1;

    /* Wrong length for metadata offset */
    write_cols[5].len = 1;
    ret = kastore_open(&store, _tmp_file_name, "w", 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_table_cols(&store, write_cols, sizeof(write_cols) / sizeof(*write_cols));
    ret = kastore_close(&store);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&tables, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, TSK_ERR_FILE_FORMAT);
    ret = tsk_table_collection_free(&tables);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    write_cols[5].len = 2;
}

static void
test_example_round_trip(void)
{
    int ret;
    tsk_treeseq_t *ts1 = caterpillar_tree(5, 3, 3);
    tsk_treeseq_t ts2;
    tsk_table_collection_t t1, t2;
    FILE *f;

    ret = tsk_treeseq_copy_tables(ts1, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_dump(&t1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));

    /* Reading multiple times into the same tables with TSK_NO_INIT is supported. */
    ret = tsk_table_collection_load(&t2, _tmp_file_name, TSK_NO_INIT);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    tsk_table_collection_free(&t2);

    /* Do the same thing with treeseq API */
    remove(_tmp_file_name);
    ret = tsk_treeseq_dump(ts1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_treeseq_load(&ts2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, ts2.tables, 0));
    tsk_treeseq_free(&ts2);

    /* Use loadf form */
    f = fopen(_tmp_file_name, "w+");
    ret = tsk_table_collection_dumpf(&t1, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    fseek(f, 0, SEEK_SET);
    ret = tsk_table_collection_loadf(&t2, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    tsk_table_collection_free(&t2);
    fclose(f);

    /* Do the same thing with treeseq API */
    f = fopen(_tmp_file_name, "w+");
    ret = tsk_treeseq_dumpf(ts1, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    fseek(f, 0, SEEK_SET);
    ret = tsk_treeseq_loadf(&ts2, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, ts2.tables, 0));
    tsk_treeseq_free(&ts2);

    fclose(f);
    tsk_table_collection_free(&t1);
    tsk_treeseq_free(ts1);
    free(ts1);
}

static void
test_multiple_round_trip(void)
{
    int ret;
    tsk_size_t j;
    tsk_size_t num_examples = 10;
    tsk_treeseq_t *ts;
    tsk_table_collection_t in_tables[num_examples];
    tsk_table_collection_t out_tables;
    FILE *f = fopen(_tmp_file_name, "w+");

    CU_ASSERT_NOT_EQUAL_FATAL(f, NULL);

    for (j = 0; j < num_examples; j++) {
        ts = caterpillar_tree(5 + j, 3 + j, 3 + j);
        ret = tsk_treeseq_copy_tables(ts, &in_tables[j], 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        ret = tsk_treeseq_dumpf(ts, f, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        tsk_treeseq_free(ts);
        free(ts);
    }

    fseek(f, 0, SEEK_SET);
    for (j = 0; j < num_examples; j++) {
        ret = tsk_table_collection_loadf(&out_tables, f, 0);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_TRUE(tsk_table_collection_equals(&in_tables[j], &out_tables, 0));
        tsk_table_collection_free(&out_tables);
    }

    /* Can do the same with the same set of previously init'd tables. */
    ret = tsk_table_collection_init(&out_tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    fseek(f, 0, SEEK_SET);
    for (j = 0; j < num_examples; j++) {
        ret = tsk_table_collection_loadf(&out_tables, f, TSK_NO_INIT);
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_TRUE(tsk_table_collection_equals(&in_tables[j], &out_tables, 0));
    }
    tsk_table_collection_free(&out_tables);

    /* Can also read until EOF to do the same thing */
    ret = tsk_table_collection_init(&out_tables, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    fseek(f, 0, SEEK_SET);
    j = 0;
    while (true) {
        ret = tsk_table_collection_loadf(&out_tables, f, TSK_NO_INIT);
        if (ret == TSK_ERR_EOF) {
            break;
        }
        CU_ASSERT_EQUAL_FATAL(ret, 0);
        CU_ASSERT_TRUE(tsk_table_collection_equals(&in_tables[j], &out_tables, 0));
        j++;
    }
    tsk_table_collection_free(&out_tables);
    CU_ASSERT_EQUAL_FATAL(j, num_examples);

    for (j = 0; j < num_examples; j++) {
        tsk_table_collection_free(&in_tables[j]);
    }
    fclose(f);
}

static void
test_copy_store_drop_columns(void)
{
    int ret;
    tsk_treeseq_t *ts = caterpillar_tree(5, 3, 3);
    tsk_table_collection_t t1, t2;

    ret = tsk_treeseq_copy_tables(ts, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    /* Dropping no columns should have no effect on the data */
    copy_store_drop_columns(ts, 0, NULL, _tmp_file_name);
    ret = tsk_table_collection_load(&t2, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));

    tsk_table_collection_free(&t1);
    tsk_table_collection_free(&t2);
    tsk_treeseq_free(ts);
    free(ts);
}

static void
test_skip_tables(void)
{
    int ret;
    tsk_treeseq_t *ts1 = caterpillar_tree(5, 3, 3);
    tsk_treeseq_t ts2;
    tsk_table_collection_t t1, t2;
    FILE *f;

    ret = tsk_treeseq_dump(ts1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(&t1, _tmp_file_name, TSK_LOAD_SKIP_TABLES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, ts1->tables, TSK_CMP_IGNORE_TABLES));
    CU_ASSERT_EQUAL(t1.individuals.num_rows, 0);
    CU_ASSERT_EQUAL(t1.nodes.num_rows, 0);
    CU_ASSERT_EQUAL(t1.edges.num_rows, 0);
    CU_ASSERT_EQUAL(t1.migrations.num_rows, 0);
    CU_ASSERT_EQUAL(t1.sites.num_rows, 0);
    CU_ASSERT_EQUAL(t1.mutations.num_rows, 0);
    CU_ASSERT_EQUAL(t1.provenances.num_rows, 0);

    /* Test _loadf code path as well */
    f = fopen(_tmp_file_name, "r+");
    ret = tsk_table_collection_loadf(&t2, f, TSK_LOAD_SKIP_TABLES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    fclose(f);
    tsk_table_collection_free(&t2);

    /* Without TSK_LOAD_SKIP_TABLES we reach end of file */
    f = fopen(_tmp_file_name, "r+");
    ret = tsk_table_collection_loadf(&t2, f, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL(fgetc(f), EOF);
    fclose(f);
    tsk_table_collection_free(&t2);

    /* Setting TSK_LOAD_SKIP_TABLES only reads part of the file */
    f = fopen(_tmp_file_name, "r+");
    ret = tsk_table_collection_loadf(&t2, f, TSK_LOAD_SKIP_TABLES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_NOT_EQUAL(fgetc(f), EOF);
    fclose(f);
    tsk_table_collection_free(&t2);

    /* We should be able to make a tree sequence */
    ret = tsk_treeseq_init(&ts2, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts2);

    /* Do the same thing with treeseq API */
    ret = tsk_treeseq_load(&ts2, _tmp_file_name, TSK_LOAD_SKIP_TABLES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, ts2.tables, 0));
    tsk_treeseq_free(&ts2);

    f = fopen(_tmp_file_name, "r+");
    ret = tsk_treeseq_loadf(&ts2, f, TSK_LOAD_SKIP_TABLES);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, ts2.tables, 0));
    fclose(f);
    tsk_treeseq_free(&ts2);

    tsk_table_collection_free(&t1);
    tsk_treeseq_free(ts1);
    free(ts1);
}

static void
test_skip_reference_sequence(void)
{
    int ret;
    tsk_treeseq_t *ts1 = caterpillar_tree(5, 3, 3);
    tsk_treeseq_t ts2;
    tsk_table_collection_t t1, t2;
    FILE *f;

    CU_ASSERT_TRUE(tsk_treeseq_has_reference_sequence(ts1));

    ret = tsk_treeseq_dump(ts1, _tmp_file_name, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    ret = tsk_table_collection_load(
        &t1, _tmp_file_name, TSK_LOAD_SKIP_REFERENCE_SEQUENCE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);

    CU_ASSERT_FALSE(tsk_table_collection_equals(&t1, ts1->tables, 0));
    CU_ASSERT_TRUE(tsk_table_collection_equals(
        &t1, ts1->tables, TSK_CMP_IGNORE_REFERENCE_SEQUENCE));
    CU_ASSERT_FALSE(tsk_table_collection_has_reference_sequence(&t1));

    /* Test _loadf code path as well */
    f = fopen(_tmp_file_name, "r+");
    ret = tsk_table_collection_loadf(&t2, f, TSK_LOAD_SKIP_REFERENCE_SEQUENCE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, &t2, 0));
    fclose(f);
    tsk_table_collection_free(&t2);

    /* Setting TSK_LOAD_SKIP_REFERENCE_SEQUENCE only reads part of the file */
    f = fopen(_tmp_file_name, "r+");
    ret = tsk_table_collection_loadf(&t2, f, TSK_LOAD_SKIP_REFERENCE_SEQUENCE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_NOT_EQUAL(fgetc(f), EOF);
    fclose(f);
    tsk_table_collection_free(&t2);

    /* We should be able to make a tree sequence */
    ret = tsk_treeseq_init(&ts2, &t1, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    tsk_treeseq_free(&ts2);

    /* Do the same thing with treeseq API */
    ret = tsk_treeseq_load(&ts2, _tmp_file_name, TSK_LOAD_SKIP_REFERENCE_SEQUENCE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, ts2.tables, 0));
    tsk_treeseq_free(&ts2);

    f = fopen(_tmp_file_name, "r+");
    ret = tsk_treeseq_loadf(&ts2, f, TSK_LOAD_SKIP_REFERENCE_SEQUENCE);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_TRUE(tsk_table_collection_equals(&t1, ts2.tables, 0));
    fclose(f);
    tsk_treeseq_free(&ts2);

    tsk_table_collection_free(&t1);
    tsk_treeseq_free(ts1);
    free(ts1);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        { "test_format_data_load_errors", test_format_data_load_errors },
        { "test_missing_indexes", test_missing_indexes },
        { "test_malformed_indexes", test_malformed_indexes },
        { "test_missing_reference_sequence", test_missing_reference_sequence },
        { "test_bad_column_types", test_bad_column_types },
        { "test_missing_required_columns", test_missing_required_columns },
        { "test_missing_optional_column_pairs", test_missing_optional_column_pairs },
        { "test_missing_required_column_pairs", test_missing_required_column_pairs },
        { "test_bad_offset_columns", test_bad_offset_columns },
        { "test_force_offset_64", test_force_offset_64 },
        { "test_metadata_schemas_optional", test_metadata_schemas_optional },
        { "test_load_node_table_errors", test_load_node_table_errors },
        { "test_load_bad_file_formats", test_load_bad_file_formats },
        { "test_load_errors", test_load_errors },
        { "test_load_eof", test_load_eof },
        { "test_dump_errors", test_dump_errors },
        { "test_example_round_trip", test_example_round_trip },
        { "test_multiple_round_trip", test_multiple_round_trip },
        { "test_copy_store_drop_columns", test_copy_store_drop_columns },
        { "test_skip_tables", test_skip_tables },
        { "test_skip_reference_sequence", test_skip_reference_sequence },
        { NULL, NULL },
    };

    return test_main(tests, argc, argv);
}
