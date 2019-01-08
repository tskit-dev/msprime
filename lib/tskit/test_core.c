#include "testlib.h"
#include "tsk_core.h"

#include <unistd.h>

static void
test_strerror(void)
{
    int j;
    const char *msg;
    int max_error_code = 1024; /* totally arbitrary */

    for (j = 0; j < max_error_code; j++) {
        msg = tsk_strerror(-j);
        CU_ASSERT_FATAL(msg != NULL);
        CU_ASSERT(strlen(msg) > 0);
    }
}

static void
test_strerror_kastore(void)
{
    int kastore_errors[] = {KAS_ERR_NO_MEMORY, KAS_ERR_IO, KAS_ERR_KEY_NOT_FOUND};
    size_t j;
    int err;

    for (j = 0; j < sizeof(kastore_errors) / sizeof(*kastore_errors); j++) {
        err = tsk_set_kas_error(kastore_errors[j]);
        CU_ASSERT_TRUE(tsk_is_kas_error(err));
        CU_ASSERT_STRING_EQUAL(tsk_strerror(err), kas_strerror(kastore_errors[j]));
    }
}

static void
test_generate_uuid(void)
{
    size_t uuid_size = 36;
    char uuid[uuid_size + 1];
    char other_uuid[uuid_size + 1];
    int ret;

    ret = tsk_generate_uuid(uuid, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(strlen(uuid), uuid_size);
    CU_ASSERT_EQUAL(uuid[8], '-');
    CU_ASSERT_EQUAL(uuid[13], '-');
    CU_ASSERT_EQUAL(uuid[18], '-');
    CU_ASSERT_EQUAL(uuid[23], '-');

    ret = tsk_generate_uuid(other_uuid, 0);
    CU_ASSERT_EQUAL_FATAL(ret, 0);
    CU_ASSERT_EQUAL_FATAL(strlen(other_uuid), uuid_size);
    CU_ASSERT_STRING_NOT_EQUAL(uuid, other_uuid);
}

int
main(int argc, char **argv)
{
    CU_TestInfo tests[] = {
        {"test_strerror", test_strerror},
        {"test_strerror_kastore", test_strerror_kastore},
        {"test_generate_uuid", test_generate_uuid},
        {NULL},
    };

    return test_main(tests, argc, argv);
}
