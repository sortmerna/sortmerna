/*
 * test_harness.h -- minimal C test macros
 * Pure C89, no dependencies. Single-TU design: include this header in
 * exactly one .c file that contains main(). Other test files should be
 * #include'd into that file (not compiled separately).
 */
#ifndef TEST_HARNESS_H
#define TEST_HARNESS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int _tests_run = 0;
static int _tests_failed = 0;
static int _test_current_failed = 0;

#define TEST(name) static void name(void)

/* Set SMR_TEST_FILTER to a substring to run only the matching tests. */
#define RUN_TEST(name) do { \
    const char *_filter = getenv("SMR_TEST_FILTER"); \
    if (_filter && !strstr(#name, _filter)) break; \
    printf("  %-60s", #name); \
    _test_current_failed = 0; \
    name(); \
    if (!_test_current_failed) { \
        printf(" PASS\n"); \
    } \
    _tests_run++; \
} while(0)

#define ASSERT_TRUE(expr) do { \
    int _val = !!(expr); \
    if (!_val) { \
        printf(" FAIL\n    %s:%d: %s\n", __FILE__, __LINE__, #expr); \
        _tests_failed++; _test_current_failed = 1; return; \
    } \
} while(0)

#define ASSERT_EQ_INT(a, b) do { \
    int _a = (int)(a); int _b = (int)(b); \
    if (_a != _b) { \
        printf(" FAIL\n    %s:%d: %d != %d\n", __FILE__, __LINE__, _a, _b); \
        _tests_failed++; _test_current_failed = 1; return; \
    } \
} while(0)

#define ASSERT_EQ_U64(a, b) do { \
    unsigned long long _a = (unsigned long long)(a); \
    unsigned long long _b = (unsigned long long)(b); \
    if (_a != _b) { \
        printf(" FAIL\n    %s:%d: %llu != %llu\n", __FILE__, __LINE__, _a, _b); \
        _tests_failed++; _test_current_failed = 1; return; \
    } \
} while(0)

#define ASSERT_EQ_SZ(a, b) do { \
    size_t _a = (size_t)(a); size_t _b = (size_t)(b); \
    if (_a != _b) { \
        printf(" FAIL\n    %s:%d: %zu != %zu\n", __FILE__, __LINE__, _a, _b); \
        _tests_failed++; _test_current_failed = 1; return; \
    } \
} while(0)

#define ASSERT_NE(a, b) ASSERT_TRUE((a) != (b))

#define ASSERT_STR_EQ(a, b) do { \
    const char *_a = (a); const char *_b = (b); \
    if (strcmp(_a, _b) != 0) { \
        printf(" FAIL\n    %s:%d: \"%s\" != \"%s\"\n", __FILE__, __LINE__, _a, _b); \
        _tests_failed++; _test_current_failed = 1; return; \
    } \
} while(0)

#define ASSERT_NULL(p) do { \
    const void *_p = (const void *)(p); \
    if (_p != NULL) { \
        printf(" FAIL\n    %s:%d: expected NULL\n", __FILE__, __LINE__); \
        _tests_failed++; _test_current_failed = 1; return; \
    } \
} while(0)

#define ASSERT_NOT_NULL(p) do { \
    const void *_p = (const void *)(p); \
    if (_p == NULL) { \
        printf(" FAIL\n    %s:%d: unexpected NULL\n", __FILE__, __LINE__); \
        _tests_failed++; _test_current_failed = 1; return; \
    } \
} while(0)

/* Bitwise-equality check for doubles. Use when two code paths should produce
 * the same value with no floating-point reassociation difference expected. */
#define ASSERT_DOUBLE_BITEQ(a, b) do { \
    double _a = (double)(a); double _b = (double)(b); \
    if (memcmp(&_a, &_b, sizeof(double)) != 0) { \
        printf(" FAIL\n    %s:%d: %.17g !bit== %.17g\n", __FILE__, __LINE__, _a, _b); \
        _tests_failed++; _test_current_failed = 1; return; \
    } \
} while(0)

#define ASSERT_DOUBLE_NEAR(a, b, eps) do { \
    double _a = (double)(a); double _b = (double)(b); double _e = (double)(eps); \
    double _d = (_a > _b) ? (_a - _b) : (_b - _a); \
    if (_d > _e) { \
        printf(" FAIL\n    %s:%d: |%.6g - %.6g| = %.6g > %.6g\n", __FILE__, __LINE__, _a, _b, _d, _e); \
        _tests_failed++; _test_current_failed = 1; return; \
    } \
} while(0)

#define TEST_MAIN_BEGIN() int main(void) { printf("Running tests:\n");
#define TEST_MAIN_END() \
    printf("\n%d tests, %d failures\n", _tests_run, _tests_failed); \
    return _tests_failed > 0 ? 1 : 0; }

#endif /* TEST_HARNESS_H */
