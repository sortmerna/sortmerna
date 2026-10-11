/*
 * test_link.c -- verify the bundled libsmr_api.a is self-contained.
 *
 * This links against the merged archive alone, with no help from the rest of
 * the build tree, and so fails to build if a static dependency is missing from
 * the bundle. All of the public entry points live in the single smr_api.cpp
 * translation unit, so calling one of them pulls smr_api.o in and with it the
 * whole dependency closure - running the program is not what does the testing,
 * linking it is.
 */
#include <stdio.h>

#include "smr_api.h"

int main(void)
{
    const char *version = smr_version();

    if (version == NULL || version[0] == '\0') {
        fprintf(stderr, "smr_version() returned no version\n");
        return 1;
    }

    printf("%s\n", version);
    return 0;
}
