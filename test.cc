#include <stdio.h>
#include <stdbool.h>

/*
 * Simple testing framework
 */

#define TEST(a, str) {                          \
    if (a)                                      \
        fprintf(stdout, "\x1b[32m(P) %s\x1b[0m\n", str);                \
    else                                        \
        fprintf(stderr, "\x1b[31m(F) %s@%d\x1b[0m\n", str, __LINE__);   \
}                                               \

int
main(void)
{
    TEST(false, "failing test");
    TEST(true, "passing test");
}
