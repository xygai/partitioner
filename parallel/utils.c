#include "utils.h"
#include <stdlib.h>

// random in within range [0, upper - 1]
int RandomInt(int upper)
{
    int num = rand() % upper;
    return num;
}

// set n elements of target to val
int* IntSet(int n, int val, int* target)
{
    int i;
    for (i = 0; i<n; i++)
    {
        target[i] = val;
    }
    return target;
}

// copy elements of source to dest
void IntCopy(int n, int* source, int*dest)
{
    int i;
    for (i = 0; i < n; i++)
    {
        dest[i] = source[i];
    }
}
