
#ifndef SERIAL_UTILS_H
#define SERIAL_UTILS_H

int RandomInt(int upper);

int* IntSet(int n, int val, int* target);

void IntCopy(int n, int* source, int*dest);

void FreePtrList(void *ptr);
#endif //SERIAL_UTILS_H
