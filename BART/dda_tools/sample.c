#include "sample.h"

int dJCompare(const void *s1, const void *s2){
    return ((Sample *)s1)->dJ > ((Sample *)s2)->dJ ? 1 : -1;
}
