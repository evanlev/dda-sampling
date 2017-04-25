//#ifndef SAMPLE_H
//#define SAMPLE_H
#pragma once

// Sample class    
typedef struct sample_s{
    int kt;
    double dJ;
} Sample;

int dJCompare(const void *s1, const void *s2);

//#endif

