#ifndef __SAMPLE_H
#define __SAMPLE_H

// Sample class    
typedef struct sample_s{
    int kt;
    double dJ;
} Sample;

extern int dJCompare(const void *s1, const void *s2);

#endif	// __SAMPLE_H
