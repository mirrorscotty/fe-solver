#ifndef OUTPUT_H
#define OUTPUT_H

#include "finite-element1d.h"

void CSVOutFixedNode(struct fe1d*, int, char*);
void CSVOutFixedNode2(struct fe1d*, int, char*);
void CSVOutAvg(struct fe1d*, int, char*);
void CSVOutFixedTime(struct fe1d*, int, char*);

#endif

