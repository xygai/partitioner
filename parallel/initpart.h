
#ifndef SERIAL_INITPART_H
#define SERIAL_INITPART_H

#include "graph.h"


void AllocateBisectMemory(graph_t *graph);
void BisectGGP(graph_t *graph);
void ComputeBoundary(graph_t* graph);
int ComputeEdgeCut(graph_t *graph);

#endif //SERIAL_INITPART_H
