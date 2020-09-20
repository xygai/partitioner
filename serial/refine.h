
#ifndef SERIAL_REFINE_H
#define SERIAL_REFINE_H

#include "graph.h"
#include "utils.h"
#include "initpart.h"

void AllocateRefineMemory(graph_t *graph);
void FreeRefineMemory(graph_t *graph);

void ComputeRefineParams(graph_t *graph);

void Refine(graph_t *graph, int * where, int nparts);


#endif //SERIAL_REFINE_H
