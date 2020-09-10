
#ifndef SERIAL_REFINE_H
#define SERIAL_REFINE_H

#include "graph.h"
#include "utils.h"
#include "initpart.h"

#include "mpi.h"

void AllocateRefineMemory(graph_t *graph);
void FreeRefineMemory(graph_t *graph);
int decomp1d(int n, int p, int rank, int *s, int *e);

void ComputeRefineParams(graph_t *graph);
void Refine(graph_t *graph, int * where, int nparts,  MPI_Comm comm);
#endif //SERIAL_REFINE_H
