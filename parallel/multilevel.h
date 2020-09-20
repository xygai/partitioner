

#ifndef SERIAL_MULTILEVEL_H
#define SERIAL_MULTILEVEL_H

#include "graph.h"
#include "coarsen.h"
#include "initpart.h"
#include "refine.h"
#include "utils.h"
#include "mpi.h"

int MlevelBisect(graph_t *graph);
int MlevelRecursiveBisect(graph_t *graph, int nparts, int *part, int fpart);
void PartGraphRecursive(int nvtxs, int *xadj, int *adjncy, int *vwgt, int *adjwgt, int nparts, int *part);

void SplitGraph(graph_t *graph, graph_t **r_lgraph, graph_t **r_rgraph);
void InitPartParallel(graph_t * graph, int nparts, int *where1, int *edgecut, MPI_Comm comm);
void KeepPart(graph_t *graph, int *part, int mypart);
int CountPart(int partID, int nvtxs, int *part);
#endif //SERIAL_MULTILEVEL_H
