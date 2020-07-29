

#ifndef SERIAL_MULTILEVEL_H
#define SERIAL_MULTILEVEL_H

#include "graph.h"
#include "coarsen.h"
#include "initpart.h"
#include "refine.h"
#include "utils.h"

int MlevelBisect(graph_t *graph);
int MlevelRecursiveBisect(graph_t *graph, int nparts, int *part, int fpart);

void SplitGraph(graph_t *graph, graph_t **r_lgraph, graph_t **r_rgraph);
#endif //SERIAL_MULTILEVEL_H
