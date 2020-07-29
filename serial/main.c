#include <stdio.h>
#include <stdlib.h>

#include "graph.h"
#include "multilevel.h"


int main(int argc, char* argv[]) {

    int nparts;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s <InFile> <OutFile> <nparts>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    nparts = atoi(argv[3]);

    graph_t *graph = ReadGraph(argv[1]);

    /* print the original graph */
    //PrintGraph(graph);
    printf("total vertex weight = %d, total edge weight = %d\n", TotalVertexWeight(graph), TotalEdgeWeight(graph));

    int nvtxs = graph->nvtxs;
    int *part = malloc(sizeof(int) * nvtxs);

    /* partition */
    //BisectGGP(graph);
    int edgecut = MlevelRecursiveBisect(graph, nparts, part,0);

    //TestSplit(graph);

    WritePartition(argv[2], nvtxs, part);
    printf("total edge cut = %d\n", edgecut);

    free(part);
    //FreeGraph(&graph);

    return 0;
}


/* test the splitgraph function in multilevel.c */
void TestSplit(graph_t *graph){
    graph_t *lgraph = NULL, *rgraph = NULL;
    SplitGraph(graph, &lgraph, &rgraph);
    printf("l_nvtxs = %d l_nedges = %d\n", lgraph->nvtxs, lgraph->nedges);
    printf("r_nvtxs = %d r_nedges = %d\n", rgraph->nvtxs, rgraph->nedges);

}