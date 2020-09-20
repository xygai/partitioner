#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "graph.h"
#include "multilevel.h"


int main(int argc, char* argv[]) {

    int nparts;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s <InFile> <OutFile> <nparts>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* time the program */
    struct timeval start, end, rstart, rend;
    gettimeofday(&start, NULL);

    nparts = atoi(argv[3]);

    graph_t *graph = ReadGraph(argv[1]);

    /* print the original graph */
    //PrintGraph(graph);
    printf("vertex count:    %d\n", TotalVertexWeight(graph));
    printf("edge count:      %d\n", TotalEdgeWeight(graph));
    printf("nparts:          %d\n", nparts);

    int nvtxs = graph->nvtxs;
    int *part = malloc(sizeof(int) * nvtxs);

    /* partition */
    //BisectGGP(graph);
    int edgecut = MlevelRecursiveBisect(graph, nparts, part,0);

    //TestSplit(graph);

    /* end time */
    gettimeofday(&end, NULL);
    double time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/ 1000000.0;
    printf("partition  time: %f\n", time_spent);

    WritePartition(argv[2], nvtxs, part);


    /* test for refinement */
    graph = ReadGraph(argv[1]);

    gettimeofday(&rstart, NULL);
    Refine(graph, part, nparts);
    gettimeofday(&rend, NULL);
    double rtime_spent = (rend.tv_sec - rstart.tv_sec) + (rend.tv_usec - rstart.tv_usec)/ 1000000.0;

    int redgecut = ComputeEdgeCut(graph);

    printf("refinement time: %f\n", rtime_spent);
    printf("execution  time: %f\n", time_spent+rtime_spent);
    printf("initial edgecut: %d\nupdated edgecut: %d\ndiff in edgecut: %d\n", edgecut, redgecut ,edgecut - redgecut);

    WritePartition("refine.piece.part", nvtxs, part);

    free(part);
    FreeGraph(&graph);

    return 0;
}


/* test the splitgraph function in multilevel.c */
void TestSplit(graph_t *graph){
    graph_t *lgraph = NULL, *rgraph = NULL;
    SplitGraph(graph, &lgraph, &rgraph);
    printf("l_nvtxs = %d l_nedges = %d\n", lgraph->nvtxs, lgraph->nedges);
    printf("r_nvtxs = %d r_nedges = %d\n", rgraph->nvtxs, rgraph->nedges);

}
