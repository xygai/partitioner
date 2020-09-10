#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"
#include "graph.h"
#include "multilevel.h"


int main(int argc, char* argv[]) {

    int nparts;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s <InFile> <OutFile> <nparts>\n", argv[0]);
        exit(EXIT_FAILURE);
    }


    MPI_Comm comm;
    int nprocs, rank;
    double start, end;

    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    //printf("this is process %d\n", rank);

    nparts = atoi(argv[3]);

    graph_t *graph = ReadGraph(argv[1]);

    /* print the original graph */
    //PrintGraph(graph);
    //printf("total vertex weight = %d, total edge weight = %d\n", TotalVertexWeight(graph), TotalEdgeWeight(graph));

    int nvtxs = graph->nvtxs;
    int *part = malloc(sizeof(int) * nvtxs);
    IntSet(nvtxs, 0, part);

    /* start timer */
    start = MPI_Wtime();

    /* partition */
    InitPartParallel(graph, nparts, part, comm);

    /* end timer */
    end = MPI_Wtime();

    if(rank == 0){
        printf("run time = %f\n", end - start);
        WritePartition(argv[2], nvtxs, part);
    }

    /* test for refinement */
    graph = ReadGraph(argv[1]);

    Refine(graph, part, nparts, comm);

    free(part);
    FreeGraph(&graph);

    MPI_Finalize();

    return 0;
}


/* test the splitgraph function in multilevel.c */
void TestSplit(graph_t *graph){
    graph_t *lgraph = NULL, *rgraph = NULL;
    SplitGraph(graph, &lgraph, &rgraph);
    printf("l_nvtxs = %d l_nedges = %d\n", lgraph->nvtxs, lgraph->nedges);
    printf("r_nvtxs = %d r_nedges = %d\n", rgraph->nvtxs, rgraph->nedges);

}
