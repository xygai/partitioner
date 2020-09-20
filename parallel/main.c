#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"
#include "graph.h"
#include "multilevel.h"


int main(int argc, char* argv[]) {

    int nparts, edgecut, redgecut;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s <InFile> <OutFile> <nparts>\n", argv[0]);
        exit(EXIT_FAILURE);
    }


    MPI_Comm comm;
    int nprocs, rank;
    double start, end, rstart, rend;

    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    //printf("this is process %d\n", rank);

    /* start timer */
    start = MPI_Wtime();

    nparts = atoi(argv[3]);

    graph_t *graph = ReadGraph(argv[1]);

    //if(rank == 0){
    //    printf("vertex count:    %d\nedge count:      %d\n", TotalVertexWeight(graph), TotalEdgeWeight(graph));
    //}

    /* print the original graph */
    //PrintGraph(graph);
    //printf("total vertex weight = %d, total edge weight = %d\n", TotalVertexWeight(graph), TotalEdgeWeight(graph));

    int nvtxs = graph->nvtxs;
    int *part = malloc(sizeof(int) * nvtxs);
    IntSet(nvtxs, 0, part);


    /* partition */
    InitPartParallel(graph, nparts, part, &edgecut, comm);

    /* end timer */
    end = MPI_Wtime();

    if(rank == 0){
        printf("nprocs = %d, nparts = %d, ngroups = %d\n", nprocs, nparts, (nprocs/(nparts/2)) < 1 ? 1: nprocs/(nparts/2));
        printf("initial edgecut: %d\n", edgecut);
        WritePartition(argv[2], nvtxs, part);
    }

    /* test for refinement */
    graph = ReadGraph(argv[1]);

    /* start refinement timer */
    rstart = MPI_Wtime();

    Refine(graph, part, nparts, comm);

    /* end refinement timer */
    rend = MPI_Wtime();

    if(rank == 0){
        redgecut = ComputeEdgeCut(graph);
        printf("updated edgecut: %d\n", redgecut);
        printf("diff in edgecut: %d\n", edgecut - redgecut);
        printf("partition  time: %f\n", end - start);
        printf("refinement time: %f\n", rend - rstart);
        printf("execution  time: %f\n", (rend - rstart) + (end - start));
        WritePartition("refine.piece.part", nvtxs, part);
    }

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
