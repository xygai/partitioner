
#include "initpart.h"
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>


//allocate memory for bisection
void AllocateBisectMemory(graph_t *graph)
{
    graph->where  = malloc(sizeof(int) * graph->nvtxs);
    graph->bndptr = malloc(sizeof(int) * graph->nvtxs);
    graph->bndind = malloc(sizeof(int) * graph->nvtxs);

    IntSet(graph->nvtxs, -1, graph->where);
    IntSet(graph->nvtxs, -1, graph->bndptr);
    IntSet(graph->nvtxs, -1, graph->bndind);
}


// bisect with Graph Growing Algorithm(GGP)
// partition a graph into 2 parts with strictly equal size
void BisectGGP(graph_t *graph){

    int i, j, k;

    AllocateBisectMemory(graph);

    int nvtxs       = graph->nvtxs;
    int *xadj       = graph->xadj;
    int *adjncy     = graph->adjncy;
    int *where      = graph->where;
    IntSet(nvtxs, 1, where);            // assign all vertices to part 1, then we grow part 0

    /* allocate memory for queue and visit status */
    int *queue, *visited;
    queue       = malloc(nvtxs * sizeof(int));
    visited     = malloc(nvtxs * sizeof(int));
    IntSet(nvtxs, 0, visited);          // visited = 1; not yet visited = 0

    int front = 0, back = 1;                // track the front/back of the queue
    queue[front] = RandomInt(nvtxs);;       // randomly pick a seed to grow
    int nleft  = nvtxs - 1;


    /* Breadth first search */
    //TODO: compute edge cut and return, do multiple times and select the one with minimum edge cut.
    do {
        if (front == back){
            if (nleft == 0)
                printf("empty graph\n");
            else
                printf("disconnected graph\n");
            break;
        }

        i = queue[front++];
        where[i] = 0;                                   // move i from part 1 to part 0

        for (j = xadj[i]; j < xadj[i+1]; j++) {         // j is index of neighbour vertices of i
            k = adjncy[j];                              // k is a neighbour to i
            if (visited[k] == 0) {                      // if k is not visited before
                queue[back++] = k;                      // enqueue k
                visited[k] = 1;                         // mark k as visited
                where[k] = 0;                           // move k to part 0 FIXME
                nleft--;                                // update total number of unvisited vertices
            }
        }
    } while(nleft > nvtxs/2);                           //TODO: use weights instead of vertex count

    //TODO: add 2way Refinement

    free(queue);
    free(visited);
}

void ComputeBoundary(graph_t* graph)
{
    int i, j,  istart, iend, sum_vwgt, sum_adjewgt, mypart;
    int *where = graph->where;

    IntSet(graph->nvtxs, -1, graph->bndptr);        // -1: not on boundary

    /* Compute the boundary info */
    int nbnd = 0;

    for (i = 0; i < graph->nvtxs; i++)
    {
        istart = graph->xadj[i];
        iend   = graph->xadj[i+1];

        mypart = where[i];
        sum_vwgt = sum_adjewgt = 0;                     // total weights of a vertex and edges adjacent to the vertex

        for (j=istart; j<iend; j++)
        {
            if (mypart == where[graph->adjncy[j]])      // if the adjacent vertex is in the same partition
                sum_vwgt += graph->adjwgt[j];
            else
                sum_adjewgt += graph->adjwgt[j];
        }

        if (sum_adjewgt > 0 || istart == iend)          // at least one neighbour is in different partition
        {
            graph->bndind[nbnd] = i;                    // vertex i is on boundary
            graph->bndptr[i] = nbnd++;
        }
    }

    graph->nbnd = nbnd;
}

// compute edgecut of a partitioned graph
int ComputeEdgeCut(graph_t *graph)
{
    int i,j,k;
    int sum = 0;

    for(i=0; i< graph->nvtxs; i++)
    {
        for (j = graph->xadj[i]; j < graph->xadj[i+1]; j++)
        {
            k = graph->adjncy[j];                       // k is a neighbour to i
            if (graph->where[k] != graph->where[i])
            {
                sum += graph->adjwgt[j];                // add weight of edge e(i,k) to sum
            }
        }
    }

    return sum/2;
}