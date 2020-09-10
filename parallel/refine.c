
#include <stdio.h>
#include <stdlib.h>
#include "refine.h"

void AllocateRefineMemory(graph_t *graph){
    int nvtxs;
    nvtxs = graph->nvtxs;

    graph->id     = malloc(nvtxs* sizeof(int));
    graph->ed     = malloc(nvtxs* sizeof(int));
    graph->to     = malloc(nvtxs* sizeof(int));
    IntSet(nvtxs, 0, graph->id);
    IntSet(nvtxs, 0, graph->ed);
    IntSet(nvtxs, -1, graph->to);
}

void FreeRefineMemory(graph_t *graph){
    free(graph->id);
    free(graph->ed);
    free(graph->to);
}


int decomp1d(int n, int p, int rank, int *s, int *e){
    int count = n/p;

    int rem = n%p;

    if (rank < rem){
        *s = rank * (count + 1);
        *e = *s + count;
    } else {
        *s = rem * (count + 1) + (rank - rem) * count;
        *e = *s + count - 1;
    }

}


void ComputeRefineParams(graph_t *graph)
{
    int i, j, nvtxs, nbnd, mincut, istart, iend, tid, ted, me, tto;
    int *xadj, *adjncy, *adjwgt;
    int *where, *bndptr, *bndind, *id, *ed, *to;

    nvtxs  = graph->nvtxs;
    xadj   = graph->xadj;
    adjncy = graph->adjncy;
    adjwgt = graph->adjwgt;

    where  = graph->where;
    id     = graph->id;
    ed     = graph->ed;
    to     = graph->to;

    bndptr = IntSet(nvtxs, -1, graph->bndptr);
    bndind = graph->bndind;


    /* Compute the required info for refinement  */
    for (nbnd=0, mincut=0, i=0; i<nvtxs; i++) {
        istart = xadj[i];
        iend   = xadj[i+1];

        me = where[i];
        tto = -1;
        tid = ted = 0;

        for (j=istart; j<iend; j++) {
            if (me == where[adjncy[j]])
            {
                tid += adjwgt[j];
            }
            else {
                ted += adjwgt[j];
                tto  = where[adjncy[j]];    // partition of the last neighbour
            }
        }
        id[i] = tid;
        ed[i] = ted;
        if(ted > tid)
            to[i] = tto;

        if (ted > 0 || istart == iend) {
            bndind[nbnd] = i;
            bndptr[i] = nbnd++;
            mincut += ted;
        }
    }

    graph->mincut = mincut/2;
    graph->nbnd   = nbnd;

}


void Refine(graph_t *graph, int * part, int nparts, MPI_Comm comm){

    int nprocs, rank;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    AllocateBisectMemory(graph);
    AllocateRefineMemory(graph);

    int vs, ve; // vertex start, vertex end
    decomp1d(graph->nvtxs,  nprocs,  rank,  &vs, & ve);
    //printf("[ID: %d] start %d end %d total %d\n", rank, vs, ve, graph->nvtxs);

    /* add the parititon info */
    IntCopy(graph->nvtxs, part, graph->where);

    /* each processor compute boundary vertices of its own part of the graph */
    ComputeRefineParams(graph);

    printf("no of boundary vtx = %d, ec = %d\n", graph->nbnd, graph->mincut);

    int i,j,k, count = 0;
    for(i = 0 ; i < graph->nvtxs; i++){
        if(graph->ed[i]>graph->id[i])
            count++;
            //printf("id = %d, ed = %d, from = %d, to = %d\n",graph->id[i], graph->ed[i], graph->where[i], graph->to[i]);
    }
    printf("count = %d\n",count);

    int *totalgain = malloc(sizeof(int) * nparts * nparts);      // stores total gain by perform the proposed move from i to j, assume maximum 64 partitions

    IntSet(nparts * nparts, 0, totalgain);                       // initialize

    /* calculate total gain of the proposed moves */
    int me, mypart, topart;
    for (k = 0; k < graph->nbnd; k++){
        me     = graph->bndind[k];
        mypart = graph->where[me];
        topart = graph->to[me];
        if(topart > -1)
            totalgain[mypart * nparts + topart] += (graph->ed[me] - graph->id[me]);
    }

    /* if totalgain[i][j] > totalgain[j][i], all moves from i to j are allowed */
    for (k = 0; k < graph->nbnd; k++){
        me     = graph->bndind[k];
        mypart = graph->where[me];
        topart = graph->to[me];
        if(topart > -1 && totalgain[mypart*nparts + topart] > totalgain[topart*nparts + mypart])
            graph->where[me] = graph->to[me];
    }

    printf("new total edge cut = %d\n", ComputeEdgeCut(graph));

    /* Free refine memory */
    free(totalgain);
    //Bisect memory freed by FreeGraph?
    //FreeRefineMemory(graph);
}
