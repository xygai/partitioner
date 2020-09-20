
#include "refine.h"
#include <stdio.h>
#include <stdlib.h>

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


void Refine(graph_t *graph, int * part, int nparts){


    /* allocate memory for refinement */
    AllocateBisectMemory(graph);
    AllocateRefineMemory(graph);

    /* add the parititon info */
    IntCopy(graph->nvtxs, part, graph->where);

    /* each processor compute boundary vertices of its own part of the graph */
    ComputeRefineParams(graph);

    int i,j,k, count = 0;
    for(i = 0 ; i < graph->nvtxs; i++){
        if(graph->ed[i]>graph->id[i])
            count++;
    }
    //printf("count = %d, ec = %d\n",count, graph->mincut);

    int *totalgain = malloc(sizeof(int) * nparts * nparts);      // stores local total gain by perform the proposed move from i to j, assume maximum 64 partitions
    IntSet(nparts * nparts, 0, totalgain);

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


    /* Free refine memory */
    free(totalgain);
    //Bisect memory freed by FreeGraph?
    FreeRefineMemory(graph);
}
