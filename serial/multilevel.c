
#include "multilevel.h"
#include <stdio.h>
#include <stdlib.h>

//#define debug 0

#ifdef debug
int level = 0;
#endif

// multilevel bisection (not recursive)
int MlevelBisect(graph_t *graph)
{
    //TODO: repeat the steps below a few times and find the optimal one

    //TODO: step 1. Coarsen

    //step 2: init bisection
    BisectGGP(graph);
    ComputeBoundary(graph);     //TODO: Greedy GGP will need info of boundary vertices

    //TODO: step 3 Refine;

    int EdgeCut = ComputeEdgeCut(graph);

    return EdgeCut;
}

// do multilevel recursive bisection
// fpart: base partition id
int MlevelRecursiveBisect(graph_t *graph, int nparts, int *part, int fpart){

    graph_t *lgraph = NULL, *rgraph = NULL;
    int *label, *where;
    int i, nvtxs, total_edge_cut;

    nvtxs   = graph->nvtxs;

#ifdef debug
    if(nparts>2)
        printf("partition a graph at level %d\n", level++);
#endif

    /* perform the bisection */
    total_edge_cut = MlevelBisect(graph);

    where   = graph->where;
    label   = graph->label;

    for (i = 0; i < nvtxs; i++)
        part[label[i]] = where[i] + fpart;          // update partition id (for recursive call)

    if(nparts > 2)
        SplitGraph(graph, &lgraph, &rgraph);

#ifdef debug
    if(lgraph != NULL){
        printf("lnvtxs %d rnvtxs %d\n", lgraph->nvtxs, rgraph->nvtxs);

        printf("lgraph:\n");
        PrintGraph(lgraph);

        printf("lgraph:\n");
        PrintGraph(rgraph);
    }
#endif

    /* Free the memory of the top level graph */
    FreeGraph(&graph);

    /* Do the recursive call */
    if (nparts > 2) {
        total_edge_cut += MlevelRecursiveBisect(lgraph, nparts/2, part, fpart);
        total_edge_cut += MlevelRecursiveBisect(rgraph, nparts/2, part, fpart + nparts/2);
    }

    return total_edge_cut;
}

// SplitGraph into two parts: lgraph(0), rgraph(1)
// METIS SplitGraphPart & graph partitioning book p18
void SplitGraph(graph_t *graph, graph_t **r_lgraph, graph_t **r_rgraph)
{
    graph_t *lgraph, *rgraph;
    int *xadj, *adjncy, *adjwgt, *where, *label;
    int *sxadj[2], *svwgt[2], *sadjncy[2], *sadjwgt[2];
    int i, j, k, curr_nedge, istart, iend, mypart, nvtxs, snvtxs[2], snedges[2], *slabel[2];
    int *start_sadjncy, *start_adjwgt;                           // point to the start of sadjncy[0/1] and sadjwgt[0/1]

    nvtxs   = graph->nvtxs;
    xadj    = graph->xadj;
    adjncy  = graph->adjncy;
    adjwgt  = graph->adjwgt;
    where   = graph->where;
    label   = graph->label;

    int * new_label = malloc(sizeof(int) * nvtxs);           // stores new label of a vertex in the sgraph(split graph) to which the vertex belongs

    snvtxs[0]  = snvtxs[1]  = 0;
    snedges[0] = snedges[1] = 0;

    /* split edges of graph */
    for (i=0; i<nvtxs; i++) {
        mypart = where[i];
        new_label[i] = snvtxs[mypart]++;                          // calculate new label of a vertex in sgraph; increase vertex count of sgraph
        snedges[mypart] += xadj[i+1] - xadj[i];                   // increase edge count of sgraph
    }

    /* allocate memory for split graphs */
    lgraph      = SetupSplitGraph(graph, snvtxs[0], snedges[0]);
    sxadj[0]    = lgraph->xadj;
    svwgt[0]    = lgraph->vwgt;
    sadjncy[0]  = lgraph->adjncy;
    sadjwgt[0]  = lgraph->adjwgt;
    slabel[0]   = lgraph->label;

    rgraph      = SetupSplitGraph(graph, snvtxs[1], snedges[1]);
    sxadj[1]    = rgraph->xadj;
    svwgt[1]    = rgraph->vwgt;
    sadjncy[1]  = rgraph->adjncy;
    sadjwgt[1]  = rgraph->adjwgt;
    slabel[1]   = rgraph->label;

    /* split the graph */
    snvtxs[0]   = snvtxs[1]     = 0;
    snedges[0]  = snedges[1]    = 0;
    sxadj[0][0] = sxadj[1][0]   = 0;

    for (i=0; i<nvtxs; i++) {

        mypart  = where[i];
        istart  = xadj[i];
        iend    = xadj[i+1];

        start_sadjncy   = sadjncy[mypart];
        start_adjwgt    = sadjwgt[mypart];

        curr_nedge = snedges[mypart];                           // current number of edges in sgraph before i is distributed

        for (j=istart; j<iend; j++) {
            k = adjncy[j];
            if (where[k] == mypart) {
                start_sadjncy[curr_nedge] = k;
                start_adjwgt[curr_nedge++] = adjwgt[j];
            }
        }
        snedges[mypart] = curr_nedge;                           // update edge count in the sgraph

        /* copy vertex weights */
        slabel[mypart][snvtxs[mypart]]   = label[i];            // Eg., vertex i is the snvtxs[0] th vertex in partition 0
        sxadj[mypart][++snvtxs[mypart]]  = snedges[mypart];     // update xadj of sgraph
    }

    for (mypart=0; mypart<2; mypart++) {
        iend          = snedges[mypart];
        start_sadjncy = sadjncy[mypart];
        for (i=0; i<iend; i++)
            start_sadjncy[i] = new_label[start_sadjncy[i]];     // update adjncy of split graph (from old label in graph to new label in sgraph)
    }

    lgraph->nedges = snedges[0];
    rgraph->nedges = snedges[1];

    *r_lgraph = lgraph;
    *r_rgraph = rgraph;

    free(new_label);
}