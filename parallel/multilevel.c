
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
    if(nparts > 2)  // added for parallelization
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

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/* partition the graph with adjacency information, this is the entrance of the partitioning scheme */
void PartGraphRecursive(int nvtxs, int *xadj, int *adjncy, int *vwgt, int *adjwgt, int nparts, int *edgecut, int *part)
{

    graph_t *graph;

    /* set up the graph */
    graph = SetupGraph(nvtxs, xadj, adjncy, vwgt, adjwgt);

    /* start the partitioning */
    *edgecut = MlevelRecursiveBisect(graph, nparts, part, 0);

}



void InitPartParallel(graph_t * graph, int nparts, int *where1,  MPI_Comm comm)
{
    int nprocs, rank, ngroups;
    MPI_Comm newcomm;
    int newnprocs, newrank;
    int fproc;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    ngroups = MIN(nprocs, 1);


    // make a copy of the graph
    int *xadj, *adjncy, *vwgt, *adjwgt;
    int edgecut = 0;
    int gnvtxs = graph->nvtxs;
    int gnedges = graph->nedges;

    xadj = malloc(sizeof(int) * (gnvtxs + 1));
    vwgt = malloc(sizeof(int) * (gnvtxs + 1));
    adjncy = malloc(sizeof(int) * gnedges);
    adjwgt = malloc(sizeof(int) * gnedges);

    int i ;
    for(i = 0; i < graph->nvtxs + 1; i++){
        xadj[i] = graph->xadj[i];
        vwgt[i] = graph->vwgt[i];
    }
    for(i = 0; i < graph->nedges; i++){
        adjncy[i] = graph->adjncy[i];
        adjwgt[i] = graph->adjwgt[i];
    }

    // create processor groups
    MPI_Comm_split(comm, rank%ngroups, 0, &newcomm);
    MPI_Comm_rank(newcomm, &newrank);
    MPI_Comm_size(newcomm, &newnprocs);

    // rb
    int lnparts = nparts;
    int fpart = fproc = 0;
    int lnprocs = newnprocs;
    int *part, *where0;
    part = malloc(sizeof(int) * gnvtxs);

    IntSet(gnvtxs, 0, part);

    while(lnprocs > 1 && lnparts > 1){

        //IntSet(gnvtxs, 0, part);
        //MlevelRecursiveBisect(graph, 2, part, 0);   //take parameters as inputs instead of graphs
        PartGraphRecursive(graph->nvtxs, graph->xadj, graph->adjncy, graph->vwgt, graph->adjwgt, 2, &edgecut, part);

        //int count0 = CountPart(0, graph->nvtxs, part);
        //int count1 = CountPart(1, graph->nvtxs, part);
        //printf("[ID:%d] [newID:%d] lnprocs = %d, lnparts = %d, nvtxs = %d, count0 = %d, count1 = %d\n", rank, newrank, lnprocs, lnparts, graph->nvtxs, count0, count1);

        /* pick one part from the bisection */
        if(newrank < fproc + lnprocs/2){
            KeepPart(graph, part, 0);   // this shrinks the graph on the fly
            lnprocs = lnprocs/2;
            lnparts = lnparts/2;
        } else {
            KeepPart(graph, part, 1);
            fpart = fpart + lnparts/2;
            fproc = fproc + lnprocs/2;
            lnprocs = lnprocs - lnprocs/2;
            lnparts = lnparts - lnparts/2;
        }

        //printf("myid = %d, mynewid = %d, lnprocs = %d, lnparts = %d, nvtxs = %d\n", rank, newrank, lnprocs, lnparts,  graph->nvtxs);
    }

    /* check processor groups */
    //printf("myid = %d, mynewid = %d, nvtxs = %d\n", rank, newrank, graph->nvtxs);

    where0 = malloc(sizeof(int) * gnvtxs);
    IntSet(gnvtxs, 0, where0);

    /* assume nprocs == nparts, which indicates lnpart == 1 */
    if(lnparts == 1){
        if(newrank == fproc){
            for(i=0; i< graph->nvtxs; i++)
                where0[graph->label[i]] = fpart;    // graph already modified
            //printf("[ID: %d] [newID: %d] end here!\n", rank, newrank);

        }
    } else {
        //printf("[ID: %d] [newID: %d] bisect again!\n", rank, newrank);
        //IntSet(gnvtxs, -1, part);
        //MlevelRecursiveBisect(graph, lnparts, part, 0);
        PartGraphRecursive(graph->nvtxs, graph->xadj, graph->adjncy, graph->vwgt, graph->adjwgt, lnparts, &edgecut, part);
        for(i=0; i< graph->nvtxs; i++)
            where0[graph->label[i]] = fpart + part[i];    // graph already modified

    }

        //int count0 = CountPart(0, graph->nvtxs, where0);
        //int count1 = CountPart(1, graph->nvtxs, where0);
        //int count2 = CountPart(2, graph->nvtxs, where0);
        //int count3 = CountPart(3, graph->nvtxs, where0);
        //printf("[ID:%d] [newID:%d] lnprocs = %d, lnparts = %d, nvtxs = %d, count0 = %d, count1 = %d, count2 = %d, count3 = %d,fpart = %d, fpe = %d\n", rank, newrank, lnprocs, lnparts, graph->nvtxs, count0, count1, count2, count3, fpart, fproc);

    MPI_Allreduce(where0, where1, gnvtxs, MPI_INT, MPI_SUM, newcomm);

    int totalcut = 0;
    MPI_Allreduce(&edgecut, &totalcut, 1, MPI_INT, MPI_SUM, newcomm);
    if(rank == 0){
        printf("total cut = %d\n", totalcut);
    }

    FreeGraph(&graph);
    free(where0);
    free(part);
    MPI_Comm_free(&newcomm);
    MPI_Barrier(comm);

}


void KeepPart(graph_t *graph, int *part, int mypart){
    int h, i, j, k;
    int nvtxs, mynvtxs, mynedges;
    int *xadj, *vwgt, *adjncy, *adjwgt, *label;
    int *rename;

    nvtxs  = graph->nvtxs;
    xadj   = graph->xadj;
    vwgt   = graph->vwgt;
    adjncy = graph->adjncy;
    adjwgt = graph->adjwgt;
    label  = graph->label;

    rename = malloc(sizeof(int)* nvtxs);
    //IntSet(nvtxs, -1, rename);

    for (mynvtxs=0, i=0; i<nvtxs; i++) {
        if (part[i] == mypart)
            rename[i] = mynvtxs++;
    }

    for (mynvtxs=0, mynedges=0, j=xadj[0], i=0; i<nvtxs; i++) {
        if (part[i] == mypart) {
            for (; j<xadj[i+1]; j++) {
                k = adjncy[j];
                if (part[k] == mypart) {
                    adjncy[mynedges] = rename[k];
                    adjwgt[mynedges++] = adjwgt[j];
                }
            }
            j = xadj[i+1];  /* Save xadj[i+1] for later use */

            vwgt[mynvtxs+h] = vwgt[i+h];

            label[mynvtxs] = label[i];
            xadj[++mynvtxs] = mynedges;
        }
        else {
            j = xadj[i+1];  /* Save xadj[i+1] for later use */
        }
    }

    graph->nvtxs  = mynvtxs;
    graph->nedges = mynedges;
    free(rename);
}

int CountPart(int partID, int nvtxs, int *part){
    int count = 0;
    int i;

    for(i = 0; i < nvtxs;i++){
        if(part[i]== partID)
            count++;
    }
    return count;
}
