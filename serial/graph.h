
#ifndef SERIAL_GRAPH_H
#define SERIAL_GRAPH_H

typedef struct graph_t {
    int nvtxs, nedges;	    /* The number of vertices and edges in the graph */

    int *xadj;		        /* CSR format Pointers to the locally stored vertices */
    int *vwgt;		        /* Vertex weights */

    int *adjncy;            /* CSR format adjacency lists of nvtxs */
    int *adjwgt;            /* CSR format wights of edges. adjwgt[i] is weight of adjncy[i] */

    int *label;             /* vertex label of the top level graph, which remains unchanged in split graphs */

    /* Partition parameters */
    int *where;             /* to which partition a vertex belongs */

    /* Refine parameters */
    int nbnd;               /* number of vertices on the boundary of a partition */
    int *bndptr;            /* bndptr[10] = 5, the 10th vertex of the original input graph is the 5th boundary vertex */
    int *bndind;            /* bndind[5] = 10; the 5th boundary vertex is the 10th vertex of the original input graph */

} graph_t;

/* setup graph */
graph_t *CreateGraph();
void InitGraph(graph_t *graph);
graph_t *SetupSplitGraph(graph_t *graph, int snvtxs, int snedges);
//graph_t *SetupGraph(int nvtxs, int *xadj, int *adjncy, int *vwgt, int *adjwgt, int *where);


/* Free allocated memory */
void FreeGraph(graph_t **r_graph);


/* Read a METIS format graph */
graph_t *ReadGraph(char *filename);
void PrintGraph(graph_t *graph);


/* write partition to file */
void WritePartition(char *filename, int nvtxs, int *part);


int TotalVertexWeight(graph_t *graph);
int TotalEdgeWeight(graph_t *graph);


#endif //SERIAL_GRAPH_H
