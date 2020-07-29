
#include <stdlib.h>
#include <stdio.h>
#include "graph.h"
#include "utils.h"

// create a graph
graph_t *CreateGraph()
{
    graph_t *graph = malloc(sizeof(graph_t)) ;

    InitGraph(graph);

    return graph;
}

// initialize empty graph
void InitGraph(graph_t *graph)
{
    /* graph parameters */
    graph->nvtxs    = -1;
    graph->nedges   = -1;
    graph->xadj     = NULL;
    graph->vwgt     = NULL;
    graph->adjncy   = NULL;
    graph->adjwgt   = NULL;
    graph->label    = NULL;

    /* partition parameters */
    graph->where    = NULL;
    graph->bndptr   = NULL;
    graph->bndind   = NULL;
    graph->nbnd     = -1;
}

// free graph
void FreeGraph(graph_t **r_graph)
{
    graph_t *graph = *r_graph;
    free(graph->xadj);
    free(graph->vwgt);
    free(graph->adjncy);
    free(graph->adjwgt);
    free(graph->label);

    /* free partition info */
    free(graph->where);
    free(graph->bndptr);
    free(graph->bndind);

    free(graph);
    *r_graph = NULL;
}

//read graph from input file
graph_t *ReadGraph(char *filename)
{
    int nfields;
    graph_t *graph;
    char *line      = NULL;
    size_t len      = 0;
    int vcount      = 0;
    int ecount      = 0;
    int edge;
    char *newstr, *curstr;

    graph = CreateGraph();

    /* open the file */
    FILE *input_file = fopen(filename, "r");
    if(input_file == NULL)
    {
        fprintf(stderr, "Error opening the input graph file!\n");
        exit(EXIT_FAILURE);
    }

    /* Read the first line for number of vertices and edges */
    getline (&line, &len, input_file);
    nfields = sscanf(line, "%d %d", &graph->nvtxs, &graph->nedges);

    if (nfields < 2)
        fprintf(stderr, "The input file does not specify the number of vertices or edges.\n");

    graph->nedges *= 2;

    /* allocate memory for xadj and adjncy */
    graph->xadj             = malloc(sizeof(int) * (graph->nvtxs + 1));
    graph->adjncy           = malloc(sizeof(int) * graph->nedges);
    graph->label            = malloc(sizeof(int) * graph->nvtxs);
    graph->xadj[vcount]     = 0;

    int i = 0;
    for(i = 0; i < graph->nvtxs; i++)
    {
        graph->label[i] = i;
    }

    /* Read the following lines  */
    for(vcount = 1; vcount < graph->nvtxs + 1; vcount++)
    {
        getline (&line, &len, input_file);
        curstr = line;
        newstr = NULL;

        while (1) {
            edge = strtol(curstr, &newstr, 10);

            if (newstr == curstr)
                break;  /* End of line */

            curstr = newstr;

            if (edge < 1 || edge > graph->nvtxs)
                fprintf(stderr, "Edge %d for vertex %d is out of bound\n", edge, vcount);

            if (ecount == graph->nedges)
                fprintf(stderr, "There are more edges in the file than the %d specified.\n", graph->nedges/2);

            graph->adjncy[ecount] = edge - 1;    //edge - 1 for C nubmering

            ecount++;
        }
        graph->xadj[vcount] = ecount;
    }

    if(ecount != graph->nedges){
        fprintf(stderr, "less edges found than specified in the first line\n");
        exit(EXIT_FAILURE);
    }

    /* allocate memory for weights */
    graph->vwgt     = malloc(sizeof(int) * graph->nvtxs);
    graph->adjwgt   = malloc(sizeof(int) * graph->nedges);
    IntSet(graph->nvtxs,  1, graph->vwgt);                  // init all vertex weight to 1
    IntSet(graph->nedges, 1, graph->adjwgt);                // init all edge weigt to 1

    fclose(input_file);
    free(line);

    return graph;
}


// print graph in CSR format
void PrintGraph(graph_t *graph){

    printf("xadj: ");
    for(int i = 0; i != graph->nvtxs+1; ++i)
    {
        printf("%d " ,graph->xadj[i]);
    }
    printf("\n");

    if(graph->where != NULL){
        printf("part: ");

        for(int i =0; i != graph->nvtxs; ++i)
        {
            printf("%d ", graph->where[i]);     // label back to fortran numbering
        }
        printf("\n");
    }

    if(graph->label != NULL){
        printf("label: ");

        for(int i =0; i != graph->nvtxs; ++i)
        {
            printf("%d ", graph->label[i] + 1);     // label back to fortran numbering
        }
        printf("\n");
    }

    printf("adjncy: ");
    for(int i =0; i != graph->nedges; ++i)
    {
        printf("%d ", graph->adjncy[i] + 1);        // adjncy back to fortran numbering
    }
    printf("\n");

}

// write partition to file
void WritePartition(char *filename, int nvtxs, int *part)
{
    FILE *out_file = fopen(filename, "w");

    if(out_file == NULL){
        fprintf(stderr, "Error opening the input graph file!\n");
        exit(EXIT_FAILURE);
    }

    int i;
    for(i = 0; i < nvtxs; i++)
    {
        fprintf(out_file, "%d\n", part[i]);
    }

    fclose(out_file);
}


// setup the split graph
graph_t *SetupSplitGraph(graph_t *graph, int snvtxs, int snedges)
{
    graph_t *sgraph;

    sgraph = CreateGraph();

    sgraph->nvtxs       = snvtxs;
    sgraph->nedges      = snedges;

    sgraph->xadj        = malloc(sizeof(int) * (snvtxs+1));
    sgraph->vwgt        = malloc(sizeof(int) * snvtxs);
    sgraph->adjncy      = malloc(sizeof(int) * snedges);
    sgraph->adjwgt      = malloc(sizeof(int) * snedges);
    sgraph->label       = malloc(sizeof(int) * snvtxs);

    return sgraph;
}


int TotalVertexWeight(graph_t *graph)
{
    int i;
    int sum_vwgt = 0;

    if (graph->vwgt == NULL){
        sum_vwgt = graph->nvtxs;
    } else {
        for(i = 0; i < graph->nvtxs; i++)
        {
            sum_vwgt += graph->vwgt[i];
        }
    }
    return sum_vwgt;
}

int TotalEdgeWeight(graph_t *graph)
{
    int i;
    int sum_ewgt = 0;

    if (graph->adjwgt == NULL){
        sum_ewgt = graph->nedges;
    } else {
        for(i = 0; i < graph->nedges; i++)
        {
            sum_ewgt += graph->adjwgt[i];
        }
    }
    return sum_ewgt;
}