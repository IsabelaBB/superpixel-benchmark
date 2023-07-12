
#ifndef _GFT_GRAPH_H_
#define _GFT_GRAPH_H_

#include "gft_common.h"
#include "gft_stack.h"

namespace gft{
  namespace Graph{

    typedef struct _Node {
      int outdegree;
      int arraysize;
      /*Adjacency lists (lists of nodes that are adjacent to a given node)*/
      int   *adjList;
      float *Warcs;
    } GraphNode;

    typedef struct _Graph {
      int nnodes;
      GraphNode *nodes;
      float     *Wnodes;
    } Graph;

    /* A utility function that creates a graph of 'nnodes' vertices */
    Graph *Create(int nnodes, int outdegree, float *Wnodes);
    Graph *Clone(Graph *graph);
    Graph *Transpose(Graph *graph);
    
    void   Destroy(Graph **graph);

    float GetArcWeight(Graph *graph, int src, int dest);
    
    /* Adds an edge to an undirected graph */
    void AddEdge(Graph *graph, int src, int dest, float w);
    
    /* Adds an arc to a directed graph */
    void AddDirectedEdge(Graph *graph, int src, int dest, float w);

    void UpdateEdge(Graph *graph, int src, int dest, float w);
    void UpdateDirectedEdge(Graph *graph, int src, int dest, float w);

    void UpdateDirectedEdgeIfHigher(Graph *graph, int src, int dest, float w);

    float GetMaximumArc(Graph *graph);
    
    /*Tarjan's algorithm is a procedure for finding strongly connected components of a directed graph.*/
    int *Tarjan(Graph *graph);

    int *Tarjan(Graph *graph, int p);

    int *Tarjan(Graph *graph, int *V, int n);
    
  } //end Graph namespace
} //end gft namespace

    
#endif

