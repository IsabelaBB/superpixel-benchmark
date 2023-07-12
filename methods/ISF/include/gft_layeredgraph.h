
#ifndef _GFT_LAYEREDGRAPH_H_
#define _GFT_LAYEREDGRAPH_H_

#include "gft_common.h"
#include "gft_stack.h"
#include "gft_graph.h"
#include "gft_sparsegraph.h"


namespace gft{
  namespace LayeredGraph{

    typedef struct _LayeredGraph {
      int nlayers;
      int ncols;
      int nrows;
      gft::Graph::Graph *graph;
    } LayeredGraph;

    LayeredGraph *Create(int nlayers, int ncols, int nrows);
    void Destroy(LayeredGraph **lg);
    
    void SetArcs(LayeredGraph *lg, gft::SparseGraph::SparseGraph *sg, int l);
    void SetArcs(LayeredGraph *lg,
		 int l_orig, int l_dest,
		 float w, float r);
    
    void TransposeLayer(LayeredGraph *lg, int l);
    

  } //end LayeredGraph namespace
} //end gft namespace

    
#endif


    
