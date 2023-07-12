
#include "gft_layeredgraph.h"

namespace gft{
  namespace LayeredGraph{

    LayeredGraph *Create(int nlayers, int ncols, int nrows){
      LayeredGraph *lg;
      int n;
      n = ncols*nrows;
      lg = (LayeredGraph *)calloc(1, sizeof(LayeredGraph));
      if(lg == NULL)
	gft::Error((char *)MSG1,(char *)"LayeredGraph::Create");
      lg->nlayers = nlayers;
      lg->ncols = ncols;
      lg->nrows = nrows;
      lg->graph = gft::Graph::Create(nlayers*n, 8+8*nlayers, NULL);
      return lg;
    }

    
    void Destroy(LayeredGraph **lg){
      LayeredGraph *tmp;
      if(lg == NULL) return;
      tmp = *lg;
      if(tmp != NULL){
	gft::Graph::Destroy(&tmp->graph);
	free(tmp);
      }
      *lg = NULL;
    }

    
    void SetArcs(LayeredGraph *lg, gft::SparseGraph::SparseGraph *sg, int l){
      int p, q, i, n, src, dest;
      int p_x, p_y, q_x, q_y;
      float wpq;
      n = sg->ncols*sg->nrows;
      for(p = 0; p < n; p++){
	p_x = p % sg->ncols;
	p_y = p / sg->ncols;
	for(i = 1; i < sg->A->n; i++){
	  q_x = p_x + sg->A->dx[i];
	  q_y = p_y + sg->A->dy[i];

          if ((q_x >= 0)&&(q_x < lg->ncols)&&
	  (q_y >= 0)&&(q_y < lg->nrows)){

	  	q = q_x + q_y*sg->ncols;

	  	wpq = (float)((sg->n_link[p])[i]);
          	//printf("Wpq %f\n",wpq);

	  	src = p + l*n;
	  	dest = q + l*n;
	  	gft::Graph::AddDirectedEdge(lg->graph, src, dest, wpq);
	  }
	  /*
	  else{ 
		q = q_x + q_y*sg->ncols;
		src = p + l*n;
	  	dest = q + l*n;
		gft::Graph::AddDirectedEdge(lg->graph, src, dest, NIL);
	  }
	  */
	}
      }
    }

    
    void SetArcs(LayeredGraph *lg,
		 int l_orig, int l_dest,
		 float w, float r){
      gft::AdjRel::AdjRel *A;
      int n, p, q, i, src, dest;
      int p_x, p_y, q_x, q_y;
      n = lg->ncols*lg->nrows;
      A = gft::AdjRel::Circular(r);
       
      //printf("Entro set Arcs 2 \n");

      for(p = 0; p < n; p++){
	p_x = p % lg->ncols;
	p_y = p / lg->ncols;

	for(i = 0; i < A->n; i++){
	  q_x = p_x + A->dx[i];
	  q_y = p_y + A->dy[i];
          //printf("p_x, p_y = %d, %d \n",p_x,p_y);
          //printf("q_x, q_y = %d, %d \n",q_x,q_y);

	  if ((q_x >= 0)&&(q_x < lg->ncols)&&
	  (q_y >= 0)&&(q_y < lg->nrows)){

	  	q = q_x + q_y*lg->ncols;

	  	src  = p + l_orig*n;
	  	dest = q + l_dest*n;

          	//printf("Wpq_L %f\n",w);
	  	gft::Graph::AddDirectedEdge(lg->graph, src, dest, w);
	  }
	  /*
	  else{ 
		q = q_x + q_y*lg->ncols;

	  	src  = p + l_orig*n;
	  	dest = q + l_dest*n;
		gft::Graph::AddDirectedEdge(lg->graph, src, dest, NIL);
	  }
	  */
	}
      }
      gft::AdjRel::Destroy(&A);
    }

    
    void TransposeLayer(LayeredGraph *lg, int l){
      int p, q, n, i;
      float wpq, wqp;
      n = lg->ncols*lg->nrows;
      for(p = l*n; p < (l+1)*n; p++){
	for(i = 0; i < lg->graph->nodes[p].outdegree; i++){
	  q = lg->graph->nodes[p].adjList[i];
	  if(q >= l*n && q < (l+1)*n){
	    if(q > p){
	      wpq = lg->graph->nodes[p].Warcs[i];
	      wqp = gft::Graph::GetArcWeight(lg->graph, q, p);
	      gft::Graph::UpdateDirectedEdge(lg->graph, q, p, wpq);
	      gft::Graph::UpdateDirectedEdge(lg->graph, p, q, wqp);
	    }
	  }
	}
      }
      
    }
    
    
  } //end LayeredGraph namespace
} //end gft namespace
