
#include "gft_ift.h"

namespace gft{
  namespace ift{

    /*
    void DrawSeeds(Image32::Image32 *img, int *S, int val){
      int nS,i;
      if(S == NULL) 
	return;
      nS = S[0];
      for(i=1; i<=nS; i++){
	img->data[S[i]] = val;
      }
    }
    */

    
    /*Weighted Distance Transform.*/
    Image32::Image32 *pred_IFTSUM(SparseGraph::SparseGraph *sg,
				  int *S,
				  Image32::Image32 *label,
				  float power, int type){
      Heap::Heap *Q=NULL;
      Image32::Image32 *pred;
      int i,p,q,n;
      float edge,tmp;
      float *cost=NULL;
      int u_x,u_y,v_x,v_y;
      AdjRel::AdjRel *A;
      float *Dpq;
      
      n    = sg->ncols*sg->nrows;
      pred = Image32::Create(sg->ncols, sg->nrows);
      cost = gft::AllocFloatArray(n);
      Q = Heap::Create(n, cost);
      A = sg->A;

      //--------------------
      Dpq = (float *)malloc(A->n*sizeof(float));
      for(i=1; i<A->n; i++){
	Dpq[i] = sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i]);
      }
      //--------------------
      
      Image32::Set(pred, NIL);
      for(p=0; p<n; p++){
	if(label != NULL && label->data[p] != NIL) cost[p] = 0.0;
	else cost[p] = FLT_MAX;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++){
	  cost[S[i]] = 0.0;
	  Heap::Insert_MinPolicy(Q, S[i]);
	}
      }
      
      while(!Heap::IsEmpty(Q)){
	Heap::Remove_MinPolicy(Q, &p);
	u_x = p%sg->ncols; 
	u_y = p/sg->ncols; 
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x >= 0 && v_x < sg->ncols &&
	     v_y >= 0 && v_y < sg->nrows){
	    q = v_x + sg->ncols*v_y;
	    if(Q->color[q] != BLACK){
	      
	      edge = (sg->n_link[p])[i];
	      if(type == 1)
		tmp  = cost[p] + powf(edge, power);
	      else
		tmp  = cost[p] + powf(MAX(edge,1.0), power) - 1.0 + Dpq[i];
	      
	      if(tmp < cost[q]){
		Heap::Update_MinPolicy(Q, q, tmp);
		if(label != NULL)
		  label->data[q] = label->data[p];
		pred->data[q] = p;
	      }
	    }
	  }
	}
      }
      free(Dpq);
      gft::FreeFloatArray(&cost);
      Heap::Destroy(&Q);
      return pred;
    }

    
    void method_IFTW_FIFO_InnerCut(SparseGraph::SparseGraph *sg,
				   int *S,
				   Image32::Image32 *label){
      PQueue32::PQueue32 *Q=NULL;
      int i,j,p,q,n;
      int edge,tmp;
      Image32::Image32 *value;
      int u_x,u_y,v_x,v_y;
      AdjRel::AdjRel *A;
      
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->n;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);
      A = sg->A;
      
      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }

      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }
      
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(Image32::IsValidPixel(label,v_x,v_y)){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p] != 0){
		j = SparseGraph::get_edge_index(q, p, sg);
		edge = (sg->n_link[q])[j] + 1;
	      }
	      else
		edge = (sg->n_link[p])[i] + 1;
	      
	      tmp  = edge;
	      
	      if(tmp < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = tmp;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
    }
    
    
    void method_IFTW_FIFO_OuterCut(SparseGraph::SparseGraph *sg,
				   int *S,
				   Image32::Image32 *label){
      PQueue32::PQueue32 *Q=NULL;
      int i,j,p,q,n;
      int edge,tmp;
      Image32::Image32 *value;
      int u_x,u_y,v_x,v_y;
      AdjRel::AdjRel *A;
      
      value = Image32::Create(sg->ncols,
			      sg->nrows);
      n = label->ncols*label->nrows;
      Q = PQueue32::Create(sg->Wmax+2,n,value->data);
      A = sg->A;

      for(p=0; p<n; p++){
	if(label->data[p]==NIL) value->data[p] = INT_MAX;
	else                    value->data[p] = 0;
      }
      
      if(S != NULL){
	for(i=1; i<=S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p=0; p<n; p++)
	  if(label->data[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);	    
      }

      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::FastRemoveMinFIFO(Q);
	u_x = p%label->ncols; //PixelX(label, p);
	u_y = p/label->ncols; //PixelY(label, p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(Image32::IsValidPixel(label,v_x,v_y)){
	    q = v_x + label->ncols*v_y;
	    if(Q->L.elem[q].color != BLACK){
	      
	      if(label->data[p]==0){
		j = SparseGraph::get_edge_index(q, p, sg);
		edge = (sg->n_link[q])[j] + 1;
	      }
	      else
		edge = (sg->n_link[p])[i] + 1;
	      
	      tmp  = edge;
	      
	      if(tmp < value->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::FastRemoveElem(Q, q);
		value->data[q] = tmp;
		label->data[q] = label->data[p];
		PQueue32::FastInsertElem(Q, q);
	      }
	    }
	  }
	}
      }
      Image32::Destroy(&value);
      PQueue32::Destroy(&Q);
    }
    

    //----------------------------------------

    /*
    Image32::Image32 *EDT(Image32::Image32 *bin, float r){
      Image32::Image32 *Dx=NULL, *Dy=NULL, *cost;
      PQueue32::PQueue32 *Q=NULL;
      int i,p,q,u_x,u_y,v_x,v_y,dx,dy,tmp;
      AdjRel::AdjRel *A;
      
      A = AdjRel::Circular(r);
      Dx = Image32::Create(bin->ncols, bin->nrows);
      Dy = Image32::Create(bin->ncols, bin->nrows);
      cost = Image32::Create(bin->ncols, bin->nrows);
      Q = PQueue32::Create(bin->ncols+bin->nrows, bin->n, cost->data);
      
      for(p = 0; p < bin->n; p++){
	if(bin->data[p] > 0){
	  cost->data[p] = 0;
	  PQueue32::InsertElem(&Q, p);
	}
	else
	  cost->data[p] = INT_MAX;
      }
      while(!PQueue32::IsEmpty(Q)) {
	p = PQueue32::RemoveMinFIFO(Q);
	u_x = p%bin->ncols;
	u_y = p/bin->ncols;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(Image32::IsValidPixel(bin,v_x,v_y)){
	    q = v_x + bin->ncols*v_y;
	    if(cost->data[p] < cost->data[q]){
	      dx  = Dx->data[p] + abs(v_x-u_x);
	      dy  = Dy->data[p] + abs(v_y-u_y);
	      tmp = dx*dx + dy*dy;
	      if (tmp < cost->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  PQueue32::RemoveElem(Q, q);
		cost->data[q] = tmp;
		Dx->data[q] = dx;
		Dy->data[q] = dy;
		PQueue32::InsertElem(&Q, q);
	      }
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
      AdjRel::Destroy(&A);
      Image32::Destroy(&Dx);
      Image32::Destroy(&Dy);
      return cost;
    }
    */


  } /*end ift namespace*/
} /*end gft namespace*/


