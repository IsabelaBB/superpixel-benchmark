
#include "gft_sparsegraph.h"

namespace gft{
  namespace SparseGraph{


    SparseGraph *ByAbsDiff(Image32::Image32 *img, float r){
      SparseGraph *sg=NULL;
      AdjRel::AdjRel *A=AdjRel::Circular(r);
      int p,q,i,n;
      int weight;
      int u_x,u_y,v_x,v_y;
      
      n  = img->ncols*img->nrows;
      sg = Create(img->ncols, 
		  img->nrows, A);
      sg->type = DISSIMILARITY;
      sg->Wmax = 0;
      
      for(p=0; p<n; p++){
	u_x = p%img->ncols; //PixelX(img,p);
	u_y = p/img->ncols; //PixelY(img,p);
	
	(sg->n_link[p])[0] = NIL;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(Image32::IsValidPixel(img,v_x,v_y)){
	    q = v_x + img->ncols*v_y;
	    
	    weight = abs(img->data[p]-img->data[q]);
	    (sg->n_link[p])[i] = weight;
	    //printf("weight %d\n",(sg->n_link[p])[i]);
	    if(weight>sg->Wmax){ 
			sg->Wmax = weight; 
			//printf("Wmax %d\n",sg->Wmax);
		}
	  }
	  else{
	    (sg->n_link[p])[i] = NIL;
	  }
	}
      }
      return sg;
    }


    SparseGraph *ByAccAbsDiff(Image32::Image32 *img, float r, float R){
      SparseGraph *sg=NULL;
      Image32::Image32 *W;
      AdjRel::AdjRel *A=AdjRel::Circular(R);
      int p,q,i;
      int weight;
      int u_x,u_y,v_x,v_y;
      
      W = Image32::Create(img->ncols, img->nrows);
      for(p = 0; p < img->n; p++){
	u_x = p%img->ncols; //PixelX(img,p);
	u_y = p/img->ncols; //PixelY(img,p);
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(Image32::IsValidPixel(img,v_x,v_y)){
	    q = v_x + img->ncols*v_y;
	    
	    weight = abs(img->data[p]-img->data[q]);
	    W->data[p] += weight;
	  }
	}
      }
      sg = ByWeightImage(W, r);
      Image32::Destroy(&W);
      AdjRel::Destroy(&A);
      return sg;
    }
    

    SparseGraph *ByAccAbsDiff(CImage::CImage *cimg, float r, float R){
      SparseGraph *sg=NULL;
      Image32::Image32 *W;
      AdjRel::AdjRel *A=AdjRel::Circular(R);
      int p,q,i;
      float weight,sum;
      int dr,dg,db;
      int u_x,u_y,v_x,v_y;
      
      W = Image32::Create(cimg->C[0]->ncols, cimg->C[0]->nrows);
      for(p = 0; p < cimg->C[0]->n; p++){
	u_x = p%cimg->C[0]->ncols; //PixelX(img,p);
	u_y = p/cimg->C[0]->ncols; //PixelY(img,p);

	sum = 0.0;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(Image32::IsValidPixel(cimg->C[0],v_x,v_y)){
	    q = v_x + cimg->C[0]->ncols*v_y;
	    
	    dr = abs(cimg->C[0]->data[p] - cimg->C[0]->data[q]);
	    dg = abs(cimg->C[1]->data[p] - cimg->C[1]->data[q]);
	    db = abs(cimg->C[2]->data[p] - cimg->C[2]->data[q]);
	    weight = sqrtf(dr*dr + dg*dg + db*db);
	    sum += weight;
	  }
	}
	W->data[p] = ROUND(sum);
      }
      sg = ByWeightImage(W, r);
      Image32::Destroy(&W);
      AdjRel::Destroy(&A);
      return sg;
    }
    
    
    SparseGraph *ByWeightImage(Image32::Image32 *W, float r){
      SparseGraph *sg=NULL;
      AdjRel::AdjRel *A=AdjRel::Circular(r);
      int p,q,i,n;
      int weight;
      int u_x,u_y,v_x,v_y;
      
      n  = W->ncols*W->nrows;
      sg = Create(W->ncols, 
		  W->nrows, A);
      sg->type = DISSIMILARITY;
      sg->Wmax = 0;
      
      for(p=0; p<n; p++){
	u_x = p%W->ncols; //PixelX(W,p);
	u_y = p/W->ncols; //PixelY(W,p);
	
	(sg->n_link[p])[0] = NIL;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(Image32::IsValidPixel(W,v_x,v_y)){
	    q = v_x + W->ncols*v_y;
	    
	    weight = W->data[p] + W->data[q];
	    (sg->n_link[p])[i] = weight;
	    
	    if(weight>sg->Wmax) sg->Wmax = weight;
	  }
	  else{
	    (sg->n_link[p])[i] = NIL;
	  }
	}
      }
      return sg;
    }


    SparseGraph   *ByHomogeneityAffinity(Image32::Image32 *img, float r){
      SparseGraph *sg=NULL;
      gft::AdjRel::AdjRel *A = gft::AdjRel::Circular(r);
      int p,q,i,n;
      int weight;
      gft::Pixel u,v;
      
      n        = img->ncols*img->nrows;
      sg       = Create(img->ncols, img->nrows, A);
      sg->type = DISSIMILARITY;
      sg->Wmax = 0;
      
      for(p=0; p<n; p++){
	u.x = p%img->ncols;
	u.y = p/img->ncols;
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  if(gft::Image32::IsValidPixel(img,v.x,v.y)){
	    q = v.x + img->ncols * v.y;
	    
	    weight = abs(img->data[p]-img->data[q]);
	    (sg->n_link[p])[i] = weight;
	    
	    if(weight>sg->Wmax) sg->Wmax = weight;
	  }
	  else{
	    (sg->n_link[p])[i] = NIL;
	  }
	}
      }
      sg->Wmax++;
      ChangeType(sg, CAPACITY);
      return sg;
    }
    
    
    SparseGraph *ByCHomogeneityAffinity(CImage::CImage *cimg, float r){
      SparseGraph *sg=NULL;
      gft::AdjRel::AdjRel *A = gft::AdjRel::Circular(r);
      int p,q,i,n,ncols,nrows;
      int weight,dr,dg,db;
      gft::Pixel u,v;
      
      ncols = cimg->C[0]->ncols;
      nrows = cimg->C[0]->nrows;
      n  = ncols*nrows;
      sg = Create(ncols, nrows, A);
      sg->Wmax = 0;
      
      for(p=0; p<n; p++){
	u.x = p%ncols; //PixelX(img,p);
	u.y = p/ncols; //PixelY(img,p);
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  if(gft::Image32::IsValidPixel(cimg->C[0],v.x,v.y)){
	    q = v.x + ncols*v.y;
	    dr = abs(cimg->C[0]->data[p] - cimg->C[0]->data[q]);
	    dg = abs(cimg->C[1]->data[p] - cimg->C[1]->data[q]);
	    db = abs(cimg->C[2]->data[p] - cimg->C[2]->data[q]);
	    
	    //weight = MAX(dr, MAX(dg, db));
	    weight = ROUND(sqrtf(dr*dr + dg*dg + db*db));
	    (sg->n_link[p])[i] = weight;
	    
	    if(weight > sg->Wmax) sg->Wmax = weight;
	  }
	  else{
	    (sg->n_link[p])[i] = NIL;
	  }
	}
      }
      sg->Wmax++;
      for(p=0; p<n; p++){
	(sg->n_link[p])[0] = sg->Wmax;
	for(i=1; i<A->n; i++){
	  if((sg->n_link[p])[i]==NIL) continue;
	  (sg->n_link[p])[i] = sg->Wmax - (sg->n_link[p])[i];
	}
      }
      sg->type = CAPACITY;
      return sg;
    }


    
    SparseGraph *WeightedMean(SparseGraph *G1, 
			      SparseGraph *G2,
			      float w2){
      SparseGraph *sg;  
      int i,p,n,ncols,nrows;
      int weight,v1,v2;
      gft::Pixel u,v;
      gft::AdjRel::AdjRel *A;
      
      if(G1->A->n != G2->A->n){
	gft::Error((char*)"Incompatible graphs",
		   (char*)"WeightedMean");
      }
      A = gft::AdjRel::Clone(G1->A);
      ncols     = G1->ncols; 
      nrows     = G1->nrows;
      n         = ncols*nrows;
      sg        = Create(ncols, nrows, A);
      sg->type  = CAPACITY;
      sg->Wmax  = 0;

      ChangeType(G1, CAPACITY);
      ChangeType(G2, CAPACITY);

      for(p=0; p<n; p++){
	u.x = p%ncols;
	u.y = p/ncols;
	weight = 0;
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  if(v.x>=0 && v.y>=0 && v.x<ncols && v.y<nrows){
	    
	    v1 = (G1->n_link[p])[i];
	    v2 = (G2->n_link[p])[i];
	    
	    weight = ROUND(w2*v2 +(1.0-w2)*v1);
	    
	    (sg->n_link[p])[i] = weight;
	    if(weight>sg->Wmax) sg->Wmax = weight;
	  }
	  else{
	    (sg->n_link[p])[i] = NIL;
	  }
	}
      }
      //sg->Wmax = MAX(sg->Wmax, _MAXWEIGHT_);
      
      return sg;
    }


    SparseGraph *ByLevel(gft::Image32::Image32 *obj, int T, float r){
      SparseGraph *sg=NULL;
      gft::AdjRel::AdjRel *A = gft::AdjRel::Circular(r);
      int p,q,i,n,max,min;
      int weight;
      gft::Pixel u,v;
      
      n  = obj->ncols*obj->nrows;
      sg = Create(obj->ncols, obj->nrows, A);
      sg->type = CAPACITY;
      sg->Wmax = 0;
      
      for(p=0; p<n; p++){
	u.x = p%obj->ncols; 
	u.y = p/obj->ncols; 
	for(i=1; i<A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  if(gft::Image32::IsValidPixel(obj,v.x,v.y)){
	    q = v.x + obj->ncols * v.y;
	    
	    max = MAX(obj->data[p], obj->data[q]);
	    min = MIN(obj->data[p], obj->data[q]);
	    if(max>=T && min<T) weight = 0;
	    else   	        weight = 1024; //_MAXWEIGHT_;
	    
	    (sg->n_link[p])[i] = weight;
	  }
	  else{
	    (sg->n_link[p])[i] = NIL;
	  }
	}
      }
      sg->Wmax = 1024; //_MAXWEIGHT_;
      return sg;
    }
    


    SparseGraph *ByExclusion(gft::Image32::Image32 *obj,
			     gft::Image32::Image32 *bkg, float r){
      SparseGraph *sg=NULL;
      gft::AdjRel::AdjRel *A = gft::AdjRel::Circular(r);
      int p,q,i,n,omin,bmin;
      int weight;
      gft::Pixel u,v;
      
      n  = obj->ncols*obj->nrows;
      sg = Create(obj->ncols, obj->nrows, A);
      sg->type = CAPACITY;
      sg->Wmax = 0;
      
      for(p = 0; p < n; p++){
	u.x = p%obj->ncols; 
	u.y = p/obj->ncols; 
	for(i = 1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];
	  if(gft::Image32::IsValidPixel(obj,v.x,v.y)){
	    q = v.x + obj->ncols * v.y;
	    
	    omin = MIN(obj->data[p], obj->data[q]);
	    bmin = MIN(bkg->data[p], bkg->data[q]);
	    weight = MAX(omin, bmin);
	    
	    (sg->n_link[p])[i] = weight;
	    //if(weight > sg->Wmax) sg->Wmax = weight;
	  }
	  else{
	    (sg->n_link[p])[i] = NIL;
	  }
	}
      }
      sg->Wmax = 1024; //_MAXWEIGHT_;
      return sg;
    }


    void LinearStretch(SparseGraph *sg, 
		       int f1, int f2, 
		       int g1, int g2){
      gft::AdjRel::AdjRel *A;
      int p,n,i,v;
      float a=1.0;
      
      n = sg->ncols*sg->nrows;
      if (f1 != f2) 
	a = (float)(g2-g1)/(float)(f2-f1);
      else
	gft::Error((char *)"Invalid input value",
		   (char *)"SparseGraph::LinearStretch");
      
      A = sg->A;
      for(p = 0; p < n; p++){
	for(i = 1; i < A->n; i++){
	  if((sg->n_link[p])[i]==NIL) continue;
	  v = (sg->n_link[p])[i];
	  
	  if(v < f1)
	    (sg->n_link[p])[i] = g1;
	  else if(v > f2)
	    (sg->n_link[p])[i] = g2;
	  else
	    (sg->n_link[p])[i] = MIN(ROUND(a*(v-f1)+g1),g2);
	}
      }
      sg->Wmax = g2;
    }

    

    int get_edge_index(int p, int q, SparseGraph *g){
      AdjRel::AdjRel *A;
      int ncols,dx,dy,i;
      
      A = g->A;
      ncols = g->ncols;
      dx = q%ncols - p%ncols;
      dy = q/ncols - p/ncols;
      
      for(i=1; i<A->n; i++){
	if(dx == A->dx[i] && dy == A->dy[i])
	  return i;
      }
      return -1;
    }
    

    SparseGraph   *Create(int ncols, int nrows, AdjRel::AdjRel *A){
      SparseGraph *sg;
      int p,n,i;
      
      n = ncols*nrows;
      sg = (SparseGraph *)calloc(1, sizeof(SparseGraph));
      if(sg == NULL)
	gft::Error((char *)MSG1,(char *)"SparseGraph::Create");
      sg->A = A;
      sg->ncols = ncols;
      sg->nrows = nrows;
      sg->n_link = (int **)calloc(n,sizeof(int *));
      if(sg->n_link == NULL)
	gft::Error((char *)MSG1,(char *)"SparseGraph::Create");
      for(p=0; p<n; p++){
	sg->n_link[p] = gft::AllocIntArray(A->n);
	for(i=0; i<A->n; i++)
	  (sg->n_link[p])[i] = NIL;
      }
      sg->Wmax = 0;
      
      return sg;
    }


    void    Destroy(SparseGraph **sg){
      SparseGraph *aux;  
      int p,n;
      
      aux = *sg;
      if(aux != NULL){
	n = aux->ncols*aux->nrows;
	AdjRel::Destroy(&(aux->A));
	for(p=0; p<n; p++){
	  free(aux->n_link[p]);
	}
	free(aux->n_link);
	free(*sg);
	*sg = NULL;
      }
    }
    

    SparseGraph   *Clone(SparseGraph *sg){
      SparseGraph *clone = NULL;
      int p,n,i;
      
      n = sg->ncols*sg->nrows;
      
      clone = Create(sg->ncols, 
		     sg->nrows, 
		     AdjRel::Clone(sg->A));
      for(p=0; p<n; p++){
	for(i=1; i<(sg->A)->n; i++){
	  (clone->n_link[p])[i] = (sg->n_link[p])[i];
	}
      }
      clone->Wmax = sg->Wmax;
      clone->type = sg->type;

      return clone;
    }


    SparseGraph *ReadFromTxt(char *filename){
      SparseGraph *sg;
      gft::AdjRel::AdjRel *A = gft::AdjRel::Circular(1.0);
      int ncols,nrows,i,j,p,q,k,val,type;
      FILE *fp;
      
      fp = fopen(filename,"r");
      if(fp == NULL)
	gft::Error((char *)"Cannot open file",
		   (char *)"ReadSparseGraphFromTxt");
      fscanf(fp,"%d %d %d\n",&type,&ncols,&nrows);
      sg = Create(ncols, nrows, A);
      sg->type = type;
      sg->Wmax = 0;
      
      for(i=0; i<nrows*2-1; i++){
	if(i%2==0) //Horizontal edges.
	  for(j=0; j<ncols-1; j++){
	    p = j + ncols*(i/2);
	    q = p + 1;
	    fscanf(fp,"%d",&val);
	    k = get_edge_index(p, q, sg);
	    (sg->n_link[p])[k] = val;
	    k = get_edge_index(q, p, sg);
	    (sg->n_link[q])[k] = val;
	    if(val>sg->Wmax)
	      sg->Wmax = val;
	  }
	else       //Vertical edges.
	  for(j=0; j<ncols; j++){
	    p = j + ncols*(i-1)/2;
	    q = j + ncols*(i+1)/2;
	    fscanf(fp,"%d",&val);
	    k = get_edge_index(p, q, sg);
	    (sg->n_link[p])[k] = val;
	    k = get_edge_index(q, p, sg);
	    (sg->n_link[q])[k] = val;
	    if(val>sg->Wmax)
	      sg->Wmax = val;
	  }
      }
      
      fclose(fp);
      
      return sg;
    }


    void         Write2Txt(SparseGraph *sg, 
			   char *filename){
      int ncols,nrows,i,j,p,q,k,val;
      FILE *fp;
      
      if(sg->A->n!=5)
	gft::Error((char *)"Current implementation only supports adjacency-4",
		   (char *)"WriteSparseGraph2Txt");
      
      fp = fopen(filename,"w");
      if(fp == NULL)
	gft::Error((char *)"Cannot open file",
		   (char *)"WriteSparseGraph2Txt");
      
      ncols = sg->ncols;
      nrows = sg->nrows;
      fprintf(fp,"%d %d %d\n",sg->type,ncols,nrows);
      
      for(i=0; i<nrows*2-1; i++){
	if(i%2==0) //Horizontal edges.
	  for(j=0; j<ncols-1; j++){
	    p = j + ncols*(i/2);
	    q = p + 1;
	    k = get_edge_index(p, q, sg);
	    val = (sg->n_link[p])[k];
	    fprintf(fp,"%d ",val);
	  }
	else       //Vertical edges.
	  for(j=0; j<ncols; j++){
	    p = j + ncols*(i-1)/2;
	    q = j + ncols*(i+1)/2;
	    k = get_edge_index(p, q, sg);
	    val = (sg->n_link[p])[k];
	    fprintf(fp,"%d ",val);
	  }
	fprintf(fp,"\n");
      }
      
      fclose(fp);
    }


    Image32::Image32 *ArcWeightImage(SparseGraph *sg){
      int i,p,n,weight;
      Image32::Image32 *arcw;
      AdjRel::AdjRel *A;
      
      A    = sg->A;
      arcw = Image32::Create(sg->ncols,sg->nrows);
      n    = arcw->ncols*arcw->nrows;
      
      if(sg->type == CAPACITY){
	for(p=0; p<n; p++){
	  arcw->data[p] = sg->Wmax;
	  for(i=1; i<A->n; i++){
	    weight = (sg->n_link[p])[i];
	    if(weight==NIL) continue;
	    arcw->data[p] = MIN(weight, arcw->data[p]);
	  }
	}
      }
      else{ //DISSIMILARITY
	for(p=0; p<n; p++){
	  arcw->data[p] = 0;
	  for(i=1; i<A->n; i++){
	    weight = (sg->n_link[p])[i];
	    if(weight==NIL) continue;
	    arcw->data[p] = MAX(weight, arcw->data[p]);
	  }
	}
      }
      return(arcw);
    }
    

    void Orient2Digraph(SparseGraph *sg, 
			Image32::Image32 *img, float per){
      AdjRel::AdjRel *A;
      int p,q,n,i,val,max = INT_MIN;
      int u_x,u_y,v_x,v_y;
      
      n = sg->ncols*sg->nrows;
      A = sg->A;
      
      for(p=0; p<n; p++){
	u_x = p%img->ncols;
	u_y = p/img->ncols;
	
	for(i=0; i<A->n; i++){
	  if((sg->n_link[p])[i]==NIL) continue;
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(!IsValidPixel(img,v_x,v_y)) continue;
	  q = v_x + img->ncols*v_y;
	  
	  val = (sg->n_link[p])[i];
	  if(img->data[p] > img->data[q]){
	    //val += ROUND(val*(per/100.0)); 
	    val = ROUND(val*(1.0 + per/100.0));
	  }
	  else if(img->data[p] < img->data[q]){
	    //val -= ROUND(val*(per/100.0)); 
	    val = ROUND(val*(1.0 - per/100.0));
	  }
	  (sg->n_link[p])[i] = val;
	  if(val > max) max = val;
	}
      }
      sg->Wmax = max;
    }
    

    void Transpose(SparseGraph *sg){
      AdjRel::AdjRel *A;
      int p,q,n,i,j,val_pq,val_qp;
      int u_x,u_y,v_x,v_y;
      
      n = sg->ncols*sg->nrows;
      A = sg->A;
      
      for(p=0; p<n; p++){
	u_x = p%sg->ncols;
	u_y = p/sg->ncols;
	
	for(i=0; i<A->n; i++){
	  if((sg->n_link[p])[i]==NIL) continue;
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  
	  if(v_x < 0 || v_x >= sg->ncols ||
	     v_y < 0 || v_y >= sg->nrows)
	    continue;
	  
	  q = v_x + sg->ncols*v_y;
	  if(p > q){
	    j = get_edge_index(q, p, sg);
	    
	    val_pq = (sg->n_link[p])[i];
	    val_qp = (sg->n_link[q])[j];
	    
	    (sg->n_link[p])[i] = val_qp;
	    (sg->n_link[q])[j] = val_pq;
	  }
	}
      }
    }
    


    void Orient2DigraphInner(SparseGraph *sg, Image32::Image32 *P_sum){
      AdjRel::AdjRel *A;
      int p,q,n,i,val,max = INT_MIN;
      int u_x,u_y,v_x,v_y;
      
      n = sg->ncols*sg->nrows;
      A = sg->A;
      
      for(p=0; p<n; p++){
	u_x = p%sg->ncols;
	u_y = p/sg->ncols;
	
	for(i=1; i<A->n; i++){
	  if((sg->n_link[p])[i]==NIL) continue;
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(!IsValidPixel(P_sum,v_x,v_y)) continue;
	  q = v_x + sg->ncols*v_y;
	  
	  if (P_sum->data[q] == p) {
	    sg->n_link[p][i] = 0;  //w(P_sum(q),q)=0
	  }
	  
	  val = (sg->n_link[p])[i]; 

	  if(val > max) max = val;
	}
      }
      sg->Wmax = max;
    }

    void Orient2DigraphOuter(SparseGraph *sg, Image32::Image32 *P_sum){
      AdjRel::AdjRel *A;
      int p,q,n,i,val,max = INT_MIN;
      int u_x,u_y,v_x,v_y;
      
      n = sg->ncols*sg->nrows;
      A = sg->A;

      for(p=0; p<n; p++){
	u_x = p%sg->ncols;
	u_y = p/sg->ncols;
	
	for(i=1; i<A->n; i++){
	  if((sg->n_link[p])[i]==NIL) continue;
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(!IsValidPixel(P_sum,v_x,v_y)) continue;
	  q = v_x + sg->ncols*v_y;
	  
	  if (P_sum->data[p] == q) {
	    //j = get_edge_index(q, p, sg);
	    sg->n_link[p][i] = 0;  // w(q,P_sum(q)) = 0
	  }
	  val = (sg->n_link[p])[i]; 
	  
	  if(val > max) max = val;
	}
      }
      sg->Wmax = max;
    }
    

    void  SuppressZeroWeightedArcs(SparseGraph *g){
      int i,p,n,wzero = 0;
      AdjRel::AdjRel *A;
      
      A = g->A;
      n = g->ncols*g->nrows;
      for(p=0; p<n && wzero==0; p++){
	for(i=1; i<A->n; i++){
	  if((g->n_link[p])[i] == NIL)
	    continue;
	  
	  if((g->n_link[p])[i] == 0){
	    wzero = 1;
	    break;
	  }
	}
      }
      
      if(wzero==1){
	for(p=0; p<n; p++){
	  for(i=1; i<A->n; i++){
	    if((g->n_link[p])[i] == NIL)
	      continue;
	    (g->n_link[p])[i] += 1; 
	  }
	}
	g->Wmax++;
      }
    }
    


    void ChangeType(SparseGraph *sg, 
		    int type){
      int i,p,n;
      AdjRel::AdjRel *A;
      
      //printf("Change from %d to type: %d\n",sg->type, type);
      
      if(sg->type==type)
	return;
      
      A = sg->A;
      n = sg->ncols*sg->nrows;
      for(p=0; p<n; p++){
	//(sg->n_link[p])[0] = sg->Wmax or 0;
	for(i=1; i<A->n; i++){
	  if((sg->n_link[p])[i]==NIL) continue;
	  (sg->n_link[p])[i] = sg->Wmax - (sg->n_link[p])[i];
	}
      }
      if(sg->type==DISSIMILARITY)
	sg->type = CAPACITY;
      else
	sg->type = DISSIMILARITY;
    }
    

    // Increasing transformation.
    void PowSparseGraph(SparseGraph *sg, 
			int power, int max){
      AdjRel::AdjRel *A;
      int p,n,i;
      long double v,m;
      
      if(sg->Wmax<=1) return;
      if(power<=0)    return;
      
      n = sg->ncols*sg->nrows;
      A = sg->A;
      
      for(p=0; p<n; p++){
	for(i=0; i<A->n; i++){
	  if((sg->n_link[p])[i]==NIL) continue;
	  v = (double)((sg->n_link[p])[i]);
	  m = powl((v/sg->Wmax),power);
	  (sg->n_link[p])[i] = ROUND(m*(max-1))+1;
	}
      }
      sg->Wmax = max;
    }
    
    
    
    SparseGraph   *Convert2Symmetric(SparseGraph *g,
				     int method){
      SparseGraph *sg;
      AdjRel::AdjRel *A;
      int i,k,p,q,n,valpq,valqp,val=0;
      int ncols,nrows;
      int u_x,u_y,v_x,v_y;
      
      A     = g->A;
      ncols = g->ncols;
      nrows = g->nrows;
      n = ncols*nrows;
      
      sg = Create(ncols, nrows,
		  AdjRel::Clone(A));
      sg->type = g->type;
      sg->Wmax = 0;
      for(p=0; p<n; p++){
	u_x = p%ncols;
	u_y = p/ncols;
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(v_x>=0 && v_x<ncols && 
	     v_y>=0 && v_y<nrows){
	    q = v_x + v_y*ncols;
	    k = get_edge_index(q, p, g);
	    
	    valpq = (g->n_link[p])[i];
	    valqp = (g->n_link[q])[k];
	    switch(method){
	    case 0:
	      val = (valpq + valqp)/2;
	      break;
	    case 1:
	      val = MAX(valpq, valqp);
	      break;
	    case 2:
	      val = MIN(valpq, valqp);
	      break;
	    }
	    (sg->n_link[p])[i] = val;
	    (sg->n_link[q])[k] = val;
	    if(val>sg->Wmax)
	      sg->Wmax = val;
	  }
	}
      }
      
      return sg;
    }
    
    
    
    
  } /*end SparseGraph namespace*/
} /*end gft namespace*/


