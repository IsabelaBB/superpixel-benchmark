
#include "gft_morphology.h"

namespace gft{
  namespace Image32{


    Image32 *Dilate(Image32 *img, gft::AdjRel::AdjRel *A){
      Image32 *dil=NULL;
      int p,q,n,max,i,xp,yp,xq,yq;
      
      n = img->ncols*img->nrows;
      dil = Create(img->ncols, img->nrows);
      for(p = 0; p < n; p++){
	xp = p%img->ncols;
	yp = p/img->ncols;
	
	max = INT_MIN;
	for(i=0; i < A->n; i++){
	  xq = xp + A->dx[i];
	  yq = yp + A->dy[i];
	  if(IsValidPixel(img, xq, yq)){
	    q = xq + yq*img->ncols;
	    if(img->data[q] > max)
	      max = img->data[q];
	  }
	}
	dil->data[p] = max;
      }
      return(dil);
    }

    
    Image32 *Erode(Image32 *img, gft::AdjRel::AdjRel *A){
      Image32 *ero=NULL;
      int p,q,n,min,i,xp,yp,xq,yq;
      
      n = img->ncols*img->nrows;
      ero = Create(img->ncols, img->nrows);
      for(p = 0; p < n; p++){
	xp = p%img->ncols;
	yp = p/img->ncols;
	
	min = INT_MAX;
	for(i=0; i < A->n; i++){
	  xq = xp + A->dx[i];
	  yq = yp + A->dy[i];
	  if(IsValidPixel(img, xq, yq)){
	    q = xq + yq*img->ncols;
	    if(img->data[q] < min)
	      min = img->data[q];
	  }
	}
	ero->data[p] = min;
      }
      return(ero);
    }



    Image32 *MorphGrad(Image32 *img, gft::AdjRel::AdjRel *A){
      Image32 *dil=NULL,*ero=NULL,*grad=NULL;
      int p,n;
      n = img->ncols*img->nrows;
      grad = Create(img->ncols, img->nrows);
      dil  = Dilate(img,A);
      ero  = Erode(img,A);
      /*grad = Diff(dil,ero);*/
      for(p = 0; p < n; p++){
	grad->data[p] = dil->data[p] -  ero->data[p];
      }
      Destroy(&dil);
      Destroy(&ero);
      return(grad);
    }



    /* It assumes that the next operation is a dilation, but it may
       be an erosion if you remove comments below. */

    Image32 *ErodeBin(Image32 *bin, gft::Set::Set **seed, float radius){
      Image32 *ero=NULL,*boundr=NULL,*dil=NULL;
      gft::Pixel u,v,w;
      Image32 *cost=NULL,*root;
      gft::PQueue32::PQueue32 *Q=NULL;
      int i,p,q,n,sz;
      int *sq=NULL,tmp=INT_MAX,dx,dy;
      gft::AdjRel::AdjRel *A=NULL;
      float dist;

      /* Compute seeds */
      
      if (*seed == NULL) {
	A      = gft::AdjRel::Circular(1.0);
	dil    = Dilate(bin,A);
	//boundr = Diff(dil,bin);
	boundr = gft::Image32::Create(bin->ncols, bin->nrows);
	n = boundr->ncols*boundr->nrows;
	for(p = 0; p < n; p++){
	  boundr->data[p] = dil->data[p] -  bin->data[p];
	}
	Destroy(&dil);
	gft::AdjRel::Destroy(&A);

	for (p=0; p < n; p++)
	  if (boundr->data[p]==1)
	    gft::Set::Insert(seed, p);
	Destroy(&boundr);    
      }
      
      /* Erode image */
      
      ero  = gft::Image32::Clone(bin);
      dist = (radius*radius);
      A  = gft::AdjRel::Circular(1.5);
      n  = MAX(ero->ncols,ero->nrows);
      sq = gft::AllocIntArray(n);
      for (i=0; i < n; i++) 
	sq[i]=i*i;
      
      cost = Create(ero->ncols,ero->nrows);
      root = Create(ero->ncols,ero->nrows);
      Set(cost, INT_MAX);
      n    = ero->ncols*ero->nrows;
      sz   = gft::AdjRel::GetFrameSize(A);
      Q    = gft::PQueue32::Create(2*sz*(sz+ero->ncols+ero->nrows), n, cost->data);
      
      while (*seed != NULL){
	p = gft::Set::Remove(seed);
	cost->data[p]=0;
	root->data[p]=p;
	gft::PQueue32::InsertElem(&Q, p);
      }
      
      while(!gft::PQueue32::IsEmpty(Q)) {
	p = gft::PQueue32::RemoveMinFIFO(Q);
	if (cost->data[p] <= dist){
	  
	  ero->data[p] = 0;

	  /* Seeds for erosion if we wanted to compute sequences of erosions
	     
	     if (((sq[Dx->data[p]+1]+sq[Dy->data[p]]) > dist)||
	     ((sq[Dx->data[p]]+sq[Dy->data[p]+1]) > dist)){
	     InsertSet(seed,p);
	     }
	     
	  */
	  
	  u.x = p%ero->ncols;
	  u.y = p/ero->ncols;
	  w.x = root->data[p]%ero->ncols;
	  w.y = root->data[p]/ero->ncols;
	  for (i=1; i < A->n; i++){
	    v.x = u.x + A->dx[i];
	    v.y = u.y + A->dy[i];
	    if (IsValidPixel(ero,v.x,v.y)){
	      q = v.x + ero->ncols*v.y;
	      if ((cost->data[p] < cost->data[q])&&(ero->data[q]==1)){
		dx  = abs(v.x-w.x);
		dy  = abs(v.y-w.y);
		tmp = sq[dx] + sq[dy];
		if (tmp < cost->data[q]){
		  if (cost->data[q] == INT_MAX){
		    cost->data[q] = tmp;
		    gft::PQueue32::InsertElem(&Q, q);
		  }
		  else
		    gft::PQueue32::UpdateElem(&Q, q, tmp);

		  root->data[q] = root->data[p];
		}
	      }
	    }
	  }
	} else {  /* Seeds for dilation */
	  gft::Set::Insert(seed, p);
	}
      }
      
      free(sq);
      gft::PQueue32::Destroy(&Q);
      Destroy(&root);
      Destroy(&cost);
      gft::AdjRel::Destroy(&A);
      
      return(ero);
    }

    
    void SupRec_Watershed(gft::AdjRel::AdjRel *A, 
			  Image32 *I, Image32 *J, 
			  Image32 *L, Image32 *V){
      gft::PQueue32::PQueue32 *Q=NULL;
      Image32 *P;
      int Jmax,n,p,q,tmp,i,l = 1;
      int xp,yp,xq,yq;
      
      Jmax = GetMaxVal(J);
      n = I->ncols*I->nrows;
      P = gft::Image32::Create(I->ncols, I->nrows);
      Q = gft::PQueue32::Create(Jmax+2, n, V->data);

      gft::Image32::Set(L, 0);
      for(p = 0; p < n; p++){
	P->data[p] = NIL;
	V->data[p] = J->data[p] + 1;
	gft::PQueue32::InsertElem(&Q, p);
      }
      
      while(!gft::PQueue32::IsEmpty(Q)) {
	p = gft::PQueue32::RemoveMinFIFO(Q);
	
	if(P->data[p]==NIL){
	  V->data[p] = J->data[p];
	  L->data[p] = l;
	  l++;
	}
	
	xp = p%I->ncols;
	yp = p/I->ncols;
	for(i = 1; i < A->n; i++){
	  xq = xp + A->dx[i];
	  yq = yp + A->dy[i];
	  if(gft::Image32::IsValidPixel(I, xq, yq)){
	    q = xq + yq*I->ncols;
	    if(Q->L.elem[q].color != BLACK){
	      tmp = MAX(V->data[p], I->data[q]);
	      
	      if(tmp < V->data[q]){
		if(Q->L.elem[q].color == GRAY)
		  gft::PQueue32::RemoveElem(Q, q);
		L->data[q] = L->data[p];
		P->data[q] = p;
		V->data[q] = tmp;
		gft::PQueue32::InsertElem(&Q, q);
	      }
	    }
	  }
	}
      }
      gft::Image32::Destroy(&P);
      gft::PQueue32::Destroy(&Q);
    }



    /* Devolve o representante de p */
    int Find_(Image32 *R, int p){
      if(R->data[p] == p)
	return p;
      else{
	R->data[p] = Find_(R, R->data[p]);
	return R->data[p];
      }
    }


    /*Removes all background connected components from the stack
      of binary images of I whose area (number of pixels) is less than 
      a threshold and outputs a simplified image.*/
    Image32 *AreaClosing(gft::AdjRel::AdjRel *A,
			 Image32 *I, int T){
      Image32 *J,*V;
      Image32 *Ar; /* Areas para os representantes */
      Image32 *R; /* Imagem de representantes */
      Image32 *P; /* Imagem de predecessores */
      gft::PQueue32::PQueue32 *Q;
      int p,q,n,Imax,rp,rq,tmp;
      int xp,yp,xq,yq,i;
      //Image32 *teste;
      
      Imax = gft::Image32::GetMaxVal(I);
      n  = I->ncols*I->nrows;
      V  = gft::Image32::Create(I->ncols, I->nrows);
      J  = gft::Image32::Create(I->ncols, I->nrows);
      Ar = gft::Image32::Create(I->ncols, I->nrows);
      R  = gft::Image32::Create(I->ncols, I->nrows);
      P  = gft::Image32::Create(I->ncols, I->nrows);
      
      /*teste  = CreateImage(I->ncols, I->nrows);*/
      
      Q  = gft::PQueue32::Create(Imax+2, n, V->data);
      //SetTieBreak(Q, LIFOBREAK);
      
      for(p = 0; p < n; p++){
	P->data[p] = NIL;
	V->data[p] = I->data[p] + 1;
	R->data[p] = p;
	Ar->data[p] = 0;
	J->data[p] = I->data[p];
	gft::PQueue32::InsertElem(&Q, p);
      }
      
      while(!gft::PQueue32::IsEmpty(Q)){
	p = gft::PQueue32::RemoveMinLIFO(Q);
	rp = Find_(R, p);
	
	if(Ar->data[rp] <= T && J->data[rp] < I->data[p])
	  J->data[rp] = I->data[p];
	
	Ar->data[rp]++;
	xp = p%I->ncols;
	yp = p/I->ncols;
	for(i = 1; i < A->n; i++){
	  xq = xp + A->dx[i];
	  yq = yp + A->dy[i];
	  if(gft::Image32::IsValidPixel(I, xq, yq)){
	    q = xq + yq*I->ncols;
	    /*
	      if(V->data[q] > V->data[p]){
	    */
	    if(Q->L.elem[q].color != BLACK){
	      tmp = MAX(V->data[p], I->data[q]);
	      if(tmp <= V->data[q]){
		gft::PQueue32::RemoveElem(Q, q);
		V->data[q] = tmp;
		P->data[q] = p;
		R->data[q] = rp;
		gft::PQueue32::InsertElem(&Q, q);
	      }
	    }
	    else{/* if(Q->L.elem[q].color != GRAY){*/
	      rq = Find_(R, q);
	      if(rp != rq){
		if(Ar->data[rq] <= T && J->data[rq] < I->data[p])
		  J->data[rq] = I->data[p];
		if(Ar->data[rp] < Ar->data[rq]){
		  tmp = rp;
		  rp = rq;
		  rq = tmp;
		}
		R->data[rq] = rp;
		Ar->data[rp] += Ar->data[rq];
		/*
		  if(P->data[p] == NIL && P->data[q] == NIL)
		  P->data[rq] = rp;
		*/
	      }
	    }
	  }
	}
      }
      gft::PQueue32::Reset(Q);
      //SetRemovalPolicy(Q, MAXVALUE);
      
      for(p = 0; p < n; p++){
	V->data[p] = J->data[p];
	if(P->data[p] == NIL){
	  gft::PQueue32::InsertElem(&Q, p);
	  //teste->data[p] = 1;
	}
      }
      while(!gft::PQueue32::IsEmpty(Q)){
	p = gft::PQueue32::RemoveMaxLIFO(Q);
	xp = p%I->ncols;
	yp = p/I->ncols;
	for(i = 1; i < A->n; i++){
	  xq = xp + A->dx[i];
	  yq = yp + A->dy[i];
	  if(gft::Image32::IsValidPixel(I, xq, yq)){
	    q = xq + yq*I->ncols;
	    if(V->data[p] > V->data[q]){
	      if(Q->L.elem[q].color == GRAY)
		gft::PQueue32::RemoveElem(Q, q);
	      V->data[q] = V->data[p];
	      gft::PQueue32::InsertElem(&Q, q);
	    }
	  }
	}
      }
      
      //WriteImage(teste, "teste.pgm");
      
      gft::Image32::Destroy(&J);
      gft::Image32::Destroy(&R);
      gft::Image32::Destroy(&P);
      gft::Image32::Destroy(&Ar);
      gft::PQueue32::Destroy(&Q);
      return V;
    }



    Image32 *CloseHoles(Image32 *img){
      gft::PQueue32::PQueue32 *Q;
      gft::Image32::Image32 *V;
      gft::AdjRel::AdjRel *A;
      int i,j, p,q, xp,yp,xq,yq,tmp,Imax;
      V = gft::Image32::Create(img);
      gft::Image32::Set(V, INT_MAX);
      Imax = gft::Image32::GetMaxVal(img);
      A = gft::AdjRel::Circular(1.0);
      Q = gft::PQueue32::Create(Imax+2, img->n, V->data);

      for(j = 0; j < img->ncols; j++){
	p = j;
	V->data[p] = img->data[p];
	gft::PQueue32::InsertElem(&Q, p);

	p = j + (img->nrows-1)*img->ncols;
	V->data[p] = img->data[p];
	gft::PQueue32::InsertElem(&Q, p);
      }
      for(i = 1; i < img->nrows-1; i++){
	p = i*img->ncols;
	V->data[p] = img->data[p];
	gft::PQueue32::InsertElem(&Q, p);

	p = i*img->ncols + img->ncols-1;
	V->data[p] = img->data[p];
	gft::PQueue32::InsertElem(&Q, p);
      }

      while(!gft::PQueue32::IsEmpty(Q)){
	p = gft::PQueue32::RemoveMinFIFO(Q);
	xp = p%img->ncols;
	yp = p/img->ncols;
	for(i = 1; i < A->n; i++){
	  xq = xp + A->dx[i];
	  yq = yp + A->dy[i];
	  if(gft::Image32::IsValidPixel(img, xq, yq)){
	    q = xq + yq*img->ncols;
	    tmp = MAX(V->data[p], img->data[q]);
	    if(tmp < V->data[q]){
	      if(Q->L.elem[q].color == GRAY)
		gft::PQueue32::RemoveElem(Q, q);
	      V->data[q] = tmp;
	      gft::PQueue32::InsertElem(&Q, q);
	    }
	  }
	}
      }
      gft::AdjRel::Destroy(&A);
      gft::PQueue32::Destroy(&Q);
      return V;
    }
    
    
  } /*end Image32 namespace*/
} /*end gft namespace*/

