
#include "gft_radiometric.h"


namespace gft{

  namespace Image32{

    
    gft::Curve::Curve *Histogram(Image32 *img){
      int i,p,n,nbins;
      gft::Curve::Curve *hist=NULL;
      
      nbins = GetMaxVal(img)+1;
      hist  = gft::Curve::Create(nbins);
      n     = img->ncols*img->nrows;
      for (p=0; p < n; p++)
	hist->Y[img->data[p]]++;
      for (i=0; i < nbins; i++) 
	hist->X[i] = i;
      
      return(hist);
    }


    gft::Curve::Curve *NormHistogram(Image32 *img){
      int i,sum;
      gft::Curve::Curve *hist=NULL,*nhist=NULL;
      
      hist  = Histogram(img);
      sum   = img->ncols*img->nrows;
      nhist = gft::Curve::Create(hist->n);
      for (i=0; i < nhist->n;i++){
	nhist->Y[i] = hist->Y[i]/sum;
	nhist->X[i] = hist->X[i];
      }
      gft::Curve::Destroy(&hist);
      return(nhist);
    }


    gft::Curve::Curve *NormalizeHistogram(gft::Curve::Curve *hist){
      gft::Curve::Curve *nhist;
      double sum;
      int i;
      
      nhist = gft::Curve::Clone(hist);
      sum = 0.0;
      for(i=0; i<nhist->n; i++)
	sum += nhist->Y[i];
      for(i=0; i<nhist->n; i++)
	nhist->Y[i] /= sum;
      
      return (nhist);
    }


    gft::Curve::Curve  *RemoveEmptyBins(gft::Curve::Curve *hist){
      gft::Curve::Curve *C;
      int n,i,j;
      n = 0;
      for(i=0; i<hist->n; i++){
	if(hist->Y[i]>0.0) n++;
      }
      j = 0;
      C = gft::Curve::Create(n);
      for(i=0; i<hist->n; i++){
	if(hist->Y[i]>0.0){
	  C->Y[j] = hist->Y[i];
	  C->X[j] = hist->X[i];      
	  j++;
	}
      }
      return C;
    }
   

    Image32 *LinearStretch(Image32 *img, int omin, int omax, int nmin, int nmax){
      Image32 *simg=NULL;
      int p,n;
      float a;
      
      simg = Create(img->ncols,img->nrows);
      n    = img->ncols*img->nrows;
      if (omin != omax) 
	a = (float)(nmax-nmin)/(float)(omax-omin);
      else
	a = INT_MAX;
      
      for (p=0; p < n; p++){
	if (img->data[p] < omin)
	  simg->data[p] = nmin;
	else 
	  if (img->data[p] > omax)
	    simg->data[p] = nmax;
	  else {
	    if (a != INT_MAX)	  
	      simg->data[p] = (int)(a*(img->data[p]-omin)+nmin);
	    else{
	      simg->data[p] = nmax;
	    }   
	  }
      }
      return(simg);
    }


  } //end Image32 namespace


} //end gft namespace


