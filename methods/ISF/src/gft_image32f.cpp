
#include "gft_image32f.h"

namespace gft{
  namespace Image32f{

    Image32f *Create(Image32f *img){
      Image32f *nimg=NULL;
      nimg = Create(img->ncols, img->nrows);
      nimg->dx = img->dx;
      nimg->dy = img->dy;
      return nimg;
    }
    
    Image32f *Create(int ncols, int nrows){
      Image32f *img=NULL;
      int i;
      img = (Image32f *) calloc(1,sizeof(Image32f));
      if (img == NULL){
	gft::Error((char *)MSG1,(char *)"Image32f::Create");
      }
      img->data  = gft::AllocFloatArray(nrows*ncols);
      img->ncols = ncols;
      img->nrows = nrows;
      img->n = nrows*ncols;
      img->dx = img->dy = 1.0;
      
      img->array = (float**)malloc(nrows*sizeof(float*));
      if(img->array == NULL){
	gft::Error((char *)MSG1,(char *)"Image32f::Create");
      }
      for(i = 0; i < nrows; i++){
	img->array[i] = (img->data + i*ncols);
      }
      return(img);
    }
    

    void Destroy(Image32f **img){
      Image32f *aux;
      if(img != NULL){
	aux = *img;
	if (aux != NULL){
	  if(aux->data !=  NULL) gft::FreeFloatArray(&aux->data);
	  if(aux->array != NULL) free(aux->array);
	  free(aux);
	  *img = NULL;
	}
      }
    }
    
    
    Image32f *Clone(Image32f *img){
      Image32f *imgc;
      imgc = Create(img->ncols,img->nrows);
      memcpy(imgc->data,img->data,img->ncols*img->nrows*sizeof(float));
      imgc->dx = img->dx;
      imgc->dy = img->dy;
      return(imgc);
    }


    Image32f *Clone(Image32::Image32 *img){
      Image32f *imgc;
      int i;
      imgc = Create(img->ncols,img->nrows);
      for (i=0; i < img->n; i++)
	imgc->data[i] = (float)img->data[i];
      imgc->dx = img->dx;
      imgc->dy = img->dy;
      return(imgc);
    }

    
    void Set(Image32f *img, float value){
      int i;
      for (i=0; i < img->n; i++)
	img->data[i] = value;
    }
    
    
    bool IsValidPixel(Image32f *img, int x, int y){
      return ((x >= 0)&&(x < img->ncols)&&
	      (y >= 0)&&(y < img->nrows));
    }

    
    Image32f *AddFrame(Image32f *img, int sz, float value){
      Image32f *fimg;
      int y,nbytes,offset;
      float *dst,*src;
      
      fimg = Create(img->ncols+(2*sz), img->nrows+(2*sz));
      Set(fimg, value);
      nbytes = sizeof(float)*img->ncols;
      offset = sz + fimg->ncols*sz;
      for (y=0,src=img->data,dst=fimg->data+offset; y < img->nrows;y++,src+=img->ncols,dst+=fimg->ncols){
	memcpy(dst, src, nbytes);
      }
      return(fimg);
    }


    Image32f *RemFrame(Image32f *fimg, int sz){
      Image32f *img;
      int y,nbytes,offset;
      float *dst,*src;
      
      img = Create(fimg->ncols-(2*sz), fimg->nrows-(2*sz));
      nbytes = sizeof(float)*img->ncols;
      offset = sz + fimg->ncols*sz;
      for (y=0,src=fimg->data+offset,dst=img->data; y < img->nrows;y++,src+=fimg->ncols,dst+=img->ncols){
	memcpy(dst, src, nbytes);
      }
      return(img);
    }

  
    

  } /*end Image32f namespace*/
} /*end gft namespace*/

