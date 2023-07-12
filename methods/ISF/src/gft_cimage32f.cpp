
#include "gft_cimage32f.h"

namespace gft{
  namespace CImage32f{


    CImage32f *Create(int ncols, int nrows){
      CImage32f *cimg=NULL;
      int i;
     
      cimg = (CImage32f *) calloc(1, sizeof(CImage32f));
      for (i=0; i < 3; i++)
	cimg->C[i] = Image32f::Create(ncols,nrows);
      return(cimg);
    }


    CImage32f *Create(CImage32f *cimg){
      return Create(cimg->C[0]->ncols, cimg->C[0]->nrows);
    }


    CImage32f *Create(Image32f::Image32f *img){
      return Create(img->ncols, img->nrows);
    }


    void    Destroy(CImage32f **cimg){
      CImage32f *tmp;
      int i;
      
      tmp = *cimg;
      if (tmp != NULL) {
	for (i=0; i < 3; i++)
	  Image32f::Destroy(&(tmp->C[i]));
	free(tmp);
	*cimg = NULL;
      }
    }

    
    CImage32f *Clone(CImage32f *cimg){
      CImage32f *imgc;
      int i;
      
      imgc = (CImage32f *) calloc(1,sizeof(CImage32f));
      if (imgc == NULL){
	gft::Error((char *)MSG1,(char *)"CImage32f::Clone");
      }
      for (i=0; i<3; i++)
	imgc->C[i] = Image32f::Clone(cimg->C[i]);
      return imgc;
    }


    CImage32f *Clone(CImage::CImage *cimg){
      CImage32f *imgc;
      int i;
      
      imgc = (CImage32f *) calloc(1,sizeof(CImage32f));
      if (imgc == NULL){
	gft::Error((char *)MSG1,(char *)"CImage32f::Clone");
      }
      for (i=0; i<3; i++)
	imgc->C[i] = Image32f::Clone(cimg->C[i]);
      return imgc;
    }


    void    Set(CImage32f *cimg, float r, float g, float b){
      Image32f::Set(cimg->C[0], r);
      Image32f::Set(cimg->C[1], g);
      Image32f::Set(cimg->C[2], b);
    }



    CImage32f *RGB2Lab(gft::CImage::CImage *cimg){
      CImage32f *cimg_lab;
      double l,a,b;
      int p,n;
      n = cimg->C[0]->n;
      cimg_lab = Create(cimg->C[0]->ncols, cimg->C[0]->nrows);

      //pragma omp parallel for private(l,a,b)
      for(p = 0; p < n; p++){
	gft::Color::RGB2Lab(cimg->C[0]->data[p],
			    cimg->C[1]->data[p],
			    cimg->C[2]->data[p],
			    l, a, b);
	cimg_lab->C[0]->data[p] = l;
	cimg_lab->C[1]->data[p] = a;
	cimg_lab->C[2]->data[p] = b;
      }
      return cimg_lab;
    }



    CImage32f *AddFrame(CImage32f *cimg, int sz, float r, float g, float b){
      CImage32f *fcimg;
      fcimg = (CImage32f *) calloc(1,sizeof(CImage32f));
      if (fcimg == NULL)
	gft::Error((char *)MSG1,(char *)"CImage32f::AddFrame");
      fcimg->C[0] = gft::Image32f::AddFrame(cimg->C[0], sz, r);
      fcimg->C[1] = gft::Image32f::AddFrame(cimg->C[1], sz, g);
      fcimg->C[2] = gft::Image32f::AddFrame(cimg->C[2], sz, b);
      return fcimg;
    }

    
    CImage32f *RemFrame(CImage32f *cimg, int sz){
      CImage32f *fcimg;
      fcimg = (CImage32f *) calloc(1,sizeof(CImage32f));
      if (fcimg == NULL)
	gft::Error((char *)MSG1,(char *)"CImage32f::AddFrame");
      fcimg->C[0] = gft::Image32f::RemFrame(cimg->C[0], sz);
      fcimg->C[1] = gft::Image32f::RemFrame(cimg->C[1], sz);
      fcimg->C[2] = gft::Image32f::RemFrame(cimg->C[2], sz);
      return fcimg;
    }
  
    
  } /*end CImage32f namespace*/
} /*end gft namespace*/
