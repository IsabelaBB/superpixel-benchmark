#ifndef _GFT_CIMAGE32F_H_
#define _GFT_CIMAGE32F_H_

#include "gft_image32f.h"
#include "gft_cimage.h"
#include "gft_color.h"

namespace gft{
  namespace CImage32f{

    typedef struct cimage32f {
      Image32f::Image32f *C[3];
    } CImage32f;

    CImage32f *Create(int ncols, int nrows);
    CImage32f *Create(CImage32f *cimg);
    CImage32f *Create(Image32f::Image32f *img);
    void       Destroy(CImage32f **cimg);
    CImage32f *Clone(CImage32f *cimg);
    CImage32f *Clone(CImage::CImage *cimg);

    void    Set(CImage32f *cimg, float r, float g, float b);
    
    CImage32f *RGB2Lab(gft::CImage::CImage *cimg);
 
    CImage32f *AddFrame(CImage32f *cimg, int sz, float r, float g, float b);
    CImage32f *RemFrame(CImage32f *cimg, int sz);
    
  } //end CImage32f namespace
} //end gft namespace

#endif

