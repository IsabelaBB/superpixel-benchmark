#ifndef _GFT_IMAGE32F_H_
#define _GFT_IMAGE32F_H_

#include "gft_common.h"
#include "gft_image32.h"

namespace gft{

  namespace Image32f{

    /**
     * It supports both linear and two-dimensional access 
     * (i.e., img->data[p] or img->array[y][x] for a pixel
     * (x,y) at address p=x+y*xsize).
     */
    typedef struct _image32f {
      float *data;
      float **array;
      int nrows; /* numero de linhas (altura) */
      int ncols; /* numero de colunas (largura) */
      int n;     /* numero de pixels */
      float dx;
      float dy;
    } Image32f;

    /**
     * \brief A constructor.
     */
    Image32f *Create(int ncols,int nrows);
    Image32f *Create(Image32f *img);
    
    /**
     * \brief A destructor.
     */
    void    Destroy(Image32f **img);

    /**
     * \brief A copy constructor.
     */
    Image32f *Clone(Image32f *img);
    Image32f *Clone(Image32::Image32 *img);

    void    Set(Image32f *img, float value);

    bool    IsValidPixel(Image32f *img, int x, int y);


    Image32f *AddFrame(Image32f *img, int sz, float value);
    Image32f *RemFrame(Image32f *fimg, int sz);

    
  } //end Image32f namespace
} //end gft namespace



#endif

