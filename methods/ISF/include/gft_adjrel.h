#ifndef _GFT_ADJREL_H_
#define _GFT_ADJREL_H_

#include "gft_common.h"

namespace gft{
  namespace AdjRel{

    typedef struct _adjrel {
      int *dx;
      int *dy;
      int n;
    } AdjRel;

    typedef struct _adjpxl {
      int *dp;
      int n;
    } AdjPxl;
    
    
    AdjRel *Create(int n);
    void    Destroy(AdjRel **A);
    AdjRel *Clone(AdjRel *A);
    
    AdjRel *Neighborhood_4(); /* 4-neighborhood */
    AdjRel *Neighborhood_8(); /* 8-neighborhood */
    AdjRel *Neighborhood_8_counterclockwise();
    AdjRel *Neighborhood_8_clockwise();
   
    AdjRel *Circular(float r);
    AdjRel *Box(int ncols, int nrows);

    //-----------------------------------
    int GetFrameSize(AdjRel *A);

    void Mult(AdjRel *A, int val);
    
  } /*end AdjRel namespace*/

} /*end gft namespace*/


#include "gft_image32.h"

namespace gft{

  namespace Image32{

    Image32 *Render(gft::AdjRel::AdjRel *A);

    void DrawAdjRel(Image32 *img,
		    gft::AdjRel::AdjRel *A, 
		    int p, int val);

    Image32 *GetBoundaries(Image32 *img, AdjRel::AdjRel *A);

    
  } /*end Image32 namespace*/
 
} /*end gft namespace*/




#include "gft_cimage.h"

namespace gft{
  namespace CImage{

    void DrawAdjRel(CImage *cimg,
		    gft::AdjRel::AdjRel *A, 
		    int p, int color);

    
  } //end CImage namespace
} //end gft namespace



namespace gft{
  namespace AdjRel{

    AdjPxl  *AdjPixels(AdjRel *A, int ncols);
    AdjPxl  *AdjPixels(AdjRel *A, gft::Image32::Image32 *img);

    AdjPxl  *AdjPixels(AdjRel *A, gft::CImage::CImage *cimg);

    void DestroyAdjPxl(AdjPxl **N);
    
  } /*end AdjRel namespace*/
} /*end gft namespace*/


#endif

