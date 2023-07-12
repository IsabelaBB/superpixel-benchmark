
#ifndef _GFT_MORPHOLOGY_H_
#define _GFT_MORPHOLOGY_H_

#include "gft_image32.h"
#include "gft_adjrel.h"
#include "gft_set.h"
#include "gft_pqueue32.h"


namespace gft{
  namespace Image32{

    Image32 *Dilate(Image32 *img, AdjRel::AdjRel *A);
    Image32 *Erode(Image32 *img, AdjRel::AdjRel *A);

    Image32 *ErodeBin(Image32 *bin, gft::Set::Set **seed, float radius);
    
    Image32 *MorphGrad(Image32 *img, gft::AdjRel::AdjRel *A);

    void SupRec_Watershed(gft::AdjRel::AdjRel *A, 
			  Image32 *I, Image32 *J, 
			  Image32 *L, Image32 *V);

    /*Removes all background connected components from the stack
      of binary images of I whose area (number of pixels) is less than 
      a threshold and outputs a simplified image.*/
    Image32 *AreaClosing(gft::AdjRel::AdjRel *A, Image32 *I, int T);

    Image32 *CloseHoles(Image32 *img);
    
  } //end Image32 namespace
} //end gft namespace

#endif
