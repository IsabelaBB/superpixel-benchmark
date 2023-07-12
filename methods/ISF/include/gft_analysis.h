#ifndef _GFT_ANALYSIS_H_
#define _GFT_ANALYSIS_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_adjrel.h"
#include "gft_pqueue32.h"
#include "gft_stack.h"


namespace gft{
  namespace Image32{

    Image32 *LabelBinComp(Image32 *bin, gft::AdjRel::AdjRel *A);

    //Image32  *DistTrans(Image32 *bin, gft::AdjRel::AdjRel *A, char side); 
    //Image32  *SignedDistTrans(Image32 *bin, gft::AdjRel::AdjRel *A, char side);

    Image32 *GetObjBorder(Image32 *img);
    Image32 *GetObjBorders(Image32 *img, gft::AdjRel::AdjRel *A);

    Image32 *LabelContour(Image32 *bin);
    Image32 *LabelContour(Image32 *bin, Image32 *contourid);
    
    Image32 *Mask2EDT(Image32 *bin, gft::AdjRel::AdjRel *A,
		      char side, int limit, char sign);
    
    void Mask2EDT(Image32 *bin, gft::AdjRel::AdjRel *A,
		  char side, int limit, char sign,
		  Image32 *cost, Image32 *root);

    Image32 *Multiscaleskeletons(Image32 *bin);
    
  } //end Image32 namespace
} //end gft namespace

#endif

