
#ifndef _GFT_EVALUATION_H_
#define _GFT_EVALUATION_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_scene32.h"
#include "gft_adjrel.h"
#include "gft_analysis.h"

namespace gft{
  namespace Image32{

    float    DiceSimilarity(Image32 *mask1,
			    Image32 *mask2);
    float    JaccardSimilarity(Image32 *mask1,
			       Image32 *mask2);

    //mask1: Ground Truth
    //mask2: Segmentation Result
    int     AssessTP(Image32 *mask1, Image32 *mask2);
    int     AssessFN(Image32 *mask1, Image32 *mask2);
    int     AssessFP(Image32 *mask1, Image32 *mask2);
    int     AssessTN(Image32 *mask1, Image32 *mask2);
    
    Image32 *GetObjError(Image32 *gtruth,
			 Image32 *mask);
    Image32 *GetBkgError(Image32 *gtruth,
			 Image32 *mask);

    //mask1: Ground Truth
    //mask2: Segmentation Result
    float    BoundaryError(Image32 *mask1,
			   Image32 *mask2);
    float    BoundaryFP(Image32 *mask1,
			Image32 *mask2);
    float    BoundaryFN(Image32 *mask1,
			Image32 *mask2);
    
  } //end Image32 namespace

  namespace Scene32{

    float    DiceSimilarity(Scene32 *mask1,
			    Scene32 *mask2);
    
  } //end Scene32 namespace

} //end gft namespace


#endif
