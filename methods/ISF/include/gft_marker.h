
#ifndef _GFT_MARKER_H_
#define _GFT_MARKER_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_set.h"
#include "gft_adjrel.h"
#include "gft_analysis.h"
#include "gft_morphology.h"


namespace gft{
  namespace Image32{
    
    float  MaxRadiusByErosion(Image32 *bin);
    float  MaxObjRadiusByErosion(Image32 *bin);
    float  MaxBkgRadiusByErosion(Image32 *bin);
    
    Image32 *ObjMaskByErosion(Image32 *bin, float radius);
    Image32 *BkgMaskByErosion(Image32 *bin, float radius);
    
 
   } //end Image32 namespace
} //end gft namespace

 
#endif
