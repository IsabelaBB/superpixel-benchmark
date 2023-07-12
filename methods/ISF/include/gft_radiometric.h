
#ifndef _GFT_RADIOMETRIC_H_
#define _GFT_RADIOMETRIC_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_curve.h"

namespace gft{

  namespace Image32{

    gft::Curve::Curve *Histogram(Image32 *img);
    gft::Curve::Curve *NormHistogram(Image32 *img);
    gft::Curve::Curve *NormalizeHistogram(gft::Curve::Curve *hist);
    gft::Curve::Curve *RemoveEmptyBins(gft::Curve::Curve *hist);
    
    Image32 *LinearStretch(Image32 *img, int omin, int omax, int nmin, int nmax);

    
  } //end Image32 namespace
    

} //end gft namespace

#endif


