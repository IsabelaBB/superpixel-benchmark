
#ifndef _GFT_FILTERING_H_
#define _GFT_FILTERING_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_adjrel.h"

namespace gft{

  namespace Kernel{

    typedef struct _kernel {
      float *val;
      AdjRel::AdjRel *adj;
    } Kernel;

    Kernel *Make(char *coefs);
    Kernel *Create(AdjRel::AdjRel *A);
    Kernel *Clone(Kernel *K);
    Kernel *Normalize(Kernel *K);    
    void    Destroy(Kernel **K);

    Kernel *Gaussian(AdjRel::AdjRel *A, float stddev);
    
  } //end Kernel namespace

  namespace Image32{

    Image32 *GaussianBlur(Image32 *img, float stddev);
    Image32 *GaussianBlur(Image32 *img);
    
    //void ModeFilterLabel(Image32 *label, float r);
    Image32 *SobelFilter(Image32 *img);
    Image32 *LinearFilter(Image32 *img, Kernel::Kernel *K);

    Image32 *ImageMagnitude(Image32 *imgx, Image32 *imgy);

  } //end Image32 namespace


} //end gft namespace


#endif


