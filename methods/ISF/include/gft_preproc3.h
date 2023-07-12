
#ifndef _GFT_PREPROC_H_
#define _GFT_PREPROC_H_

#include "gft_common.h"
#include "gft_adjrel3.h"
#include "gft_scene32.h"

namespace gft{
  namespace Kernel3{

    typedef struct _kernel3 {
      float *val;
      AdjRel3::AdjRel3 *adj;
      int xsize,ysize,zsize;
    } Kernel3;

    Kernel3 *Create(AdjRel3::AdjRel3 *A);
    Kernel3 *Clone(Kernel3 *K);
    void     Destroy(Kernel3 **K);

    Kernel3 *Normalize(Kernel3 *K);

    Kernel3 *SphericalGaussian(float R, float s, float f);

  } //end Kernel3 namespace
} //end gft namespace

namespace gft{
  namespace Scene32{

    void   SuppressHighIntensities(Scene32 *scn);

    Scene32 *Convolution(Scene32 *scn, Kernel3::Kernel3 *K);
    Scene32 *OptConvolution(Scene32 *scn, Kernel3::Kernel3 *K);

    Scene32 *GaussianBlur(Scene32 *scn);
    Scene32 *OptGaussianBlur(Scene32 *scn);
    Scene32 *FastGaussianBlur(Scene32 *scn);
    Scene32 *FastOptGaussianBlur(Scene32 *scn);
 
    Scene32 *Subsampling(Scene32 *scn);

    //-----------------

    Scene32  *LaplacianFilter(Scene32 *orig);
    Scene32  *SobelFilter(Scene32 *scn);
    Scene32  *SphericalGradient(Scene32 *scn, float r);

    
  } //end Scene32 namespace
} //end gft namespace

#endif
