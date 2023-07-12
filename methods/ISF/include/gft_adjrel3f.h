
#ifndef _ADJREL3F_H_
#define _ADJREL3F_H_

#include "gft_common.h"

namespace gft{
  namespace AdjRel3f{

    typedef struct _adjrel3f {
      float *dx;
      float *dy;
      float *dz;
      int n;
    } AdjRel3f;
    
    AdjRel3f *Create(int n);
    void      Destroy(AdjRel3f **A);
    AdjRel3f *Clone(AdjRel3f *A);
    
    AdjRel3f *ChangeOrientationToLPS(AdjRel3f *A,
				     char *ori);

  } /*end AdjRel3f namespace*/
} /*end gft namespace*/

#endif

