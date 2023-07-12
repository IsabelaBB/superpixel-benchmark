
#ifndef _GFT_GRADIENT3_H_
#define _GFT_GRADIENT3_H_

#include "gft_common.h"

namespace gft{
  namespace Gradient3{
    
    typedef struct _Gradient3 {
      Scene32::Scene32 *Gx;
      Scene32::Scene32 *Gy;
      Scene32::Scene32 *Gz;
      
      //Available upon request:
      //--> Must call "ComputeMagnitude".
      Scene32::Scene32 *mag;
    } Gradient3;


    Gradient3 *Create(int xsize,int ysize,int zsize);
    void       Destroy(Gradient3 **grad);
    Gradient3 *RemFrame(Gradient3 *fgrad, int sz);
    Gradient3 *LinearInterpCentr(Gradient3 *grad,
				 float dx, float dy, float dz);
    Gradient3 *ChangeOrientationToLPS(Gradient3 *grad,
				      char *ori);
    
    Gradient3 *Read(char *filename);
    void       Write(Gradient3 *grad, char *filename);
    
    Gradient3 *Spherical(Scene32::Scene32 *scn, float r);
    
    void    ComputeMagnitude(Gradient3 *grad);
    int     MaximumMag(Gradient3 *grad);
    void    Normalize(Gradient3 *grad,
		      int omin,int omax,
		      int nmin,int nmax);
    
    void    PowerEnhancement(Gradient3 *grad);


  } //end Gradient3 namespace
} //end gft namespace


#endif

