
#ifndef _GFT_RADIOMETRIC3_H_
#define _GFT_RADIOMETRIC3_H_

#include "gft_common.h"
#include "gft_scene.h"
#include "gft_curve.h"

namespace gft{


  namespace Scene32{
    Curve::Curve *Histogram(Scene32 *scn, int binwidth);

    void          LinearStretchinplace(Scene32 *scn, 
				       int omin,int omax,
				       int nmin,int nmax);
    Scene32      *LinearStretch(Scene32 *scn, 
				int omin,int omax,
				int nmin,int nmax);
  } //end Scene32 namespace


  namespace Scene16{
    Curve::Curve *Histogram(Scene16 *scn, int binwidth);

    void          LinearStretchinplace(Scene16 *scn, 
				       int omin,int omax,
				       int nmin,int nmax);
    Scene16      *LinearStretch(Scene16 *scn, 
				int omin,int omax,
				int nmin,int nmax);
  } //end Scene16 namespace


  namespace Scene8{
    Curve::Curve *Histogram(Scene8 *scn, int binwidth);

    void          LinearStretchinplace(Scene8 *scn, 
				       int omin,int omax,
				       int nmin,int nmax);
    Scene8       *LinearStretch(Scene8 *scn, 
				int omin,int omax,
				int nmin,int nmax);
  } //end Scene8 namespace


  namespace Scene{
    Curve::Curve *Histogram(Scene *scn, int binwidth);

    void          LinearStretchinplace(Scene *scn, 
				       int omin,int omax,
				       int nmin,int nmax);
    Scene        *LinearStretch(Scene *scn, 
				int omin,int omax,
				int nmin,int nmax);
  } //end Scene namespace


} //end gft namespace

#endif

