
#ifndef _GFT_CURVE_H_
#define _GFT_CURVE_H_

#include "gft_common.h"

namespace gft{
  namespace Curve{

  typedef struct _curve {
    float *X;
    float *Y;
    int n;
  } Curve;


  Curve *Create(int n);
  void   Destroy(Curve **curve);
  Curve *Clone(Curve *curve);

  /**
   * Compatible with Gnuplot.
   */
  Curve *Read(char *filename);

  /**
   * Compatible with Gnuplot.
   */
  void   Write(Curve *curve,char *filename);

  Curve *Normalize(Curve *curve);
  void   Normalizeinplace(Curve *curve);

  int    LowerPercentage(Curve *curve,
			 float perc);
  int    HigherPercentage(Curve *curve,
			  float perc);
  int    Median(Curve *curve);
  int    Otsu(Curve *curve);

  void   Sort(Curve *curve);

  void   InvertXY(Curve *curve);
  
  } //end Curve namespace
} //end gft namespace


#endif
