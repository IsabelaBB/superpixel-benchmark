
#ifndef _GFT_COLOR_H_
#define _GFT_COLOR_H_

#include "gft_common.h"

namespace gft{
  namespace Color{

    /* color-related functions */
    int Triplet(int a,int b,int c);

    /*    
    int RGB2HSV(int vi);
    int HSV2RGB(int vi);
    */

    typedef struct _ColorRGB{
      int r;
      int g;
      int b;
    } ColorRGB;
    
    typedef struct _ColorHSV{
      int h;
      int s;
      int v;
    } ColorHSV;

    /*
    typedef struct _ColorLab{
      float l;
      float a;
      float b;
    } ColorLab;
    */
    
    ColorHSV RGB2HSV(ColorRGB rgb);
    ColorRGB HSV2RGB(ColorHSV hsv);

    void RGB2Lab(const int& R, const int& G, const int& B,
		 double &l, double &a, double &b);
    void RGB2Lab(const int& R, const int& G, const int& B,
		 float &l, float &a, float &b);
    //ColorLab RGB2Lab(ColorRGB rgb);

    //ColorRGB Lab2RGB(ColorLab lab);

    int RandomColor();
    
    inline int Channel0(int c);
    inline int Channel1(int c);
    inline int Channel2(int c);

    /* merges two colors (1-ratio) of a and ratio of b, in RGB-space */
    int MergeRGB(int a,int b,float ratio);
    
    //---------inline definitions------------------

    inline int Channel0(int c){
      return ((c>>16)&0xff);
    }

    inline int Channel1(int c){
      return ((c>>8)&0xff);
    }

    inline int Channel2(int c){
      return (c&0xff);
    }

  } /*end Color namespace*/
} /*end gft namespace*/

#endif

