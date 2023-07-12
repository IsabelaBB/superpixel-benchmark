#ifndef _GFT_CIMAGE_H_
#define _GFT_CIMAGE_H_

#include "gft_image32.h"
#include "gft_color.h"

namespace gft{
  namespace CImage{

    typedef struct cimage {
      Image32::Image32 *C[3];
    } CImage;

    CImage *Create(int ncols, int nrows);
    CImage *Create(CImage *cimg);
    CImage *Create(Image32::Image32 *img);
    void    Destroy(CImage **cimg);
    CImage *Clone(CImage *cimg);
    CImage *Clone(CImage *cimg, Pixel l, Pixel h);
    CImage *Clone(Image32::Image32 *img);
    CImage *Clone(Image32::Image32 *img, int Imin, int Imax);
    void    Copy(CImage *cimg, CImage *sub, Pixel l);

    CImage *RandomColorize(gft::Image32::Image32 *img);
    
    CImage *Read(char *filename);
    void    Write(CImage *cimg, char *filename);
    void    Set(CImage *cimg, int r, int g, int b);
    void    Set(CImage *cimg, int color);
    
    CImage *ColorizeLabel(Image32::Image32 *label);

    CImage *RGB2Lab(CImage *cimg);
 
    Image32::Image32 *Lightness(CImage *cimg);

    /*The luminosity method works best overall and is the 
      default method used if you ask GIMP to change an image 
      from RGB to grayscale */
    Image32::Image32 *Luminosity(CImage *cimg);

    void MBB(CImage *cimg, int bkgcolor, Pixel *l, Pixel *h);

    CImage *AddFrame(CImage *cimg, int sz, int r, int g, int b);
    CImage *RemFrame(CImage *cimg, int sz);
    
  } //end CImage namespace
} //end gft namespace

#endif

