
#ifndef _GFT_HIGHLIGHT_H_
#define _GFT_HIGHLIGHT_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_cimage.h"
#include "gft_adjrel.h"
#include "gft_color.h"

namespace gft{
  namespace Highlight{
    
    Image32::Image32 *Wide(Image32::Image32 *img, 
			   Image32::Image32 *label, 
			   float radius, int value, bool fill);
    CImage::CImage *CWide(CImage::CImage *cimg, 
			  Image32::Image32 *label, 
			  float radius, int color, bool fill);
    
    CImage::CImage *CWideLabels(CImage::CImage *cimg, 
				Image32::Image32 *label,
				float radius, int *colormap, float fill,
				bool thickborder, bool singlecolor);

    bool    HStripedTexture(int x, int y, int w, int h);
    bool    VStripedTexture(int x, int y, int w, int h);
    bool    BackslashTexture(int x, int y, int w, int h);
    bool    SlashTexture(int x, int y, int w, int h);
    bool    GridTexture(int x, int y, int w, int h);
    bool    RGridTexture(int x, int y, int w, int h);
    
    Image32::Image32 *Texture(Image32::Image32 *img,
			      Image32::Image32 *label,
			      float radius, int value, bool fill,
			      bool (*texture)(int,int,int,int), 
			      int w, int h);
    CImage::CImage *CTexture(CImage::CImage *cimg,
			     Image32::Image32 *label,
			     float radius, int color, bool fill,
			     bool (*texture)(int,int,int,int), int w, int h);

  } //end Highlight namespace
} //end gft namespace

#endif

