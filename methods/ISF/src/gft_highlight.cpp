
#include "gft_highlight.h"

namespace gft{
  namespace Highlight{


    CImage::CImage *CWide(CImage::CImage *cimg, 
			  Image32::Image32 *label, 
			  float radius, int color, bool fill){
      CImage::CImage *hcimg;
      int p,q,i;
      AdjRel::AdjRel *A=NULL;
      int u_x,u_y,v_x,v_y;
      
      hcimg = CImage::Clone(cimg);
      A    = AdjRel::Circular(radius);
      for (u_y=0; u_y < hcimg->C[0]->nrows; u_y++){
	for (u_x=0; u_x < hcimg->C[0]->ncols; u_x++){
	  p = u_x + u_y*hcimg->C[0]->ncols;
	  if(fill){
	    if(label->data[p] > 0){
	      hcimg->C[0]->data[p] = ROUND(0.3*Color::Channel0(color) + 0.7*hcimg->C[0]->data[p]);
	      hcimg->C[1]->data[p] = ROUND(0.3*Color::Channel1(color) + 0.7*hcimg->C[1]->data[p]);
	      hcimg->C[2]->data[p] = ROUND(0.3*Color::Channel2(color) + 0.7*hcimg->C[2]->data[p]);
	    }
	  }

	  for (i=1; i < A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if (Image32::IsValidPixel(hcimg->C[0],v_x,v_y)){
	      q = v_x + v_y*hcimg->C[0]->ncols;
	      if (label->data[p] > label->data[q]){
		hcimg->C[0]->data[p] = Color::Channel0(color);
		hcimg->C[1]->data[p] = Color::Channel1(color);
		hcimg->C[2]->data[p] = Color::Channel2(color);
		break;
	      }
	    }
	  }
	}
      }
      AdjRel::Destroy(&A);
      return(hcimg);
    }


    Image32::Image32 *Wide(Image32::Image32 *img, 
			   Image32::Image32 *label, 
			   float radius, int value, bool fill){
      Image32::Image32 *himg;
      int p,q,i;
      AdjRel::AdjRel *A=NULL;
      int u_x,u_y,v_x,v_y;
      
      himg = Image32::Clone(img);
      A    = AdjRel::Circular(radius);
      for (u_y=0; u_y < himg->nrows; u_y++){
	for (u_x=0; u_x < himg->ncols; u_x++){
	  p = u_x + u_y*himg->ncols;
	  if(fill){
	    if(label->data[p] > 0)
	      himg->data[p] = ROUND(0.3*value + 0.7*himg->data[p]);
	  }
	  
	  for (i=1; i < A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if (Image32::IsValidPixel(himg,v_x,v_y)){
	      q = v_x + v_y*himg->ncols;
	      if (label->data[p] > label->data[q]){
		himg->data[p] = value;
		break;
	      }
	    }
	  }
	}
      }
      AdjRel::Destroy(&A);
      return(himg);
    }



    CImage::CImage *CWideLabels(CImage::CImage *cimg, 
				Image32::Image32 *label, 
				float radius, int *colormap, float fill, 
				bool thickborder, bool singlecolor){
      CImage::CImage *hcimg;
      int p,q,i,lb;
      AdjRel::AdjRel *A=NULL;
      int u_x,u_y,v_x,v_y, cm=0,Lmax,l,r,g,b;

      if(colormap == NULL){
	Lmax = Image32::GetMaxVal(label);
	colormap = gft::AllocIntArray(Lmax+1);
	gft::RandomSeed();
	for(l = 0; l <= Lmax; l++){
	  r = gft::RandomInteger(0, 255);
	  g = gft::RandomInteger(0, 255);
	  b = gft::RandomInteger(0, 255);
	  colormap[l] = gft::Color::Triplet(r,g,b);
	}
	cm = 1;
      }
      
      hcimg = CImage::Clone(cimg);
      A    = AdjRel::Circular(radius);
      for (u_y=0; u_y < hcimg->C[0]->nrows; u_y++){
	for (u_x=0; u_x < hcimg->C[0]->ncols; u_x++){
	  p = u_x + u_y*hcimg->C[0]->ncols;
	  if(fill > 0.0){
	    lb = label->data[p];
	    if(lb > NIL){
	      if(singlecolor) lb = 0;
	      if(colormap[lb] != NIL){
		hcimg->C[0]->data[p] = ROUND(fill*Color::Channel0(colormap[lb]) + (1.0-fill)*hcimg->C[0]->data[p]);
		hcimg->C[1]->data[p] = ROUND(fill*Color::Channel1(colormap[lb]) + (1.0-fill)*hcimg->C[1]->data[p]);
		hcimg->C[2]->data[p] = ROUND(fill*Color::Channel2(colormap[lb]) + (1.0-fill)*hcimg->C[2]->data[p]);
	      }
	    }
	  }

	  for (i=1; i < A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if (Image32::IsValidPixel(hcimg->C[0],v_x,v_y)){
	      q = v_x + v_y*hcimg->C[0]->ncols;
	      lb = label->data[p];
	      if (lb > label->data[q] || (thickborder && lb != label->data[q])){
		if(singlecolor) lb = 0;
		if(colormap[lb] != NIL){
		  hcimg->C[0]->data[p] = Color::Channel0(colormap[lb]);
		  hcimg->C[1]->data[p] = Color::Channel1(colormap[lb]);
		  hcimg->C[2]->data[p] = Color::Channel2(colormap[lb]);
		}
		break;
	      }
	    }
	  }
	}
      }
      AdjRel::Destroy(&A);
      if(cm == 1){
	gft::FreeIntArray(&colormap);
      }
      return(hcimg);
    }


    bool    HStripedTexture(int x, int y, int w, int h){
      x %= w;
      y %= h;
      if(y==h/2) return 1;
      else       return 0;
    }

    bool    VStripedTexture(int x, int y, int w, int h){
      x %= w;
      y %= h;
      if(x==w/2) return 1;
      else       return 0;
    }

    bool    BackslashTexture(int x, int y, int w, int h){
      x %= w;
      y %= h;
      if(x==y) return 1;
      else     return 0;
    }

    bool    SlashTexture(int x, int y, int w, int h){
      x %= w;
      y %= h;
      if(y==-x+h-1) return 1;
      else          return 0;
    }
    
    bool    GridTexture(int x, int y, int w, int h){
      x %= w;
      y %= h;
      return(VStripedTexture(x,y,w,h) ||
	     HStripedTexture(x,y,w,h));
    }

    bool    RGridTexture(int x, int y, int w, int h){
      x %= w;
      y %= h;
      return(BackslashTexture(x,y,w,h) ||
	     SlashTexture(x,y,w,h));
    }


    Image32::Image32 *Texture(Image32::Image32 *img, 
			      Image32::Image32 *label, 
			      float radius, int value, bool fill,
			      bool (*texture)(int,int,int,int), int w, int h){
      Image32::Image32 *himg;
      int p,q,i;
      AdjRel::AdjRel *A=NULL;
      int u_x,u_y,v_x,v_y;
      bool b;
      
      himg = Image32::Clone(img);
      A    = AdjRel::Circular(radius);
      for (u_y=0; u_y < himg->nrows; u_y++){
	for (u_x=0; u_x < himg->ncols; u_x++){
	  p = u_x + u_y*himg->ncols;
	  
	  if(label->data[p] > 0){
	    b = (*texture)(u_x,u_y,w,h);
	    if(b)
	      himg->data[p] = value;
	    else if(fill)
	      himg->data[p] = ROUND(0.3*value + 0.7*himg->data[p]);
	  }

	  for (i=1; i < A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if (Image32::IsValidPixel(himg,v_x,v_y)){
	      q = v_x + v_y*himg->ncols;
	      if (label->data[p] > label->data[q]){
		himg->data[p] = value;
		break;
	      }
	    }
	  }
	}
      }
      AdjRel::Destroy(&A);
      return(himg);
    }



    CImage::CImage *CTexture(CImage::CImage *cimg, 
			     Image32::Image32 *label, 
			     float radius, int color, bool fill,
			     bool (*texture)(int,int,int,int), int w, int h){
      CImage::CImage *hcimg;
      int p,q,i;
      AdjRel::AdjRel *A=NULL;
      int u_x,u_y,v_x,v_y;
      bool b;
      
      hcimg = CImage::Clone(cimg);
      A    = AdjRel::Circular(radius);
      for (u_y=0; u_y < hcimg->C[0]->nrows; u_y++){
	for (u_x=0; u_x < hcimg->C[0]->ncols; u_x++){
	  p = u_x + u_y*hcimg->C[0]->ncols;
	  
	  if(label->data[p] > 0){
	    b = (*texture)(u_x,u_y,w,h);
	    if(b){
	      hcimg->C[0]->data[p] = Color::Channel0(color);
	      hcimg->C[1]->data[p] = Color::Channel1(color);
	      hcimg->C[2]->data[p] = Color::Channel2(color);
	    }
	    else if(fill){
	      hcimg->C[0]->data[p] = ROUND(0.3*Color::Channel0(color) + 0.7*hcimg->C[0]->data[p]);
	      hcimg->C[1]->data[p] = ROUND(0.3*Color::Channel1(color) + 0.7*hcimg->C[1]->data[p]);
	      hcimg->C[2]->data[p] = ROUND(0.3*Color::Channel2(color) + 0.7*hcimg->C[2]->data[p]);
	    }
	  }

	  for (i=1; i < A->n; i++){
	    v_x = u_x + A->dx[i];
	    v_y = u_y + A->dy[i];
	    if (Image32::IsValidPixel(hcimg->C[0],v_x,v_y)){
	      q = v_x + v_y*hcimg->C[0]->ncols;
	      if (label->data[p] > label->data[q]){
		hcimg->C[0]->data[p] = Color::Channel0(color);
		hcimg->C[1]->data[p] = Color::Channel1(color);
		hcimg->C[2]->data[p] = Color::Channel2(color);
		break;
	      }
	    }
	  }
	}
      }
      AdjRel::Destroy(&A);
      return(hcimg);
    }



  } /*end Highlight namespace*/
} /*end gft namespace*/

