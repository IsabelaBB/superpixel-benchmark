
#include "gft_image32.h"

namespace gft{
  namespace Image32{

    Image32 *Create(Image32 *img){
      Image32 *nimg=NULL;
      nimg = Create(img->ncols, img->nrows);
      nimg->dx = img->dx;
      nimg->dy = img->dy;
      return nimg;
    }
    
    Image32 *Create(int ncols, int nrows){
      Image32 *img=NULL;
      int i;
      img = (Image32 *) calloc(1,sizeof(Image32));
      if (img == NULL){
	gft::Error((char *)MSG1,(char *)"Image32::Create");
      }
      img->data  = gft::AllocIntArray(nrows*ncols);
      img->ncols = ncols;
      img->nrows = nrows;
      img->n = nrows*ncols;
      img->dx = img->dy = 1.0;
      
      img->array = (int**)malloc(nrows*sizeof(int*));
      if(img->array == NULL){
	gft::Error((char *)MSG1,(char *)"Image32::Create");
      }
      for(i = 0; i < nrows; i++){
	img->array[i] = (img->data + i*ncols);
      }
      return(img);
    }
    

    void Destroy(Image32 **img){
      Image32 *aux;
      if(img != NULL){
	aux = *img;
	if (aux != NULL){
	  if(aux->data !=  NULL) gft::FreeIntArray(&aux->data);
	  if(aux->array != NULL) free(aux->array);
	  free(aux);
	  *img = NULL;
	}
      }
    }
    
    
    Image32 *Clone(Image32 *img){
      Image32 *imgc;
      imgc = Create(img->ncols,img->nrows);
      memcpy(imgc->data,img->data,img->ncols*img->nrows*sizeof(int));
      imgc->dx = img->dx;
      imgc->dy = img->dy;
      return(imgc);
    }


    Image32 *Clone(Image32 *img, Pixel l, Pixel h){
      int x,y,p,i;
      Image32 *sub=NULL;

      if (l.x > h.x || l.y > h.y)
	gft::Error((char *)"Invalid coordinates",(char *)"Clone");
      
      sub = Create(h.x-l.x+1, h.y-l.y+1);
      i=0;
      for (y=l.y; y<=h.y; y++){
	for (x=l.x; x<=h.x; x++){
	  if (IsValidPixel(img, x, y)){
	    p = x + y*img->ncols;
	    sub->data[i] = img->data[p];
	  }
	  i++;
	}
      }
      return(sub);
    }


    void    Copy(Image32 *img,
		 Image32 *sub,
		 Pixel l){
      int i,j,p,q;
      q = 0;
      for (i=0; i<sub->nrows; i++){
        for (j=0; j<sub->ncols; j++){
	  if (IsValidPixel(img,l.x+j,i+l.y)){
	    p = (l.x+j) + (i+l.y)*img->ncols;
	    img->data[p] = sub->data[q];
	  }
	  q++;
        }
      }
    }

    
    Image32 *Add( Image32 *img1,  Image32 *img2){
      Image32 *sum;
      int p,n;
      sum = Clone(img1);
      n = img1->n;
      for(p = 0; p < n; p++)
	sum->data[p] += img2->data[p];
      return(sum);
    }

    
    Image32 *Add( Image32 *img,   int value){
      Image32 *sum;
      int p,n;
      sum = Clone(img);
      n = img->n;
      for(p = 0; p < n; p++)
	sum->data[p] += value;
      return(sum);
    }


    Image32 *Mult(Image32 *img1, Image32 *img2){
      Image32 *mult;
      int p,n;
      mult = Clone(img1);
      n = img1->n;
      for(p = 0; p < n; p++)
	mult->data[p] *= img2->data[p];
      return(mult);
    }
   
    
    Image32 *Read(char *filename){
      FILE *fp=NULL;
      unsigned char *value=NULL;
      char type[10];
      int  i,ncols,nrows,n;
      Image32 *img=NULL;
      char z[256];
      
      fp = fopen(filename,"rb");
      if (fp == NULL){
        fprintf(stderr,"Cannot open %s\n",filename);
        exit(-1);
      }
      fscanf(fp,"%s\n",type);
      if ((strcmp(type,"P5")==0)){
	gft::NCFgets(z,255,fp);
        sscanf(z,"%d %d\n",&ncols,&nrows);
        n = ncols*nrows;
	gft::NCFgets(z,255,fp);
        sscanf(z,"%d\n",&i);

	fgetc(fp);
	//fseek(fp, -n, SEEK_END);

        value = (unsigned char *)calloc(n,sizeof(unsigned char));
        if (value != NULL){
	  fread(value,sizeof(unsigned char),n,fp);
        }
        else{
	  gft::Error((char *)MSG1,(char *)"Image32::Read");
        }
        fclose(fp);
        img = Create(ncols,nrows);
        for (i=0; i < n; i++)
	  img->data[i] = (int)value[i];
        free(value);
      }
      else{
        if ((strcmp(type,"P2")==0)){
	  gft::NCFgets(z,255,fp);
	  sscanf(z,"%d %d\n",&ncols,&nrows);
	  n = ncols*nrows;
	  gft::NCFgets(z,255,fp);
	  sscanf(z,"%d\n",&i);
	  img = Create(ncols,nrows);
	  for (i=0; i < n; i++)
	    fscanf(fp,"%d",&img->data[i]);
	  fclose(fp);
        }
        else{
	  fprintf(stderr,"Input image must be P2 or P5\n");
	  exit(-1);
        }
      }
      return(img);
    }
    
    
    void Write(Image32 *img,char *filename){
      FILE *fp;
      int i, n, Imax;
      
      fp = fopen(filename,"wb");
      if (fp == NULL){
        fprintf(stderr,"Cannot open %s\n",filename);
        exit(-1);
      }
      n    = img->ncols*img->nrows;
      if ((Imax=GetMaxVal(img))==INT_MAX){
	gft::Warning((char *)"Image with infinity values",(char *)"Image::Write");
	Imax = INT_MIN;
	for (i=0; i < n; i++)
	  if ((img->data[i] > Imax)&&(img->data[i]!=INT_MAX))
	    Imax = img->data[i];
        fprintf(fp,"P2\n");
        fprintf(fp,"%d %d\n",img->ncols,img->nrows);
        fprintf(fp,"%d\n",Imax+1);
      }
      else{
        fprintf(fp,"P2\n");
        fprintf(fp,"%d %d\n",img->ncols,img->nrows);
        if (Imax==0) Imax++;
        fprintf(fp,"%d\n",Imax);
      }
      
      for (i=0; i < n; i++){
	if (img->data[i]==INT_MAX)
	  fprintf(fp,"%d ",Imax+1);
        else
	  fprintf(fp,"%d ",img->data[i]);
	if (((i+1)%17) == 0)
	  fprintf(fp,"\n");
      }
      fclose(fp);
    }
    
    
    Image32 *ConvertToNbits(Image32 *img, int N){
      Image32 *imgN;
      int min,max,i,n,Imax;
      
      imgN = Create(img->ncols,img->nrows);
      n    = img->ncols*img->nrows;
      Imax = (int)(pow(2,N)-1);
      
      min = INT_MAX;
      max = INT_MIN;
      for (i=0; i < n; i++){
        if ((img->data[i] != INT_MIN)&&(img->data[i] != INT_MAX)){
	  if (img->data[i] > max)
	    max = img->data[i];
	  if (img->data[i] < min)
	    min = img->data[i];
        }
      }
      
      if (min != max)
        for (i=0; i < n; i++){
	  if ((img->data[i] != INT_MIN)&&(img->data[i] != INT_MAX)){
	    imgN->data[i] = (int)(((float)Imax*(float)(img->data[i] - min))/
				  (float)(max-min));
	  }
	  else{
	    if (img->data[i]==INT_MIN)
	      imgN->data[i] = 0;
	    else
	      imgN->data[i] = Imax;
	  }
        }
      return(imgN);
    }


    Image32 *Complement(Image32 *img){
      Image32 *comp;
      int p,Imax;
      Imax = GetMaxVal(img);
      comp = Create(img);
      for(p = 0; p < img->n; p++)
	comp->data[p] = Imax - img->data[p];
      return comp;
    } 

    
    int GetMinVal(Image32 *img){
      int i,min,n;
      n = img->ncols*img->nrows;
      min = img->data[0];
      for (i=1; i < n; i++)
        if (img->data[i] < min)
	  min = img->data[i];
      
      return(min);
    }
    
    int GetMaxVal(Image32 *img){
      int i,max,n;
      n = img->ncols*img->nrows;
      max = img->data[0];
      for (i=1; i < n; i++)
        if (img->data[i] > max)
	  max = img->data[i];
      
      return(max);
    }

    int GetMaxVal(Image32 *img, int ignoredvalue){
      int p,n,Imax=INT_MIN;
      n = img->ncols*img->nrows;
      for(p=0; p<n; p++){
	if(img->data[p] != ignoredvalue && img->data[p]>Imax)
	  Imax = img->data[p];
      }
      return Imax;
    }

    
    void Set(Image32 *img, int value){
      int i;
      if(value == 0)
	bzero(img->data, img->n*sizeof(int));
      else
	for (i=0; i < img->n; i++)
	  img->data[i] = value;
    }
    
    
    bool IsValidPixel(Image32 *img, int x, int y){
      return ((x >= 0)&&(x < img->ncols)&&
	      (y >= 0)&&(y < img->nrows));
    }
    

    Image32 *Threshold(Image32 *img, int L, int H){
      Image32 *bin = Clone(img);
      int p;
      for(p = 0; p<img->ncols*img->nrows; p++){
	if(img->data[p] >= L && img->data[p] <= H){
	  bin->data[p] = 1;
	}
	else
	  bin->data[p] = 0;
      }
      return bin;
    }
    
    
    void DrawRectangle(Image32 *img, 
		       int x1, int y1, 
		       int x2, int y2, int val){
      int i,j;
      for(i=y1; i<=y2; i++){
	for(j=x1; j<=x2; j++){
	  img->array[i][j] = val;
	}
      }
    }


    void DrawLineDDA(Image32 *img, 
		     int x1, int y1, 
		     int xn, int yn, int val){
      int vx, vy;
      float Dx, Dy;
      int amostras; /* numero de pontos a serem pintados */
      float m; /* coeficiente angular da reta */
      int i;
      float xk, yk;
      int p;
      
      vx = xn - x1;
      vy = yn - y1;
      
      if (vx == 0){
	Dx = 0.0;
	Dy = (float) SIGN(vy);
	amostras = abs(vy)+1;
      }
      else{
	m = ((float)vy )/((float)vx);
	if ( abs(vx) > abs(vy)){
	  Dx = (float) SIGN(vx);
	  Dy = m * Dx;
	  amostras = abs(vx)+1;
	}
	else{
	  Dy = (float) SIGN(vy);
	  Dx = Dy / m;
	  amostras = abs(vy)+1;
	}
      }
      
      xk = (float) x1;
      yk = (float) y1;
      for (i = 0; i < amostras; i++){
	if ( IsValidPixel(img, ROUND(xk), ROUND(yk)) ){
	  p = ROUND(xk)+img->ncols*(ROUND(yk));
	  img->data[p] = val;
	}
	xk += Dx;
	yk += Dy;
      }
    }
    

    void    DrawCircle(Image32 *img,
		       int x1, int y1,
		       float r,
		       int val){
      int i,j;
      for(i=y1 - ROUND(r+1); i<=y1 + ROUND(r+1); i++){
	for(j=x1 - ROUND(r+1); j<=x1 + ROUND(r+1); j++){
	  if(IsValidPixel(img, j, i)){
	    if((i-y1)*(i-y1)+(j-x1)*(j-x1) <= r*r)
	      img->array[i][j] = val;
	  }
	}
      }
    }


    //------------------------------------
    
    Image32 *AddFrame(Image32 *img, int sz, int value){
      Image32 *fimg;
      int y,*dst,*src,nbytes,offset;
      
      fimg = Create(img->ncols+(2*sz), img->nrows+(2*sz));
      Set(fimg, value);
      nbytes = sizeof(int)*img->ncols;
      offset = sz + fimg->ncols*sz;
      for (y=0,src=img->data,dst=fimg->data+offset; y < img->nrows;y++,src+=img->ncols,dst+=fimg->ncols){
	memcpy(dst, src, nbytes);
      }
      return(fimg);
    }


    Image32 *RemFrame(Image32 *fimg, int sz){
      Image32 *img;
      int y,*dst,*src,nbytes,offset;
      
      img = Create(fimg->ncols-(2*sz), fimg->nrows-(2*sz));
      nbytes = sizeof(int)*img->ncols;
      offset = sz + fimg->ncols*sz;
      for (y=0,src=fimg->data+offset,dst=img->data; y < img->nrows;y++,src+=fimg->ncols,dst+=img->ncols){
	memcpy(dst, src, nbytes);
      }
      return(img);
    }

    
    Image32 *Scale_linear(Image32 *img, float Sx, float Sy){
      float S[2][2],x,y,d1,d2,d3,d4,Ix1,Ix2,If;
      Image32 *scl;
      Pixel u,v,prev,next;
      
      if (Sx == 0.0) Sx = 1.0;
      if (Sy == 0.0) Sy = 1.0;
      
      S[0][0] = 1.0/Sx;
      S[0][1] = 0;
      S[1][0] = 0;
      S[1][1] = 1.0/Sy;
      
      scl = Create((int)(img->ncols*fabs(Sx) + 0.5),
		   (int)(img->nrows*fabs(Sy) + 0.5));
      
      for (v.y=0; v.y < scl->nrows; v.y++)
        for (v.x=0; v.x < scl->ncols; v.x++){
	  x = ((v.x-scl->ncols/2.)*S[0][0]) + img->ncols/2.;
	  y = ((v.y-scl->nrows/2.)*S[1][1]) + img->nrows/2.;
	  u.x = (int)(x+0.5);
	  u.y = (int)(y+0.5);
	  if (IsValidPixel(img,u.x,u.y)){
	    if (x < u.x){
	      next.x = u.x;
	      prev.x = u.x - 1;
	    }
	    else{
	      next.x = u.x + 1;
	      prev.x = u.x;
	    }
	    d1 = next.x - x;
	    d2 = x - prev.x;
	    if (y < u.y){
	      next.y = u.y;
	      prev.y = u.y - 1;
	    }
	    else{
	      next.y = u.y + 1;
	      prev.y = u.y;
	    }
	    d3 = next.y - y;
	    d4 = y - prev.y;
	    
	    if (IsValidPixel(img,prev.x,prev.y)&&
		IsValidPixel(img,next.x,prev.y))
	      Ix1 = (d1*img->array[prev.y][prev.x] +
		     d2*img->array[prev.y][next.x]);
	    else
	      Ix1 = img->array[u.y][u.x];
	    
	    if (IsValidPixel(img,prev.x,next.y)&&
		IsValidPixel(img,next.x,next.y))
	      Ix2 = (d1*img->array[next.y][prev.x] +
		     d2*img->array[next.y][next.x]);
	    else
	      Ix2 = img->array[u.y][u.x];
	    
	    If = d3*Ix1 + d4*Ix2;
	    
	    scl->array[v.y][v.x] = (int)If;
	  }
        }
      
      return(scl);
    }

    
    Image32 *Scale_none(Image32 *img, float Sx, float Sy){
      float S[2][2],x,y;
      Image32 *scl;
      Pixel u,v;
      
      if (Sx == 0.0) Sx = 1.0;
      if (Sy == 0.0) Sy = 1.0;
      
      S[0][0] = 1.0/Sx;
      S[0][1] = 0;
      S[1][0] = 0;
      S[1][1] = 1.0/Sy;
      
      scl = Create((int)(img->ncols*Sx),
		   (int)(img->nrows*Sy));
      
      for (v.y=0; v.y < scl->nrows; v.y++)
        for (v.x=0; v.x < scl->ncols; v.x++){
	  x = ((v.x-scl->ncols/2.)*S[0][0]) + img->ncols/2.;
	  y = ((v.y-scl->nrows/2.)*S[1][1]) + img->nrows/2.;
	  u.x = (int)(x);
	  u.y = (int)(y);
	  if (IsValidPixel(img,u.x,u.y)){
	    scl->array[v.y][v.x] = img->array[u.y][u.x];
	  }
        }
      return(scl);
    }
    
    
    Image32 *Scale(Image32 *img, float Sx, float Sy,
		   InterpolationType interpolation){
      if(interpolation == none)
	return Scale_none(img, Sx, Sy);
      else if(interpolation == linear)
	return Scale_linear(img, Sx, Sy);
      else
	return NULL;
    }


    

  } /*end Image32 namespace*/
} /*end gft namespace*/

