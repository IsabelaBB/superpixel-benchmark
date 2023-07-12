
#include "gft_scene32f.h"

namespace gft{
  namespace Scene32f{

  Scene32f *Create(int xsize,int ysize,int zsize){
    Scene32f *scn=NULL;
    float **tmp=NULL;
    int xysize;
    int i,j,p,N;

    scn = (Scene32f *) calloc(1,sizeof(Scene32f));
    if(scn == NULL) 
      gft::Error((char *)MSG1,(char *)"Scene32f::Create");
 
    scn->xsize  = xsize;
    scn->ysize  = ysize;
    scn->zsize  = zsize;
    scn->dx     = 1.0;
    scn->dy     = 1.0;
    scn->dz     = 1.0;
    scn->maxval = 0;
    scn->n      = xsize*ysize*zsize;

    //For SSE optimization we allocate a multiple of 4.
    N = scn->n;
    if(N%4!=0) N += (4-N%4);

    scn->data   = gft::AllocFltArray(N);
    if(scn->data==NULL) 
      gft::Error((char *)MSG1,(char *)"Scene32f::Create");

    scn->array = (int ***) calloc(zsize, sizeof(float **));
    if(scn->array==NULL) 
      gft::Error((char *)MSG1,(char *)"Scene32f::Create");

    tmp = (int **) calloc(zsize*ysize, sizeof(float *));
    if(tmp==NULL) 
      gft::Error((char *)MSG1,(char *)"Scene32f::Create");

    scn->array[0] = tmp;
    for(i=1; i<zsize; i++)
      scn->array[i] = scn->array[i-1] + ysize;

    xysize = xsize*ysize;
    for(i=0; i<zsize; i++){
      for(j=0; j<ysize; j++){
	p = j + i*ysize;
	tmp[p] = scn->data + xysize*i + xsize*j;
      }}
    return(scn);
  }


  Scene32f *Create(Scene32f *scn){
    Scene32f *new_scn=NULL;
    new_scn = Create(scn->xsize,
		     scn->ysize,
		     scn->zsize);
    new_scn->dx = scn->dx;
    new_scn->dy = scn->dy;
    new_scn->dz = scn->dz;
    return new_scn;
  }


  void     Destroy(Scene32f **scn){
    Scene32f *aux;
    
    aux = *scn;
    if(aux != NULL){
      if(aux->data     != NULL) gft::FreeFltArray(&aux->data);
      if(aux->array[0] != NULL) free(aux->array[0]);
      if(aux->array    != NULL) free(aux->array);
      free(aux);    
      *scn = NULL;
    }
  }

  Scene32f *Clone(Scene32f *scn){
    Scene32f *aux = Create(scn->xsize,scn->ysize,scn->zsize);
    aux->dx     = scn->dx;
    aux->dy     = scn->dy;
    aux->dz     = scn->dz;
    aux->maxval = scn->maxval;
    memcpy(aux->data,scn->data,sizeof(float)*aux->n);
    return aux;
  }

  Scene32f *SubScene(Scene32f *scn, Voxel l, Voxel h){
    return SubScene(scn,
		    l.c.x, l.c.y, l.c.z,
		    h.c.x, h.c.y, h.c.z);
  }

  Scene32f *SubScene(Scene32f *scn,
		    int xl, int yl, int zl,
		    int xh, int yh, int zh){
    Scene32f *sub=NULL;
    Voxel v;
    int i,j;
    //v4si *p1,*p2;
  
    if(!IsValidVoxel(scn,xl,yl,zl)||
       !IsValidVoxel(scn,xh,yh,zh)||
       (xl > xh)||(yl>yh)||(zl>zh))
      return NULL;

    sub = Create(xh-xl+1,yh-yl+1,zh-zl+1);
    sub->dx = scn->dx;
    sub->dy = scn->dy;
    sub->dz = scn->dz;
    j = 0;
    for(v.c.z=zl; v.c.z<=zh; v.c.z++){
      for(v.c.y=yl; v.c.y<=yh; v.c.y++){
	/*
	for(v.c.x=xl; v.c.x<=xh-3; v.c.x+=4){
	  i = GetVoxelAddress(scn, v);
	  p1 = (v4si *)(&sub->data[j]);
	  p2 = (v4si *)(&scn->data[i]);
	  *p1 = *p2;
	  j+=4;
	  }
	for(; v.c.x <= xh; v.c.x++){
	  i = GetVoxelAddress(scn, v);
	  sub->data[j] = scn->data[i];
	  j++;
	}
	*/
	v.c.x=xl;
	i = GetVoxelAddress(scn, v);
	memcpy(&sub->data[j],&scn->data[i],sizeof(float)*sub->xsize);
	j += sub->xsize;
      }
    }
    return(sub);
  }


  void     Copy(Scene32f *dest, Scene32f *src){
    if(dest->xsize!=src->xsize ||
       dest->ysize!=src->ysize ||
       dest->zsize!=src->zsize)
      gft::Error((char *)"Incompatible data",(char *)"Scene32f::Copy");
    dest->dx     = src->dx;
    dest->dy     = src->dy;
    dest->dz     = src->dz;
    dest->maxval = src->maxval;
    memcpy(dest->data,src->data,sizeof(float)*src->n);
  }


  void     Copy(Scene32f *dest, Scene32f *src, Voxel v){
    int x1,y1,z1,x2,y2,z2,dx;
    Voxel t,u;
    int p,q;
    x1 = MAX(0, -v.c.x);
    x2 = MIN(src->xsize-1, dest->xsize-1 - v.c.x);
    dx = x2 - x1 + 1;
    if(dx<=0) return;

    y1 = MAX(0, -v.c.y);
    y2 = MIN(src->ysize-1, dest->ysize-1 - v.c.y);

    z1 = MAX(0, -v.c.z);
    z2 = MIN(src->zsize-1, dest->zsize-1 - v.c.z);
    
    for(t.c.z=z1; t.c.z<=z2; t.c.z++){
      for(t.c.y=y1; t.c.y<=y2; t.c.y++){
	t.c.x=x1;
	p = GetVoxelAddress(src, t);

	u.v = v.v + t.v;
	q = GetVoxelAddress(dest, u);

	memcpy(&dest->data[q],&src->data[p],sizeof(float)*dx);
      }
    }
  }


  void     Fill(Scene32f *scn, int value){
    if(value==0)
      memset((void *)scn->data, 0, sizeof(float)*scn->n);
    else{
      int p;
      v4sf v;
      v4sf *ptr;
      ((float *)(&v))[0] = value;
      ((float *)(&v))[1] = value;
      ((float *)(&v))[2] = value;
      ((flaot *)(&v))[3] = value;
      for(p=0; p<scn->n; p+=4){
	ptr = (v4sf *)(&scn->data[p]);
	*ptr = v;
      }
    }
    scn->maxval = value;
  }


  //Scene32 *Read(char *filename){
  //}


  void     Write(Scene32f *scn, char *filename){
    FILE *fp=NULL;
    int i,n;

    // Writing the scn file
    fp = fopen(filename,"wb"); 
    if(fp == NULL) 
      gft::Error((char *)MSG2,(char *)"Scene32f::Write");

    fprintf(fp,"SCN\n");
    fprintf(fp,"%d %d %d\n",scn->xsize,scn->ysize,scn->zsize);
    fprintf(fp,"%f %f %f\n",scn->dx,scn->dy,scn->dz);
  
    n = scn->n;
    fprintf(fp,"%d\n",33);
    fwrite(scn->data,sizeof(float),n,fp);
    fclose(fp);
  }


  float GetValue_trilinear(Scene32f *scn, float x,float y, float z){
    float px,py,pz, i1,i2;
    float dx,dy,dz;
    int nbx=0,nby=0,nbz=0;
    float p1,p2,p3,p4,p5,p6,res;

    if(x>scn->xsize-1||y>scn->ysize-1||z>scn->zsize-1||x<0||y<0||z<0)
      gft::Error((char *)"Out-of-bounds",
		 (char *)"Scene32::GetValue_trilinear");

    px = (int)x;  dx=x-px;
    py = (int)y;  dy=y-py;
    pz = (int)z;  dz=z-pz;

    // If it's not on the border, it has a neighour ahead.
    if(px<scn->xsize-1) nbx=1; 
    if(py<scn->ysize-1) nby=1;
    if(pz<scn->zsize-1) nbz=1;

    // 1st: Interpolate in Z
    i1 = scn->data[GetVoxelAddress(scn,px,py,pz)];
    i2 = scn->data[GetVoxelAddress(scn,px,py,pz+nbz)];
    p1 = i1 + dz*(i2-i1);
    i1 = scn->data[GetVoxelAddress(scn,px+nbx,py,pz)];
    i2 = scn->data[GetVoxelAddress(scn,px+nbx,py,pz+nbz)];
    p2 = i1 + dz*(i2-i1);
    i1 = scn->data[GetVoxelAddress(scn,px,py+nby,pz)];
    i2 = scn->data[GetVoxelAddress(scn,px,py+nby,pz+nbz)];
    p3 = i1 + dz*(i2-i1);
    i1 = scn->data[GetVoxelAddress(scn,px+nbx,py+nby,pz)];
    i2 = scn->data[GetVoxelAddress(scn,px+nbx,py+nby,pz+nbz)];
    p4 = i1 + dz*(i2-i1);
    // 2nd: Interpolate in X
    p5 = p1 + dx*(p2-p1);
    p6 = p3 + dx*(p4-p3);
    // 3rd: Interpolate in Y
    res = p5 + dy*(p6-p5);
    return res;
  }


  // return the nearest voxel.
  float   GetValue_nn(Scene32f *scn, float x, float y, float z){
    Voxel v;
    if(x>scn->xsize-1||y>scn->ysize-1||z>scn->zsize-1||
       x<(float)0||y<(float)0||z<(float)0)
      gft::Error((char *)"Out-of-bounds",
		 (char *)"Scene32::GetValue_nn");

    v.c.x=(int)(x+0.5);
    v.c.y=(int)(y+0.5);
    v.c.z=(int)(z+0.5);
    return scn->array[v.c.z][v.c.y][v.c.x];
  }


  int          GetMaximumValue(Scene32f *scn){
    int p;
    float Imax;
    
    Imax = FLT_MIN;
    for(p=0; p<scn->n; p++) {
      if(scn->data[p] > Imax)
	Imax = scn->data[p];
    }
    scn->maxval = Imax;
    return(Imax); 
  }

  float          GetMinimumValue(Scene32f *scn){
    int p;
    float Imin;
  
    Imin = INT_MAX; 
    for(p=0; p<scn->n; p++){
      if(scn->data[p] < Imin)
	Imin = scn->data[p];
    }    
    return(Imin); 
  }


  void     MBB(Scene32f *scn, Voxel *l, Voxel *h){
    Voxel v;
    
    l->c.x  = scn->xsize-1;
    l->c.y  = scn->ysize-1;
    l->c.z  = scn->zsize-1;
    h->c.x = 0;
    h->c.y = 0;
    h->c.z = 0;
  
    for(v.c.z=0; v.c.z<scn->zsize; v.c.z++)
      for(v.c.y=0; v.c.y<scn->ysize; v.c.y++)
	for(v.c.x=0; v.c.x<scn->xsize; v.c.x++)    
	  if(scn->data[GetVoxelAddress(scn,v)] > 0){
	    if(v.c.x < l->c.x)
	      l->c.x = v.c.x;
	    if(v.c.y < l->c.y)
	      l->c.y = v.c.y;
	    if(v.c.z < l->c.z)
	      l->c.z = v.c.z;
	    if(v.c.x > h->c.x)
	      h->c.x = v.c.x;
	    if(v.c.y > h->c.y)
	      h->c.y = v.c.y;	
	    if(v.c.z > h->c.z)
	      h->c.z = v.c.z;	
	  }
  }


  Scene32f *MBB(Scene32f *scn){
    Voxel lower,higher;
    Scene32f *mbb=NULL;
    MBB(scn, &lower, &higher);
    mbb = SubScene(scn,
		   lower.c.x,  lower.c.y,  lower.c.z,
		   higher.c.x, higher.c.y, higher.c.z);
    return(mbb);
  }


  Scene32f *AddFrame(Scene32f *scn,  int sz, float value){
    Scene32f *fscn;
    int y, z,*dst,*src,nbytes,offset1, offset2;
    
    fscn = Create(scn->xsize+(2*sz),
		  scn->ysize+(2*sz), 
		  scn->zsize+(2*sz));
    fscn->dx = scn->dx;
    fscn->dy = scn->dy;
    fscn->dz = scn->dz;

    Fill(fscn,value);
    nbytes = sizeof(float)*scn->xsize;
  
    offset1 = 0;
    offset2 = GetVoxelAddress(fscn, sz, sz, sz);
    
    for(z=0; z<scn->zsize; z++){
      src = scn->data+offset1;
      dst = fscn->data+offset2;
      for(y=0; y<scn->ysize; y++){
	memcpy(dst,src,nbytes);
	src += scn->xsize;
	dst += fscn->xsize;
      }
      offset1 += scn->xsize*scn->ysize;
      offset2 += fscn->xsize*fscn->ysize;
    }
    return(fscn);
  }


  Scene32f *RemFrame(Scene32f *fscn, int sz){
    Scene32f *scn;
    int y,z,*dst,*src,nbytes,offset;
    
    scn = Create(fscn->xsize-(2*sz),
		 fscn->ysize-(2*sz),
		 fscn->zsize-(2*sz));
    scn->dx = fscn->dx;
    scn->dy = fscn->dy;
    scn->dz = fscn->dz;

    nbytes = sizeof(float)*scn->xsize;  
    offset = GetVoxelAddress(fscn, sz, sz, sz);

    src = fscn->data+offset;
    dst = scn->data;
    for(z=0; z<scn->zsize; z++,src+=2*sz*fscn->xsize) {
      for(y=0; y<scn->ysize; y++,src+=fscn->xsize,dst+=scn->xsize){
	memcpy(dst,src,nbytes);
      }
    }
    return(scn);
  }


  Scene16::Scene16 *ConvertTo16(Scene32f *scn){
    Scene16::Scene16 *scn16;
    int p;

    scn16 = Scene16::Create(scn->xsize,
			    scn->ysize,
			    scn->zsize);
    scn16->dx = scn->dx;
    scn16->dy = scn->dy;
    scn16->dz = scn->dz;
    for(p=0; p<scn->n; p++){
      scn16->data[p] = (ushort)scn->data[p];
    }
    return scn16;
  }

  Scene8::Scene8  *ConvertTo8(Scene32f *scn){
    Scene8::Scene8 *scn8;
    int p;

    scn8 = Scene8::Create(scn->xsize,
			  scn->ysize,
			  scn->zsize);
    scn8->dx = scn->dx;
    scn8->dy = scn->dy;
    scn8->dz = scn->dz;
    for(p=0; p<scn->n; p++){
      scn8->data[p] = (uchar)scn->data[p];
    }
    return scn8;
  }



  } //end Scene32f namespace
} //end gft namespace

