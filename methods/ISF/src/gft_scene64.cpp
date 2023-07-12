
#include "gft_scene64.h"

extern "C" {
  #include "gft_bzlib.h"
  #include "nifti1_io.h"
}


namespace gft{
  namespace Scene64{
    
    Scene64 *Create(int xsize,int ysize,int zsize){
      Scene64 *scn=NULL;
      long long **tmp=NULL;
      int xysize;
      int i,j,p,N;
      
      scn = (Scene64 *) calloc(1,sizeof(Scene64));
      if(scn == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene64::Create");
      
      scn->xsize  = xsize;
      scn->ysize  = ysize;
      scn->zsize  = zsize;
      scn->dx     = 1.0;
      scn->dy     = 1.0;
      scn->dz     = 1.0;
      scn->maxval = 0;
      scn->n      = xsize*ysize*zsize;
      scn->nii_hdr = NULL;
      
      N = scn->n;
      
      scn->data   = gft::AllocLongLongArray(N);
      if(scn->data==NULL) 
	gft::Error((char *)MSG1,(char *)"Scene64::Create");
      
      scn->array = (long long ***) calloc(zsize, sizeof(long long **));
      if(scn->array==NULL) 
	gft::Error((char *)MSG1,(char *)"Scene64::Create");
      
      tmp = (long long **) calloc(zsize*ysize, sizeof(long *));
      if(tmp==NULL) 
	gft::Error((char *)MSG1,(char *)"Scene64::Create");
      
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
    
    
    Scene64 *Create(Scene64 *scn){
      Scene64 *new_scn=NULL;
      new_scn = Create(scn->xsize,
		       scn->ysize,
		       scn->zsize);
      new_scn->dx = scn->dx;
      new_scn->dy = scn->dy;
      new_scn->dz = scn->dz;
      return new_scn;
    }
    
    
    void     Destroy(Scene64 **scn){
      Scene64 *aux;
      
      aux = *scn;
      if(aux != NULL){
	if(aux->data     != NULL) free(aux->data);
	if(aux->array[0] != NULL) free(aux->array[0]);
	if(aux->array    != NULL) free(aux->array);
	free(aux);    
	*scn = NULL;
      }
    }
    
    Scene64 *Clone(Scene64 *scn){
      Scene64 *aux = Create(scn->xsize,scn->ysize,scn->zsize);
      aux->dx     = scn->dx;
      aux->dy     = scn->dy;
      aux->dz     = scn->dz;
      aux->maxval = scn->maxval;
      memcpy(aux->data,scn->data,sizeof(long long)*aux->n);
      return aux;
    }

    Scene64 *SubScene(Scene64 *scn, Voxel l, Voxel h){
      return SubScene(scn,
		      l.c.x, l.c.y, l.c.z,
		      h.c.x, h.c.y, h.c.z);
    }
    
    Scene64 *SubScene(Scene64 *scn,
		      int xl, int yl, int zl,
		      int xh, int yh, int zh){
      Scene64 *sub=NULL;
      Voxel v;
      int i,j;
      
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
	  v.c.x=xl;
	  i = GetVoxelAddress(scn, v);
	  memcpy(&sub->data[j],&scn->data[i],sizeof(long long)*sub->xsize);
	  j += sub->xsize;
	}
      }
      return(sub);
    }

    
    void     Copy(Scene64 *dest, Scene64 *src){
      if(dest->xsize!=src->xsize ||
	 dest->ysize!=src->ysize ||
	 dest->zsize!=src->zsize)
	gft::Error((char *)"Incompatible data",(char *)"Scene64::Copy");
      dest->dx     = src->dx;
      dest->dy     = src->dy;
      dest->dz     = src->dz;
      dest->maxval = src->maxval;
      memcpy(dest->data,src->data,sizeof(long long)*src->n);
    }


    void     Copy(Scene64 *dest, Scene64 *src, Voxel v){
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
	  
	  memcpy(&dest->data[q],&src->data[p],sizeof(long long)*dx);
	}
      }
    }


    void     Fill(Scene64 *scn, long long value){
      if(value==0)
	memset((void *)scn->data, 0, sizeof(long long)*scn->n);
      else{
	int p;
	for(p=0; p<scn->n; p++)
	  scn->data[p] = value;
      }
      scn->maxval = value;
    }



    long long    GetMaximumValue(Scene64 *scn){
      int p;
      long long Imax;
      
      Imax = INT_MIN;
      for(p=0; p<scn->n; p++) {
	if(scn->data[p] > Imax)
	  Imax = scn->data[p];
      }
      scn->maxval = Imax;
      return(Imax); 
    }
    
    long long   GetMinimumValue(Scene64 *scn){
      int p;
      long long Imin;
      
      Imin = LLONG_MAX;
      for(p=0; p<scn->n; p++){
	if(scn->data[p] < Imin)
	  Imin = scn->data[p];
      }    
      return(Imin); 
    }
    
    
    //The dimensions of the window (i.e., xsize,ysize,zsize) 
    //should be given in millimeters.
    float ComputeWindowSum(Scene64 *Iscn,
			   float xsize, 
			   float ysize, 
			   float zsize,
			   Voxel u){
      Voxel v;
      int dx,dy,dz;
      int q;
      float sum = 0.0;
      
      if(!IsValidVoxel(Iscn, u.c.x, u.c.y, u.c.z))
	return 0.0;
      
      dx = ROUND(xsize/(Iscn->dx*2.0));
      dy = ROUND(ysize/(Iscn->dy*2.0));
      dz = ROUND(zsize/(Iscn->dz*2.0));

      //----------------------------
      v.c.x = MIN(u.c.x+dx, Iscn->xsize-1);
      v.c.y = MIN(u.c.y+dy, Iscn->ysize-1);
      v.c.z = MIN(u.c.z+dz, Iscn->zsize-1);
      q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
      sum = Iscn->data[q];
      
      v.c.x = MAX(u.c.x-dx-1, 0);
      q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
      sum -= Iscn->data[q];
      
      v.c.x = MIN(u.c.x+dx, Iscn->xsize-1);
      v.c.y = MAX(u.c.y-dy-1, 0);
      q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
      sum -= Iscn->data[q];
      
      v.c.y = MIN(u.c.y+dy, Iscn->ysize-1);
      v.c.z = MAX(u.c.z-dz-1, 0);
      q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
      sum -= Iscn->data[q];
      
      //----------------------------
      v.c.x = MAX(u.c.x-dx-1, 0);
      v.c.y = MIN(u.c.y+dy, Iscn->ysize-1);
      v.c.z = MAX(u.c.z-dz-1, 0);
      q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
      sum += Iscn->data[q];
      
      v.c.x = MAX(u.c.x-dx-1, 0);
      v.c.y = MAX(u.c.y-dy-1, 0);
      v.c.z = MIN(u.c.z+dz, Iscn->zsize-1);
      q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
      sum += Iscn->data[q];
      
      v.c.x = MIN(u.c.x+dx, Iscn->xsize-1);
      v.c.y = MAX(u.c.y-dy-1, 0);
      v.c.z = MAX(u.c.z-dz-1, 0);
      q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
      sum += Iscn->data[q];
      
      v.c.x = MAX(u.c.x-dx-1, 0);
      v.c.y = MAX(u.c.y-dy-1, 0);
      v.c.z = MAX(u.c.z-dz-1, 0);
      q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
      sum -= Iscn->data[q];
      
      return sum;
    }


    //The dimensions of the window (i.e., xsize,ysize,zsize) 
    //should be given in millimeters.
    float ComputeWindowDensity(Scene64 *Iscn,
			       float xsize, 
			       float ysize, 
			       float zsize,
			       Voxel u){
      int dx,dy,dz;
      float sum = 0.0;
      float vol;
      
      if(!IsValidVoxel(Iscn, u.c.x, u.c.y, u.c.z))
	return 0.0;
      
      dx = ROUND(xsize/(Iscn->dx*2.0));
      dy = ROUND(ysize/(Iscn->dy*2.0));
      dz = ROUND(zsize/(Iscn->dz*2.0));
      
      vol = (2.0*dx + 1)*(2.0*dy + 1)*(2.0*dz + 1);
      sum = ComputeWindowSum(Iscn, xsize, ysize, zsize, u);
      
      return (sum/vol);
    }


    //The dimensions of the window (i.e., xsize,ysize,zsize) 
    //should be given in millimeters.
    Voxel  FindWindowOfMaximumSum(Scene64 *Iscn, 
				  float xsize, 
				  float ysize, 
				  float zsize){
      Voxel P,u,v;
      int dx,dy,dz;
      int q;
      long long sum,max=INT_MIN;
      
      P.c.x = P.c.y = P.c.z = 0;
      dx = ROUND(xsize/(Iscn->dx*2.0));
      dy = ROUND(ysize/(Iscn->dy*2.0));
      dz = ROUND(zsize/(Iscn->dz*2.0));
 
      for(u.c.x = 0; u.c.x < Iscn->xsize; u.c.x++){
	for(u.c.y = 0; u.c.y < Iscn->ysize; u.c.y++){
	  for(u.c.z = 0; u.c.z < Iscn->zsize; u.c.z++){
	    v.c.x = MIN(u.c.x+dx, Iscn->xsize-1);
	    v.c.y = MIN(u.c.y+dy, Iscn->ysize-1);
	    v.c.z = MIN(u.c.z+dz, Iscn->zsize-1);
	    q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
	    sum = Iscn->data[q];
	    
	    v.c.x = MAX(u.c.x-dx-1, 0);
	    q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
	    sum -= Iscn->data[q];
	    
	    v.c.x = MIN(u.c.x+dx, Iscn->xsize-1);
	    v.c.y = MAX(u.c.y-dy-1, 0);
	    q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
	    sum -= Iscn->data[q];
	    
	    v.c.y = MIN(u.c.y+dy, Iscn->ysize-1);
	    v.c.z = MAX(u.c.z-dz-1, 0);
	    q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
	    sum -= Iscn->data[q];
	    
	    //----------------------------
	    v.c.x = MAX(u.c.x-dx-1, 0);
	    v.c.y = MIN(u.c.y+dy, Iscn->ysize-1);
	    v.c.z = MAX(u.c.z-dz-1, 0);
	    q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
	    sum += Iscn->data[q];
	    
	    v.c.x = MAX(u.c.x-dx-1, 0);
	    v.c.y = MAX(u.c.y-dy-1, 0);
	    v.c.z = MIN(u.c.z+dz, Iscn->zsize-1);
	    q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
	    sum += Iscn->data[q];
	    
	    v.c.x = MIN(u.c.x+dx, Iscn->xsize-1);
	    v.c.y = MAX(u.c.y-dy-1, 0);
	    v.c.z = MAX(u.c.z-dz-1, 0);
	    q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
	    sum += Iscn->data[q];
	    
	    v.c.x = MAX(u.c.x-dx-1, 0);
	    v.c.y = MAX(u.c.y-dy-1, 0);
	    v.c.z = MAX(u.c.z-dz-1, 0);
	    q = GetVoxelAddress(Iscn,v.c.x,v.c.y,v.c.z);
	    sum -= Iscn->data[q];
	    
	    if(sum > max){    
	      max = sum;
	      P = u;
	    }
	  }
	}
      }
      return P;
    }
    


  } //end Scene64 namespace
} //end gft namespace

