
#include "gft_scene8.h"

namespace gft{
  namespace Scene8{

  Scene8 *Create(int xsize,int ysize,int zsize){
    Scene8 *scn=NULL;
    uchar **tmp=NULL;
    int xysize;
    int i,j,p,N;

    scn = (Scene8 *) calloc(1,sizeof(Scene8));
    if(scn == NULL) 
      Error((char *)MSG1,(char *)"Scene8::Create");
 
    scn->xsize  = xsize;
    scn->ysize  = ysize;
    scn->zsize  = zsize;
    scn->dx     = 1.0;
    scn->dy     = 1.0;
    scn->dz     = 1.0;
    scn->maxval = 0;
    scn->n      = xsize*ysize*zsize;

    //For SSE optimization we allocate a multiple of 16.
    N = scn->n;
    if(N%16!=0) N += (16-N%16);

    scn->data   = gft::AllocUCharArray(N);
    if(scn->data==NULL) 
      Error((char *)MSG1,(char *)"Scene8::Create");

    scn->array = (uchar ***) calloc(zsize, sizeof(uchar **));
    if(scn->array==NULL) 
      Error((char *)MSG1,(char *)"Scene8::Create");

    tmp = (uchar **) calloc(zsize*ysize, sizeof(uchar *));
    if(tmp==NULL) 
      Error((char *)MSG1,(char *)"Scene8::Create");

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


  Scene8 *Create(Scene8 *scn){
    Scene8 *new_scn=NULL;
    new_scn = Create(scn->xsize,
		     scn->ysize,
		     scn->zsize);
    new_scn->dx = scn->dx;
    new_scn->dy = scn->dy;
    new_scn->dz = scn->dz;
    return new_scn;
  }


  void     Destroy(Scene8 **scn){
    Scene8 *aux;
    
    aux = *scn;
    if(aux != NULL){
      if(aux->data     != NULL) gft::FreeUCharArray(&aux->data);
      if(aux->array[0] != NULL) free(aux->array[0]);
      if(aux->array    != NULL) free(aux->array);
      free(aux);    
      *scn = NULL;
    }
  }

  Scene8 *Clone(Scene8 *scn){
    Scene8 *aux = Create(scn->xsize,scn->ysize,scn->zsize);
    aux->dx     = scn->dx;
    aux->dy     = scn->dy;
    aux->dz     = scn->dz;
    aux->maxval = scn->maxval;
    memcpy(aux->data,scn->data,sizeof(uchar)*aux->n);
    return aux;
  }


  Scene8 *SubScene(Scene8 *scn, Voxel l, Voxel h){
    return SubScene(scn,
		    l.c.x, l.c.y, l.c.z,
		    h.c.x, h.c.y, h.c.z);
  }


  Scene8 *SubScene(Scene8 *scn,
		   int xl, int yl, int zl,
		   int xh, int yh, int zh){
    Scene8 *sub=NULL;
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
	memcpy(&sub->data[j],&scn->data[i],sizeof(uchar)*sub->xsize);
	j += sub->xsize;
      }
    }
    return(sub);
  }


  void     Copy(Scene8 *dest, Scene8 *src){
    if(dest->xsize!=src->xsize ||
       dest->ysize!=src->ysize ||
       dest->zsize!=src->zsize)
      Error((char *)"Incompatible data",(char *)"Scene8::Copy");
    dest->dx     = src->dx;
    dest->dy     = src->dy;
    dest->dz     = src->dz;
    dest->maxval = src->maxval;
    memcpy(dest->data,src->data,sizeof(uchar)*src->n);
  }


  void     Copy(Scene8 *dest, Scene8 *src, Voxel v){
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

	memcpy(&dest->data[q],&src->data[p],sizeof(uchar)*dx);
      }
    }
  }
  

  void     Fill(Scene8 *scn, uchar value){
    memset((void *)scn->data, value, sizeof(uchar)*scn->n);
    /*
    else{
      int p,i;
      v16qi v;
      v16qi *ptr;
      for(i=0; i<16; i++)
	((uchar *)(&v))[i] = value;
      for(p=0; p<scn->n; p+=16){
	ptr = (v16qi *)(&scn->data[p]);
	*ptr = v;
      }
    }
    */
    scn->maxval = value;
  }



  void     Write(Scene8 *scn, char *filename){
    FILE *fp=NULL;
    int n;

    // Checking file type
    //int len = strlen(filename);
    /*
    if ( (len>=4) && ((strcasecmp(filename + len - 4, ".hdr")==0) || (strcasecmp(filename + len - 4, ".img")==0))) {
      WriteScene_Analyze(scn, filename);
      return;
    }
    if ( (len>=8) && (strcasecmp(filename + len - 8, ".scn.bz2")==0)) {
      WriteCompressedScene(scn, filename);
      return;
    }
    if ( (len<=4) || (strcasecmp(filename + len - 4, ".scn")!=0)) {
      Error(MSG2,"WriteScene: Invalid file name or extension.");
    }
    */
    // Writing the scn file
    fp = fopen(filename,"wb"); 
    if(fp == NULL) 
      Error((char *)MSG2,(char *)"Scene8::Write");

    fprintf(fp,"SCN\n");
    fprintf(fp,"%d %d %d\n",scn->xsize,scn->ysize,scn->zsize);
    fprintf(fp,"%f %f %f\n",scn->dx,scn->dy,scn->dz);
    n = scn->n;
    fprintf(fp,"%d\n",8);
    fwrite(scn->data,sizeof(uchar),n,fp);
    fclose(fp);
  }


  // return the nearest voxel.
  uchar   GetValue_nn(Scene8 *scn, float x, float y, float z){
    Voxel v;
    if(x>scn->xsize-1||y>scn->ysize-1||z>scn->zsize-1||
       x<(float)0||y<(float)0||z<(float)0)
      Error((char *)"Out-of-bounds",
	    (char *)"Scene8::GetValue_nn");

    v.c.x=(int)(x+0.5);
    v.c.y=(int)(y+0.5);
    v.c.z=(int)(z+0.5);
    return scn->array[v.c.z][v.c.y][v.c.x];
  }


  uchar          GetMaximumValue(Scene8 *scn){
    int p;
    uchar Imax=0;
    for(p=0; p<scn->n; p++) {
      if(scn->data[p] > Imax)
	Imax = scn->data[p];
    }
    scn->maxval = Imax;
    return(Imax); 
  }

  uchar          GetMinimumValue(Scene8 *scn){
    int p;
    uchar Imin=255;
    for(p=0; p<scn->n; p++){
      if(scn->data[p] < Imin)
	Imin = scn->data[p];
    }    
    return(Imin); 
  }


  void     MBB(Scene8 *scn, Voxel *l, Voxel *h){
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


  Scene8 *MBB(Scene8 *scn){
    Voxel lower,higher;
    Scene8 *mbb=NULL;
    MBB(scn, &lower, &higher);
    mbb = SubScene(scn,
		   lower.c.x,  lower.c.y,  lower.c.z,
		   higher.c.x, higher.c.y, higher.c.z);
    return(mbb);
  }


  Scene8 *AddFrame(Scene8 *scn,  int sz, uchar value){
    Scene8 *fscn;
    int y, z,nbytes,offset1, offset2;
    uchar *dst,*src;
    
    fscn = Create(scn->xsize+(2*sz),
		  scn->ysize+(2*sz), 
		  scn->zsize+(2*sz));
    fscn->dx = scn->dx;
    fscn->dy = scn->dy;
    fscn->dz = scn->dz;

    Fill(fscn,value);
    nbytes = sizeof(uchar)*scn->xsize;
  
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


  Scene8 *RemFrame(Scene8 *fscn, int sz){
    Scene8 *scn;
    int y,z,nbytes,offset;
    uchar *dst,*src;
    
    scn = Create(fscn->xsize-(2*sz),
		 fscn->ysize-(2*sz),
		 fscn->zsize-(2*sz));
    scn->dx = fscn->dx;
    scn->dy = fscn->dy;
    scn->dz = fscn->dz;

    nbytes = sizeof(uchar)*scn->xsize;  
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


  Scene32::Scene32 *ConvertTo32(Scene8 *scn){
    Scene32::Scene32 *scn32;
    int p;

    scn32 = Scene32::Create(scn->xsize,
			    scn->ysize,
			    scn->zsize);
    scn32->dx = scn->dx;
    scn32->dy = scn->dy;
    scn32->dz = scn->dz;
    for(p=0; p<scn->n; p++){
      scn32->data[p] = (int)scn->data[p];
    }
    return scn32;
  }


  Scene16::Scene16 *ConvertTo16(Scene8 *scn){
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



  } //end Scene8 namespace
} //end gft namespace

