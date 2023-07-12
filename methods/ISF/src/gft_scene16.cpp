
#include "gft_scene16.h"

namespace gft{
  namespace Scene16{

  Scene16 *Create(int xsize,int ysize,int zsize){
    Scene16 *scn=NULL;
    ushort **tmp=NULL;
    int xysize;
    int i,j,p,N;

    scn = (Scene16 *) calloc(1,sizeof(Scene16));
    if(scn == NULL) 
      Error((char *)MSG1,(char *)"Scene16::Create");
 
    scn->xsize  = xsize;
    scn->ysize  = ysize;
    scn->zsize  = zsize;
    scn->dx     = 1.0;
    scn->dy     = 1.0;
    scn->dz     = 1.0;
    scn->maxval = 0;
    scn->n      = xsize*ysize*zsize;

    //For SSE optimization we allocate a multiple of 8.
    N = scn->n;
    if(N%8!=0) N += (8-N%8);

    scn->data   = gft::AllocUShortArray(N);
    if(scn->data==NULL) 
      Error((char *)MSG1,(char *)"Scene16::Create");

    scn->array = (ushort ***) calloc(zsize, sizeof(ushort **));
    if(scn->array==NULL) 
      Error((char *)MSG1,(char *)"Scene16::Create");

    tmp = (ushort **) calloc(zsize*ysize, sizeof(ushort *));
    if(tmp==NULL) 
      Error((char *)MSG1,(char *)"Scene16::Create");

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


  Scene16 *Create(Scene16 *scn){
    Scene16 *new_scn=NULL;
    new_scn = Create(scn->xsize,
		     scn->ysize,
		     scn->zsize);
    new_scn->dx = scn->dx;
    new_scn->dy = scn->dy;
    new_scn->dz = scn->dz;
    return new_scn;
  }


  void     Destroy(Scene16 **scn){
    Scene16 *aux;
    
    aux = *scn;
    if(aux != NULL){
      if(aux->data     != NULL) gft::FreeUShortArray(&aux->data);
      if(aux->array[0] != NULL) free(aux->array[0]);
      if(aux->array    != NULL) free(aux->array);
      free(aux);    
      *scn = NULL;
    }
  }

  Scene16 *Clone(Scene16 *scn){
    Scene16 *aux = Create(scn->xsize,scn->ysize,scn->zsize);
    aux->dx     = scn->dx;
    aux->dy     = scn->dy;
    aux->dz     = scn->dz;
    aux->maxval = scn->maxval;
    memcpy(aux->data,scn->data,sizeof(ushort)*aux->n);
    return aux;
  }


  Scene16 *SubScene(Scene16 *scn, Voxel l, Voxel h){
    return SubScene(scn,
		    l.c.x, l.c.y, l.c.z,
		    h.c.x, h.c.y, h.c.z);
  }


  Scene16 *SubScene(Scene16 *scn,
		    int xl, int yl, int zl,
		    int xh, int yh, int zh){
    Scene16 *sub=NULL;
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
	memcpy(&sub->data[j],&scn->data[i],sizeof(ushort)*sub->xsize);
	j += sub->xsize;
      }
    }
    return(sub);
  }


  void     Copy(Scene16 *dest, Scene16 *src){
    if(dest->xsize!=src->xsize ||
       dest->ysize!=src->ysize ||
       dest->zsize!=src->zsize)
      Error((char *)"Incompatible data",(char *)"Scene16::Copy");
    dest->dx     = src->dx;
    dest->dy     = src->dy;
    dest->dz     = src->dz;
    dest->maxval = src->maxval;
    memcpy(dest->data,src->data,sizeof(ushort)*src->n);
  }


  void     Copy(Scene16 *dest, Scene16 *src, Voxel v){
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

	memcpy(&dest->data[q],&src->data[p],sizeof(ushort)*dx);
      }
    }
  }


  void     Fill(Scene16 *scn, ushort value){
    if(value==0)
      memset((void *)scn->data, 0, sizeof(ushort)*scn->n);
    else{
      int p,i;
      v8hi v;
      v8hi *ptr;
      for(i=0; i<8; i++)
	((ushort *)(&v))[i] = value;
      for(p=0; p<scn->n; p+=8){
	ptr = (v8hi *)(&scn->data[p]);
	*ptr = v;
      }
    }
    scn->maxval = value;
  }



  void     Write(Scene16 *scn, char *filename){
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
      Error((char *)MSG2,(char *)"Scene16::Write");

    fprintf(fp,"SCN\n");
    fprintf(fp,"%d %d %d\n",scn->xsize,scn->ysize,scn->zsize);
    fprintf(fp,"%f %f %f\n",scn->dx,scn->dy,scn->dz);
  
    n = scn->n;
    fprintf(fp,"%d\n",16);
    fwrite(scn->data,sizeof(ushort),n,fp);
    fclose(fp);
  }


  // return the nearest voxel.
  ushort GetValue_nn(Scene16 *scn, float x, float y, float z){
    Voxel v;
    if(x>scn->xsize-1||y>scn->ysize-1||z>scn->zsize-1||
       x<(float)0||y<(float)0||z<(float)0)
      Error((char *)"Out-of-bounds",
	    (char *)"Scene16::GetValue_nn");

    v.c.x=(int)(x+0.5);
    v.c.y=(int)(y+0.5);
    v.c.z=(int)(z+0.5);
    return scn->array[v.c.z][v.c.y][v.c.x];
  }


  ushort    GetMaximumValue(Scene16 *scn){
    int p;
    ushort Imax = 0;
    
    for(p=0; p<scn->n; p++) {
      if(scn->data[p] > Imax)
	Imax = scn->data[p];
    }
    scn->maxval = Imax;
    return(Imax); 
  }

  ushort    GetMinimumValue(Scene16 *scn){
    int p;
    ushort Imin = USHRT_MAX;

    for(p=0; p<scn->n; p++){
      if(scn->data[p] < Imin)
	Imin = scn->data[p];
    }    
    return(Imin); 
  }


  void     MBB(Scene16 *scn, Voxel *l, Voxel *h){
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


  Scene16 *MBB(Scene16 *scn){
    Voxel lower,higher;
    Scene16 *mbb=NULL;
    MBB(scn, &lower, &higher);
    mbb = SubScene(scn,
		   lower.c.x,  lower.c.y,  lower.c.z,
		   higher.c.x, higher.c.y, higher.c.z);
    return(mbb);
  }


  Scene16 *AddFrame(Scene16 *scn,  int sz, ushort value){
    Scene16 *fscn;
    int y, z,nbytes,offset1, offset2;
    ushort *dst,*src;
    
    fscn = Create(scn->xsize+(2*sz),
		  scn->ysize+(2*sz), 
		  scn->zsize+(2*sz));
    fscn->dx = scn->dx;
    fscn->dy = scn->dy;
    fscn->dz = scn->dz;

    Fill(fscn,value);
    nbytes = sizeof(ushort)*scn->xsize;
  
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


  Scene16 *RemFrame(Scene16 *fscn, int sz){
    Scene16 *scn;
    int y,z,nbytes,offset;
    ushort *dst,*src;
    
    scn = Create(fscn->xsize-(2*sz),
		 fscn->ysize-(2*sz),
		 fscn->zsize-(2*sz));
    scn->dx = fscn->dx;
    scn->dy = fscn->dy;
    scn->dz = fscn->dz;

    nbytes = sizeof(ushort)*scn->xsize;  
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


  Scene32::Scene32 *ConvertTo32(Scene16 *scn){
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


  Scene8::Scene8   *ConvertTo8(Scene16 *scn){
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



  } //end Scene16 namespace
} //end gft namespace

