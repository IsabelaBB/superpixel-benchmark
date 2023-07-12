
#include "gft_segmentation3.h"

namespace gft{

  namespace Scene8{

  Scene8   *Threshold(Scene8 *scn, int lower, int higher){
    Scene8 *bin=NULL;
    int p;
    bin = Create(scn);
    for(p=0; p<scn->n; p++){
      if((scn->data[p]>=lower)&&(scn->data[p]<=higher))
	bin->data[p]=1;
      else
	bin->data[p]=0;
    }
    bin->maxval = 1;
    return(bin);
  }


  Scene8 *GetBoundaries(Scene8 *scn, AdjRel3::AdjRel3 *A){
    Scene8 *hscn=NULL;
    int p,q,i;
    Voxel u,v;
    hscn = Create(scn);
    for(u.c.z=0; u.c.z<hscn->zsize; u.c.z++){
      for(u.c.y=0; u.c.y<hscn->ysize; u.c.y++){
	for(u.c.x=0; u.c.x<hscn->xsize; u.c.x++){
	  p = GetVoxelAddress(hscn,u);
	  if(scn->data[p] != 0){
	    for(i=1; i<A->n; i++){
	      v.v = u.v + A->d[i].v;
	      if(IsValidVoxel(hscn,v)){
		q = GetVoxelAddress(hscn,v);
		if(scn->data[p] != scn->data[q]){
		  hscn->data[p] = scn->data[p];
		  break;
		}
	      } 
	      else {
		hscn->data[p] = scn->data[p];
		break;
	      }
	    }
	  }
	}
      }
    }
    return(hscn);
  }


  } //end Scene8 namespace


  namespace Scene16{

  Scene8::Scene8   *Threshold(Scene16 *scn, int lower, int higher){
    Scene8::Scene8 *bin=NULL;
    int p;
    bin = Scene8::Create(scn->xsize,scn->ysize,scn->zsize);
    bin->dx = scn->dx;
    bin->dy = scn->dy;
    bin->dz = scn->dz;
    for(p=0; p<scn->n; p++){
      if((scn->data[p]>=lower)&&(scn->data[p]<=higher))
	bin->data[p]=1;
      else
	bin->data[p]=0;
    }
    bin->maxval = 1;
    return(bin);
  }

  } //end Scene16 namespace


  namespace Scene32{

  Scene8::Scene8 *Threshold(Scene32 *scn, int lower, int higher){
    Scene8::Scene8 *bin=NULL;
    int p;
    bin = Scene8::Create(scn->xsize,scn->ysize,scn->zsize);
    bin->dx = scn->dx;
    bin->dy = scn->dy;
    bin->dz = scn->dz;
    for(p=0; p<scn->n; p++){
      if((scn->data[p]>=lower)&&(scn->data[p]<=higher))
	bin->data[p]=1;
      else
	bin->data[p]=0;
    }
    bin->maxval = 1;
    return(bin);
  }

  } //end Scene32 namespace



  namespace Scene{

  Scene8::Scene8 *Threshold(Scene *scn, int lower, int higher){
    switch(scn->nbits){
    case  8:
      return Scene8::Threshold(scn->ptr.scn8, lower,higher);
    case 16:
      return Scene16::Threshold(scn->ptr.scn16, lower,higher);
    case 32:
      return Scene32::Threshold(scn->ptr.scn32, lower,higher);
    }
    return NULL;
  }


  } //end Scene namespace



} //end gft namespace

