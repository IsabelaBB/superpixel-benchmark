
#include "gft_seedmap3.h"

namespace gft{
  namespace AdjSeedmap3{

    AdjSeedmap3 *Create(int nobjs){
      AdjSeedmap3 *asmap=NULL;
      int l;
      asmap = (AdjSeedmap3 *) calloc(1,sizeof(AdjSeedmap3));
      if(asmap == NULL)
	gft::Error((char *)MSG1,
		   (char *)"AdjSeedmap3::Create");
      
      asmap->disp = AdjRel3::Create(nobjs+1);
      asmap->object = ((AdjRegion3::AdjRegion3 **)
		       calloc(nobjs+1,
			      sizeof(AdjRegion3::AdjRegion3 *)));
      asmap->obj_border = ((AdjRegion3::AdjRegion3 **)
			   calloc(nobjs+1,
				  sizeof(AdjRegion3::AdjRegion3 *)));
      asmap->bkg_border = ((AdjRegion3::AdjRegion3 **)
			   calloc(nobjs+1,
				  sizeof(AdjRegion3::AdjRegion3 *)));
      asmap->uncertainty = ((AdjRegion3::AdjRegion3 **)
			    calloc(nobjs+1,
				   sizeof(AdjRegion3::AdjRegion3 *)));
      
      if(asmap->object==NULL || 
	 asmap->obj_border==NULL ||
	 asmap->bkg_border==NULL ||
	 asmap->uncertainty==NULL)
	gft::Error((char *)MSG1,
		   (char *)"AdjSeedmap3::Create");
      
      for(l=0; l<=nobjs; l++){
	asmap->object[l] = NULL;
	asmap->obj_border[l] = NULL;
	asmap->bkg_border[l] = NULL;
	asmap->uncertainty[l] = NULL;
      }
      asmap->nobjs = nobjs;
      return asmap;
    }


    void Destroy(AdjSeedmap3 **asmap){
      AdjSeedmap3 *aux;
      int i;
      aux = *asmap;
      if(aux != NULL){
	if(aux->disp!=NULL)
	  AdjRel3::Destroy(&aux->disp);
	if(aux->uncertainty!=NULL){
	  for(i=0; i<=aux->nobjs; i++)
	    if(aux->uncertainty[i]!=NULL)
	      AdjRegion3::Destroy(&aux->uncertainty[i]);
	  free(aux->uncertainty);
	}
	if(aux->object!=NULL){
	  for(i=0; i<=aux->nobjs; i++)
	    if(aux->object[i]!=NULL)
	      AdjRegion3::Destroy(&aux->object[i]);
	  free(aux->object);
	}
	if(aux->obj_border!=NULL){
	  for(i=0; i<=aux->nobjs; i++)
	    if(aux->obj_border[i]!=NULL)
	      AdjRegion3::Destroy(&aux->obj_border[i]);
	  free(aux->obj_border);
	}
	if(aux->bkg_border!=NULL){
	  for(i=0; i<=aux->nobjs; i++)
	    if(aux->bkg_border[i]!=NULL)
	      AdjRegion3::Destroy(&aux->bkg_border[i]);
	  free(aux->bkg_border);
	}
	free(aux);
	*asmap = NULL;
      }
    }


    AdjSeedmap3 *Create(RegionCloud3::RegionCloud3 *rcloud){
      AdjSeedmap3 *asmap=NULL;
      Scene8::Scene8 *mask=NULL,*fmask=NULL;
      Scene8::Scene8 *bmask=NULL,*fbmask=NULL;
      Scene32::Scene32 *prob=NULL;
      int nimgs,nobjs,l,p;
      AdjRel3::AdjRel3 *A=NULL;
      Voxel Ref;
      
      A = AdjRel3::Spheric(1.0);
      nobjs = rcloud->nobjs;
      nimgs = rcloud->nimages;
      asmap = Create(nobjs);
      
      for(l=1; l<=nobjs; l++){
	//-----------------------------------------
	prob = Scene32::Create((rcloud->prob[l])->xsize,
			       (rcloud->prob[l])->ysize,
			       (rcloud->prob[l])->zsize);
	for(p=0; p<(rcloud->prob[l])->n; p++)
	  prob->data[p] = (rcloud->prob[l])->data[p];
	//-----------------------------------------

	(asmap->disp)->d[l].axis.x = (rcloud->disp)->dx[l];
	(asmap->disp)->d[l].axis.y = (rcloud->disp)->dy[l];
	(asmap->disp)->d[l].axis.z = (rcloud->disp)->dz[l];
	
	mask = Scene32::Threshold(prob, MAX_PROB, MAX_PROB);
	Ref.c.x = mask->xsize/2;
	Ref.c.y = mask->ysize/2;
	Ref.c.z = mask->zsize/2;
	asmap->object[l] = AdjRegion3::Create(mask, Ref);
	
	bmask = Scene8::GetBoundaries(mask, A);
	asmap->obj_border[l] = AdjRegion3::Create(bmask, Ref);
	Scene8::Destroy(&bmask);
	Scene8::Destroy(&mask);
	
	mask = Scene32::Threshold(prob, 1, MAX_PROB-1);
	asmap->uncertainty[l] = AdjRegion3::Create(mask, Ref);
	Scene8::Destroy(&mask);
	
	mask  = Scene32::Threshold(prob, 0, 0);
	fmask = Scene8::AddFrame(mask, 2, 1);
	fbmask = Scene8::GetBoundaries(fmask, A);
	bmask = Scene8::RemFrame(fbmask, 1);
	Ref.c.x += 1;
	Ref.c.y += 1;
	Ref.c.z += 1;
	asmap->bkg_border[l] = AdjRegion3::Create(bmask, Ref);
	Scene8::Destroy(&fbmask);
	Scene8::Destroy(&bmask);
	Scene8::Destroy(&fmask);
	Scene8::Destroy(&mask);
	//-------------------------------
	Scene32::Destroy(&prob);
      }
      AdjRel3::Destroy(&A);
      
      return asmap;
    }

    
    
    void DrawObject(Scene8::Scene8 *scn,
		    AdjSeedmap3 *asmap,
		    Voxel u,
		    int l,
		    uchar val){
      Voxel du;
      du.v = u.v + (asmap->disp)->d[l].v;
      AdjRegion3::Draw(asmap->object[l],
		       scn, du, val);
    }
    
    
    void DrawObjBorder(Scene8::Scene8 *scn,
		       AdjSeedmap3 *asmap,
		       Voxel u,
		       int l,
		       uchar val){
      Voxel du;
      du.v = u.v + (asmap->disp)->d[l].v;
      AdjRegion3::Draw(asmap->obj_border[l],
		       scn, du, val);
    }
    
    
    void DrawBkgBorder(Scene8::Scene8 *scn,
		       AdjSeedmap3 *asmap,
		       Voxel u,
		       int l,
		       uchar val){
      Voxel du;
      du.v = u.v + (asmap->disp)->d[l].v;
      AdjRegion3::Draw(asmap->bkg_border[l],
		       scn, du, val);
    }
    

    void DrawUncertainty(Scene8::Scene8 *scn,
			 AdjSeedmap3 *asmap,
			 Voxel u,
			 int l,
			 uchar val){
      Voxel du;
      int p;
      du.v = u.v + (asmap->disp)->d[l].v;
      if(!AdjRegion3::FitInside(asmap->bkg_border[l],
				du, scn, 1)){
	AdjRegion3::Draw(asmap->uncertainty[l],
			 scn, du, val);
      }
      else{
	p = Scene8::GetVoxelAddress(scn,du);
	AdjRegion3::DrawOpt(asmap->uncertainty[l],
			    scn, p, val);
      }
    }


    void CopyUncertainty(Scene8::Scene8 *dest,
			 Scene8::Scene8 *src,
			 AdjSeedmap3 *asmap,
			 Voxel u,
			 int l){
      AdjRegion3::AdjRegion3 *S=NULL;
      Voxel du,v;
      int i,p,q;
      du.v = u.v + (asmap->disp)->d[l].v;
      if(!AdjRegion3::FitInside(asmap->bkg_border[l],
				du, dest, 1)){
	S = asmap->uncertainty[l];
	for(i=0; i<S->n; i++){
	  v.v = du.v + S->d[i].v;
	  if(Scene8::IsValidVoxel(dest,v)){
	    p = Scene8::GetVoxelAddress(dest,v);
	    dest->data[p] = src->data[p];
	  }
	}
      }
      else{ //Cloud3FitInside
	p = Scene8::GetVoxelAddress(dest,du);
	S = asmap->uncertainty[l];
	AdjRegion3::Optimize(S, dest);
	for(i=0; i<S->n; i++){
	  q = p + S->dp[i];
	  dest->data[q] = src->data[q];
	}
      }
    }
    
    
    void AddUncertainty(Scene8::Scene8 *dest,
			Scene8::Scene8 *src,
			AdjSeedmap3 *asmap,
			Voxel u,
			int l){
      AdjRegion3::AdjRegion3 *S=NULL;
      Voxel du,v;
      int i,p,q;
      du.v = u.v + (asmap->disp)->d[l].v;
      if(!AdjRegion3::FitInside(asmap->bkg_border[l],
				du, dest, 1)){
	S = asmap->uncertainty[l];
	for(i=0; i<S->n; i++){
	  v.v = du.v + S->d[i].v;
	  if(Scene8::IsValidVoxel(dest,v)){
	    p = Scene8::GetVoxelAddress(dest,v);
	    dest->data[p] += src->data[p];
	  }
	}
      }
      else{ //Cloud3FitInside
	p = Scene8::GetVoxelAddress(dest,du);
	S = asmap->uncertainty[l];
	AdjRegion3::Optimize(S, dest);
	for(i=0; i<S->n; i++){
	  q = p + S->dp[i];
	  dest->data[q] += src->data[q];
	}
      }
    }
    
    
    /*
      void CloudArcWeight3(Scene *arcw, 
      Scene *grad,
      Voxel u,
      BorderCloud3 *bcloud,
      Seedmap3 *smap,
      AdjSeedmap3 *asmap,
      int l,
      float w){
      Cloud3 *S=NULL;
      Scene *prob=NULL;
      Voxel v,b,du;
      int i,p,q,run;
      int dx,dy,dz;
      float v1,v2;
      
      du.x = u.x + (asmap->disp)->dx[l];
      du.y = u.y + (asmap->disp)->dy[l];
      du.z = u.z + (asmap->disp)->dz[l];
      for(run=1; run<=3; run++){
      if(run==1) S = asmap->uncertainty[l];
      if(run==2) S = asmap->obj_border[l];
      if(run==3) S = asmap->bkg_border[l];
      
      for(i=0; i<S->n; i++){
      dx = S->dx[i];
      dy = S->dy[i];
      dz = S->dz[i];
      
      v.x = du.x + dx;
      v.y = du.y + dy;
      v.z = du.z + dz;

      if(ValidVoxel(arcw,v.x,v.y,v.z)){
      p = VoxelAddress(arcw,v.x,v.y,v.z);
      
	v2 = (float)grad->data[p];
	
	prob = bcloud->prob[l];
	b.x = prob->xsize/2 + dx;
	b.y = prob->ysize/2 + dy;
	b.z = prob->zsize/2 + dz;
	if(ValidVoxel(prob,b.x,b.y,b.z)){
	q = VoxelAddress(prob,b.x,b.y,b.z);
	
	v1 = (float)prob->data[q];
	}
	else v1 = 0.0;
	
	arcw->data[p] = ROUND(w*v1+(1.0-w)*v2);
	}
	}
	}
	}
    */


    void CloudArcWeight(Scene16::Scene16 *arcw,
			Scene16::Scene16 *wobj,
			Gradient3::Gradient3 *grad,
			Voxel u,
			BorderCloud3::BorderCloud3 *bcloud,
			AdjSeedmap3 *asmap,
			int l, float w){
      AdjRegion3::AdjRegion3 *S=NULL;
      Gradient3::Gradient3 *prob=NULL;
      Voxel v,b,du;
      int i,p,q,run;
      int dx,dy,dz;
      float v1,v2,dotproduct,G;//,cos,gr;
      
      prob = bcloud->prob[l];
      Gradient3::ComputeMagnitude(prob);
      
      G = (((float)(prob->mag)->maxval)*
	   ((float)(grad->mag)->maxval));
      
      if(G==0.0) G = 0.01;
      du.v = u.v + (asmap->disp)->d[l].v;
      for(run=1; run<=3; run++){
	if(run==1) S = asmap->uncertainty[l];
	if(run==2) S = asmap->obj_border[l];
	if(run==3) S = asmap->bkg_border[l];
	
	for(i=0; i<S->n; i++){
	  dx = S->d[i].axis.x;
	  dy = S->d[i].axis.y;
	  dz = S->d[i].axis.z;
	  
	  v.v = du.v + S->d[i].v;
	  
	  if(Scene16::IsValidVoxel(arcw,v)){
	    p = Scene16::GetVoxelAddress(arcw,v);

	    v2 = (float)wobj->data[p];
	    
	    b.c.x = (prob->Gx)->xsize/2 + dx;
	    b.c.y = (prob->Gx)->ysize/2 + dy;
	    b.c.z = (prob->Gx)->zsize/2 + dz;
	    if(Scene32::IsValidVoxel(prob->Gx,b.c.x,b.c.y,b.c.z)){
	      q = Scene32::GetVoxelAddress(prob->Gx,b.c.x,b.c.y,b.c.z);
	      v1 = (float)(prob->mag)->data[q];
	      
	      //-------------------
	      //PowerEnhancement:
	      v1 = (v1*v1)/((float)(prob->mag)->maxval);
	      //-------------------
	      
	      dotproduct  = grad->Gx->data[p]*prob->Gx->data[q];
	      dotproduct += grad->Gy->data[p]*prob->Gy->data[q];
	      dotproduct += grad->Gz->data[p]*prob->Gz->data[q];
	      //gr = (float)(grad->mag)->data[p];
	      //if(v1*gr>0.0) cos = dotproduct/(v1*gr);
	      //else          cos = 0.0;
	    }
	    else{ 
	      v1  = 0.0;
	      dotproduct = 0.0;
	      //cos = 0.0;
	      //gr  = 0.0;
	    }
	    
	    //if(cos>0.0) cos = 0.0;
	    if(dotproduct>0.0) dotproduct = 0.0;
	    arcw->data[p] = ROUND((w*v1+(1.0-w)*v2)/(1.0-2.*dotproduct/G));
	  }
	}
      }
    }


  } //end AdjSeedmap3 namespace
} //end gft namespace



