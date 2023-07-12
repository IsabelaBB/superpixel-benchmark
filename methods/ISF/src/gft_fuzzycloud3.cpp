
#include "fuzzycloud3.h"


namespace gft{
  namespace RegionCloud3{


    RegionCloud3 *ByLabelList(FileList::FileList *L){
      int nimages,i,l,Lmax,p,q;
      RegionCloud3 *rcloud=NULL;
      Scene32::Scene32 *label;
      Scene8::Scene8 *mask;
      int Sx,Sy,Sz,s;     //sizes.
      gft::Voxel **C = NULL;   //centroids.
      gft::Voxel **MBBl = NULL;
      gft::Voxel **MBBh = NULL;
      gft::Voxel Vl,Vh,u,v,c;
      float dx=1.0,dy=1.0,dz=1.0;
      float sum_x,sum_y,sum_z,sum;
      
      //Find the maximum label:
      Lmax = 0;
      nimages = L->n;
      for(i=0; i<nimages; i++){
	label = Scene32::Read( FileList::GetFile(L, i) );
	Lmax = MAX(Lmax, Scene32::GetMaximumValue(label));
	dx = label->dx;
	dy = label->dy;
	dz = label->dz;
	Scene32::Destroy(&label);
      }
      
      //Initialize structures:
      C = (gft::Voxel **) calloc(nimages, sizeof(gft::Voxel *));
      MBBl = (gft::Voxel **) calloc(nimages, sizeof(gft::Voxel *));
      MBBh = (gft::Voxel **) calloc(nimages, sizeof(gft::Voxel *));
      if(C==NULL) gft::Error((char *)MSG1,
			     (char *)"RegionCloud3::ByLabelList");
      for(i=0; i<nimages; i++){
	C[i] = (gft::Voxel *) calloc(Lmax+1, sizeof(gft::Voxel));
	if(C[i]==NULL) gft::Error((char *)MSG1,
				  (char *)"RegionCloud3::ByLabelList");
	
	MBBl[i] = (gft::Voxel *) calloc(Lmax+1, sizeof(gft::Voxel));
	if(MBBl[i]==NULL) gft::Error((char *)MSG1,
				     (char *)"RegionCloud3::ByLabelList");
	
	MBBh[i] = (gft::Voxel *) calloc(Lmax+1, sizeof(gft::Voxel));
	if(MBBh[i]==NULL) gft::Error((char *)MSG1,
				     (char *)"RegionCloud3::ByLabelList");
      }
      
      //Compute centroids and MBBs:
      for(i=0; i<nimages; i++){
	label = Scene32::Read( FileList::GetFile(L, i) );
	
	for(l=1; l<=Lmax; l++){
	  mask = Scene32::Threshold(label, l, l);

	  c = Scene8::Centroid(mask);
	  Scene8::MBB(mask, &Vl, &Vh);

	  C[i][l] = c;
	  MBBl[i][l] = Vl;
	  MBBh[i][l] = Vh;
	  Scene8::Destroy(&mask);
	}
	mask = Scene32::Threshold(label, 1, INT_MAX);
	C[i][0] = Scene8::Centroid(mask);

	Scene8::Destroy(&mask);
	Scene32::Destroy(&label);
      }
      
      //Create RegionCloud3:
      rcloud = (RegionCloud3 *)calloc(1,sizeof(RegionCloud3));
      if(rcloud==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::ByLabelList");
      rcloud->prob = (Scene **)calloc(Lmax+1,sizeof(Scene *));
      if(rcloud->prob==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::LabelList");
      rcloud->disp  = AdjRel3::Create(Lmax+1);
      rcloud->fdisp = AdjRel3f::Create(Lmax+1);
      rcloud->nobjs = Lmax;
      rcloud->nimages = nimages;
      
      for(l=1; l<=Lmax; l++){
	sum_x = sum_y = sum_z = 0.0;
	for(i=0; i<nimages; i++){
	  sum_x += (C[i][l].x - C[i][0].x);
	  sum_y += (C[i][l].y - C[i][0].y);
	  sum_z += (C[i][l].z - C[i][0].z);
	}
	(rcloud->fdisp)->dx[l] = sum_x/(float)nimages;
	(rcloud->fdisp)->dy[l] = sum_y/(float)nimages;
	(rcloud->fdisp)->dz[l] = sum_z/(float)nimages;
	
	(rcloud->disp)->dx[l] = ROUND(sum_x/(float)nimages);
	(rcloud->disp)->dy[l] = ROUND(sum_y/(float)nimages);
	(rcloud->disp)->dz[l] = ROUND(sum_z/(float)nimages);
	
	Sx = Sy = Sz = 0;
	for(i=0; i<nimages; i++){
	  c.x = C[i][l].x;
	  c.y = C[i][l].y;
	  c.z = C[i][l].z;
	  Vl = MBBl[i][l];
	  Vh = MBBh[i][l];
	  
	  s = MAX(c.x-Vl.x,Vh.x-c.x)*2+1;
	  if(s > Sx) Sx = s;
	  s = MAX(c.y-Vl.y,Vh.y-c.y)*2+1;
	  if(s > Sy) Sy = s;
	  s = MAX(c.z-Vl.z,Vh.z-c.z)*2+1;
	  if(s > Sz) Sz = s;
	}
	rcloud->prob[l] = Scene32::Create(Sx, Sy, Sz);
	rcloud->prob[l]->dx =  dx;
	rcloud->prob[l]->dy =  dy;
	rcloud->prob[l]->dz =  dz;
      }
      
      //Compute Sum:
      for(i=0; i<nimages; i++){
	label = Scene32::Read( FileList::GetFile(L, i) );
	
	for(p=0; p<label->n; p++){
	  l = label->data[p];
	  if(l==0) continue;
	  u.x = Scene32::GetAddressX(label, p);
	  u.y = Scene32::GetAddressY(label, p);
	  u.z = Scene32::GetAddressZ(label, p);
	  
	  u.x -= C[i][l].x;
	  u.y -= C[i][l].y;
	  u.z -= C[i][l].z;
	  
	  v.x = u.x + (rcloud->prob[l])->xsize/2;
	  v.y = u.y + (rcloud->prob[l])->ysize/2;
	  v.z = u.z + (rcloud->prob[l])->zsize/2;

	  if(Scene32::IsValidVoxel(rcloud->prob[l], v)){
	    q = Scene32::GetVoxelAddress(rcloud->prob[l], v);
	    (rcloud->prob[l])->data[q] += 1;
	  }
	}
	Scene32::Destroy(&label);
      }
      
      //Sum to Probability:
      for(l=1; l<=Lmax; l++){
	for(p=0; p<(rcloud->prob[l])->n; p++){
	  sum = (float)(rcloud->prob[l])->data[p];
	  (rcloud->prob[l])->data[p] = ROUND(MAX_PROB*(sum/nimages));
	}
      }
      
      //Free memory:
      for(i=0; i<nimages; i++){
	free(C[i]);
	free(MBBl[i]);
	free(MBBh[i]);
      }
      free(C);
      free(MBBl);
      free(MBBh);
      
      return rcloud;
    }
    


    void Destroy(RegionCloud3 **rcloud){
      RegionCloud3 *aux;
      int i;
      
      aux = *rcloud;
      if(aux != NULL){
	if(aux->disp!=NULL)
	  AdjRel3::Destroy(&aux->disp);
	if(aux->fdisp!=NULL)
	  AdjRel3f::Destroy(&aux->fdisp);
	
	if(aux->prob!=NULL){
	  for(i=1; i<=aux->nobjs; i++)
	    if(aux->prob[i]!=NULL)
	      Scene32::Destroy(&aux->prob[i]);
	  free(aux->prob);
	}
	free(aux);
	*rcloud = NULL;
      }
    }
    

    void GetVoxelSize(RegionCloud3 *rcloud,
		      float *dx, 
		      float *dy, 
		      float *dz){
      *dx = (rcloud->prob[1])->dx;
      *dy = (rcloud->prob[1])->dy;
      *dz = (rcloud->prob[1])->dz;
    }
    
    
    RegionCloud3 *Subsampling(RegionCloud3 *rcloud){
      RegionCloud3 *sub=NULL;
      Scene32::Scene32 *tmp=NULL,*blur=NULL;
      int nobjs,l;
      
      nobjs = rcloud->nobjs;
      
      //Create RegionCloud3:
      sub = (RegionCloud3 *)calloc(1,sizeof(RegionCloud3));
      if(sub==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::Subsampling");
      sub->prob = (Scene32::Scene32 **)calloc(nobjs+1,sizeof(Scene32::Scene32 *));
      if(sub->prob==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::Subsampling");
      sub->disp  = AdjRel3::Create(nobjs+1);
      sub->fdisp = AdjRel3f::Create(nobjs+1);
      sub->nobjs = rcloud->nobjs;
      sub->nimages = rcloud->nimages;
      
      for(l=1; l<=nobjs; l++){
	(sub->fdisp)->dx[l] = (rcloud->fdisp)->dx[l]/2.;
	(sub->fdisp)->dy[l] = (rcloud->fdisp)->dy[l]/2.;
	(sub->fdisp)->dz[l] = (rcloud->fdisp)->dz[l]/2.;
	
	(sub->disp)->dx[l] = ROUND((sub->fdisp)->dx[l]);
	(sub->disp)->dy[l] = ROUND((sub->fdisp)->dy[l]);
	(sub->disp)->dz[l] = ROUND((sub->fdisp)->dz[l]);
	
	tmp  = Scene32::AddFrame(rcloud->prob[l], 1, 0);
	blur = FastGaussianBlur3(tmp);
	sub->prob[l] = Subsampling3(blur);
	Scene32::Destroy(&blur);
	Scene32::Destroy(&tmp);
      }
      return sub;
    }


    RegionCloud3 *LinearInterp(RegionCloud3 *rcloud,
			       float dx,float dy,float dz){
      RegionCloud3 *interp=NULL;
      float odx,ody,odz;
      int nobjs,l;
      
      nobjs = rcloud->nobjs;
      
      //Create RegionCloud3:
      interp = (RegionCloud3 *)calloc(1,sizeof(RegionCloud3));
      if(interp==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::LinearInterp");
      interp->prob = (Scene **)calloc(nobjs+1,sizeof(Scene *));
      if(interp->prob==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::LinearInterp");
      interp->disp  = AdjRel3::Create(nobjs+1);
      interp->fdisp = AdjRel3f::Create(nobjs+1);
      interp->nobjs = rcloud->nobjs;
      interp->nimages = rcloud->nimages;
      
      for(l=1; l<=nobjs; l++){
	odx = rcloud->prob[l]->dx;
	ody = rcloud->prob[l]->dy;
	odz = rcloud->prob[l]->dz;
	(interp->fdisp)->dx[l] = (rcloud->fdisp)->dx[l]*odx/dx;
	(interp->fdisp)->dy[l] = (rcloud->fdisp)->dy[l]*ody/dy;
	(interp->fdisp)->dz[l] = (rcloud->fdisp)->dz[l]*odz/dz;
	
	(interp->disp)->dx[l] = ROUND((interp->fdisp)->dx[l]);
	(interp->disp)->dy[l] = ROUND((interp->fdisp)->dy[l]);
	(interp->disp)->dz[l] = ROUND((interp->fdisp)->dz[l]);
	
	interp->prob[l] = FastLinearInterpCentr3(rcloud->prob[l], 
						 dx,dy,dz);
      }
      return interp;
    }


    RegionCloud3 *GaussianBlur(RegionCloud3 *rcloud){
      RegionCloud3 *blur=NULL;
      Scene32 *tmp=NULL;
      int nobjs,l;
      
      nobjs = rcloud->nobjs;
      
      //Create RegionCloud3:
      blur = (RegionCloud3 *)calloc(1,sizeof(RegionCloud3));
      if(blur==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::GaussianBlur");
      blur->prob = (Scene32 **)calloc(nobjs+1,sizeof(Scene32 *));
      if(blur->prob==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::GaussianBlur");
      blur->disp  = AdjRel3::Create(nobjs+1);
      blur->fdisp = AdjRel3f::Create(nobjs+1);
      blur->nobjs = rcloud->nobjs;
      blur->nimages = rcloud->nimages;
      
      for(l=1; l<=nobjs; l++){
	(blur->fdisp)->dx[l] = (rcloud->fdisp)->dx[l];
	(blur->fdisp)->dy[l] = (rcloud->fdisp)->dy[l];
	(blur->fdisp)->dz[l] = (rcloud->fdisp)->dz[l];
	
	(blur->disp)->dx[l] = ROUND((blur->fdisp)->dx[l]);
	(blur->disp)->dy[l] = ROUND((blur->fdisp)->dy[l]);
	(blur->disp)->dz[l] = ROUND((blur->fdisp)->dz[l]);
	
	tmp = Scene32::AddFrame(rcloud->prob[l], 1, 0);
	blur->prob[l] = Scene32::FastGaussianBlur(tmp);
	Scene32::Destroy(&tmp);
      }
      return blur;
    }
    
    
    RegionCloud3 *ChangeOrientationToLPS(RegionCloud3 *rcloud,
					 char *ori){
      RegionCloud3 *lps=NULL;
      int nobjs,l;
      
      nobjs = rcloud->nobjs;
      
      //Create RegionCloud3:
      lps = (RegionCloud3 *)calloc(1,sizeof(RegionCloud3));
      if(lps==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::ChangeOrientationToLPS");
      lps->prob = (Scene32::Scene32 **)calloc(nobjs+1,sizeof(Scene32::Scene32 *));
      if(lps->prob==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::ChangeOrientationToLPS");
      lps->disp  = AdjRel3::Create(nobjs+1);
      lps->fdisp = AdjRel3f::ChangeOrientationToLPS(rcloud->fdisp, ori);
      lps->nobjs = rcloud->nobjs;
      lps->nimages = rcloud->nimages;
      
      for(l=1; l<=nobjs; l++){
	(lps->disp)->dx[l] = ROUND((lps->fdisp)->dx[l]);
	(lps->disp)->dy[l] = ROUND((lps->fdisp)->dy[l]);
	(lps->disp)->dz[l] = ROUND((lps->fdisp)->dz[l]);
	
	lps->prob[l] = Scene32::ChangeOrientationToLPS(rcloud->prob[l], ori);
      }
      return lps;
    }



    RegionCloud3 *Read(char *filename){
      RegionCloud3 *rcloud=NULL;
      char tmp[512],more[512];
      int nimages,nobjs,l;
      FILE *fp;
      
      fp = fopen(filename,"rb");
      if(fp == NULL)
	gft::Error((char *)MSG2,
		   (char *)"RegionCloud3::Read");
      
      fread(&nimages,sizeof(int),1,fp);
      fread(&nobjs,  sizeof(int),1,fp);
      
      //Create RegionCloud3:
      rcloud = (RegionCloud3 *)calloc(1,sizeof(RegionCloud3));
      if(rcloud==NULL)
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::Read");
      rcloud->prob = (Scene32::Scene32 **)calloc(nobjs+1,sizeof(Scene32::Scene32 *));
      if(rcloud->prob==NULL)
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::Read");
      rcloud->disp  = AdjRel3::Create(nobjs+1);
      rcloud->fdisp = AdjRel3f::Create(nobjs+1);
      rcloud->nobjs   = nobjs;
      rcloud->nimages = nimages;

      fread((rcloud->fdisp)->dx,sizeof(float),nobjs+1,fp);
      fread((rcloud->fdisp)->dy,sizeof(float),nobjs+1,fp);
      fread((rcloud->fdisp)->dz,sizeof(float),nobjs+1,fp);
      fclose(fp);
      
      for(l=1; l<=nobjs; l++){
	(rcloud->disp)->dx[l] = ROUND((rcloud->fdisp)->dx[l]);
	(rcloud->disp)->dy[l] = ROUND((rcloud->fdisp)->dy[l]);
	(rcloud->disp)->dz[l] = ROUND((rcloud->fdisp)->dz[l]);
	
	strcpy(tmp, filename);
	FileList::RemoveFileExtension(tmp);
	sprintf(more, "_%02d.scn.bz2", l);
	strcat(tmp, more);
	rcloud->prob[l] = Scene32::Read(tmp);
      }
      return rcloud;
    }
    

    void Write(RegionCloud3 *rcloud, 
	       char *filename){
      char tmp[512],more[512];
      int l;
      FILE *fp;
      
      fp = fopen(filename,"wb");
      if(fp == NULL)
	gft::Error((char *)MSG2,
		   (char *)"RegionCloud3::Write");
      
      fwrite(&rcloud->nimages,sizeof(int),1,fp);
      fwrite(&rcloud->nobjs,  sizeof(int),1,fp);
      fwrite((rcloud->fdisp)->dx,sizeof(float),rcloud->nobjs+1,fp);
      fwrite((rcloud->fdisp)->dy,sizeof(float),rcloud->nobjs+1,fp);
      fwrite((rcloud->fdisp)->dz,sizeof(float),rcloud->nobjs+1,fp);
      fclose(fp);
      
      for(l=1; l<=rcloud->nobjs; l++){
	strcpy(tmp, filename);
	FileList::RemoveFileExtension(tmp);
	sprintf(more, "_%02d.scn.bz2", l);
	strcat(tmp, more);
	Scene32::Write(rcloud->prob[l], tmp);
      }
    }
    

    void RemoveElem(RegionCloud3 *rcloud,
		    Scene32::Scene32 *label){
      Scene8::Scene8 *mask=NULL;
      gft::Voxel u,v;
      int p,q,l;
      float prob,sum;
      gft::Voxel *C = NULL;   //centroids.
      
      C = (gft::Voxel *) calloc(rcloud->nobjs+1, sizeof(gft::Voxel));
      if(C==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"RegionCloud3::RemoveElem");
      
      //Compute centroids:
      for(l=1; l<=rcloud->nobjs; l++){
	mask = Scene32::Threshold(label, l, l);
	C[l] = Scene8::Centroid(mask);
	Scene8::Destroy(&mask);
      }

      //Probability to Sum:
      for(l=1; l<=rcloud->nobjs; l++){
	for(p=0; p<(rcloud->prob[l])->n; p++){
	  prob  = (float)(rcloud->prob[l])->data[p];
	  prob /= (float)MAX_PROB;
	  (rcloud->prob[l])->data[p] = ROUND(prob*rcloud->nimages);
	}
      }

      //Remove label:
      for(p=0; p<label->n; p++){
	l = label->data[p];
	if(l==0) continue;
	u.x = Scene32::GetAddressX(label, p);
	u.y = Scene32::GetAddressY(label, p);
	u.z = Scene32::GetAddressZ(label, p);
	
	u.x -= C[l].x;
	u.y -= C[l].y;
	u.z -= C[l].z;
	
	v.x = u.x + (rcloud->prob[l])->xsize/2;
	v.y = u.y + (rcloud->prob[l])->ysize/2;
	v.z = u.z + (rcloud->prob[l])->zsize/2;

	if(Scene32::IsValidVoxel(rcloud->prob[l], v)){
	  q = Scene32::GetVoxelAddress(rcloud->prob[l], v);
	  (rcloud->prob[l])->data[q] -= 1;
	}
      }

      //Sum to Probability:
      rcloud->nimages--;
      for(l=1; l<=rcloud->nobjs; l++){
	for(p=0; p<(rcloud->prob[l])->n; p++){
	  sum = (float)(rcloud->prob[l])->data[p];
	  (rcloud->prob[l])->data[p] = ROUND(MAX_PROB*(sum/rcloud->nimages));
	}
      }
      free(C);
    }



  } //end RegionCloud3 namespace
} //end gft namespace
    //-----------------------------------------

namespace gft{
  namespace BorderCloud3{
    
    BorderCloud3 *ByRegionCloud(RegionCloud3::RegionCloud3 *rcloud){
      BorderCloud3 *bcloud=NULL;
      Scene32::Scene32 *prob=NULL;
      Gradient3::Gradient3 *grad=NULL;
      float dx,dy,dz,r;
      int rv,Imax,l,nobjs = rcloud->nobjs;
      
      //Create BorderCloud3:
      bcloud = (BorderCloud3 *)calloc(1,sizeof(BorderCloud3));
      if(bcloud==NULL)
	gft::Error((char *)MSG1,
		   (char *)"BorderCloud3::ByRegionCloud");
      bcloud->prob = (ScnGradient **)calloc(nobjs+1,sizeof(ScnGradient *));
      if(bcloud->prob==NULL)
	gft::Error((char *)MSG1,
		   (char *)"BorderCloud3::ByRegionCloud");
      bcloud->disp  = AdjRel3::Create(nobjs+1);
      bcloud->fdisp = AdjRel3f::Create(nobjs+1);
      bcloud->nobjs = nobjs;
      
      for(l=1; l<=rcloud->nobjs; l++){
	dx = (rcloud->prob[l])->dx;
	dy = (rcloud->prob[l])->dy;
	dz = (rcloud->prob[l])->dz;
	
	//----------------------------------------
	r = 3.0;
	rv = MAX(ROUND(r/dx), 1);
	rv = MAX(ROUND(r/dy), rv);
	rv = MAX(ROUND(r/dz), rv);
	prob = Scene32::AddFrame(rcloud->prob[l], rv*2, 0);
	grad = Gradient3::Spherical(prob, r);
	Imax = Gradient3::MaximumMag(grad);
	Gradient3::Normalize(grad, 0, Imax, 0, MAX_PROB);
	bcloud->prob[l] = Gradient3::RemFrame(grad, rv);
	//Gradient3::PowerEnhancement(bcloud->prob[l]);
	//----------------------------------------
	
	Scene32::Destroy(&prob);
	Gradient3::Destroy(&grad);
	(bcloud->disp)->dx[l] = (rcloud->disp)->dx[l];
	(bcloud->disp)->dy[l] = (rcloud->disp)->dy[l];
	(bcloud->disp)->dz[l] = (rcloud->disp)->dz[l];
	
	(bcloud->fdisp)->dx[l] = (rcloud->fdisp)->dx[l];
	(bcloud->fdisp)->dy[l] = (rcloud->fdisp)->dy[l];
	(bcloud->fdisp)->dz[l] = (rcloud->fdisp)->dz[l];
      }
      return bcloud;
    }

    
    void Destroy(BorderCloud3 **bcloud){
      BorderCloud3 *aux;
      int i;
      aux = *bcloud;
      if(aux != NULL){
	if(aux->disp!=NULL)
	  AdjRel3::Destroy(&aux->disp);
	if(aux->fdisp!=NULL)
	  AdjRel3f::Destroy(&aux->fdisp);
	
	if(aux->prob!=NULL){
	  for(i=1; i<=aux->nobjs; i++)
	    if(aux->prob[i]!=NULL)
	      Gradient3::Destroy(&aux->prob[i]);
	  free(aux->prob);
	}
	free(aux);
	*bcloud = NULL;
      }
    }


    void GetVoxelSize(BorderCloud3 *bcloud,
		      float *dx, 
		      float *dy, 
		      float *dz){
      *dx = ((bcloud->prob[1])->Gx)->dx;
      *dy = ((bcloud->prob[1])->Gx)->dy;
      *dz = ((bcloud->prob[1])->Gx)->dz;
    }



    BorderCloud3 *Subsampling(BorderCloud3 *bcloud){
      BorderCloud3 *sub=NULL;
      Scene32::Scene32 *tmp=NULL;
      int l;
      
      //Create BorderCloud3:
      sub = (BorderCloud3 *)calloc(1,sizeof(BorderCloud3));
      if(sub==NULL)
	gft::Error((char *)MSG1,
		   (char *)"BorderCloud3::Subsampling");
      sub->prob = (Gradient3 **)calloc(bcloud->nobjs+1,sizeof(Gradient3 *));
      if(sub->prob==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"BorderCloud3::Subsampling");
      sub->disp  = AdjRel3::Create(bcloud->nobjs+1);
      sub->fdisp = AdjRel3f::Create(bcloud->nobjs+1);
      sub->nobjs = bcloud->nobjs;
      
      for(l=1; l<=bcloud->nobjs; l++){
	(sub->fdisp)->dx[l] = (bcloud->fdisp)->dx[l]/2.;
	(sub->fdisp)->dy[l] = (bcloud->fdisp)->dy[l]/2.;
	(sub->fdisp)->dz[l] = (bcloud->fdisp)->dz[l]/2.;
	
	(sub->disp)->dx[l] = ROUND((sub->fdisp)->dx[l]);
	(sub->disp)->dy[l] = ROUND((sub->fdisp)->dy[l]);
	(sub->disp)->dz[l] = ROUND((sub->fdisp)->dz[l]);
	
	sub->prob[l] = (Gradient3 *)calloc(1,sizeof(Gradient3));
	if(sub->prob[l] == NULL)
	  gft::Error((char *)MSG1,
		     (char *)"BorderCloud3::Subsampling");
	(sub->prob[l])->mag = NULL;
	
	tmp  = Scene32::FastGaussianBlur((bcloud->prob[l])->Gx);
	(sub->prob[l])->Gx = Scene32::Subsampling(tmp);
	Scene32::Destroy(&tmp);

	tmp  = Scene32::FastGaussianBlur((bcloud->prob[l])->Gy);
	(sub->prob[l])->Gy = Scene32::Subsampling(tmp);
	Scene32::Destroy(&tmp);
	
	tmp  = Scene32::FastGaussianBlur((bcloud->prob[l])->Gz);
	(sub->prob[l])->Gz = Scene32::Subsampling(tmp);
	Scene32::Destroy(&tmp);
      }
      return sub;
    }


    BorderCloud3 *LinearInterp(BorderCloud3 *bcloud,
			       float dx,float dy,float dz){
      BorderCloud3 *interp=NULL;
      float odx,ody,odz;
      int nobjs,l;
      
      nobjs = bcloud->nobjs;
      
      //Create BorderCloud3:
      interp = (BorderCloud3 *)calloc(1,sizeof(BorderCloud3));
      if(interp==NULL)
	gft::Error((char *)MSG1,
		   (char *)"BorderCloud3::LinearInterp");
      interp->prob = (Gradient3::Gradient3 **)calloc(nobjs+1,
						     sizeof(Gradient3::Gradient3 *));
      if(interp->prob==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"BorderCloud3::LinearInterp");
      interp->disp  = AdjRel3::Create(nobjs+1);
      interp->fdisp = AdjRel3f::Create(nobjs+1);
      interp->nobjs = bcloud->nobjs;
      for(l=1; l<=nobjs; l++){
	odx = ((bcloud->prob[l])->Gx)->dx;
	ody = ((bcloud->prob[l])->Gx)->dy;
	odz = ((bcloud->prob[l])->Gx)->dz;
	(interp->fdisp)->dx[l] = (bcloud->fdisp)->dx[l]*odx/dx;
	(interp->fdisp)->dy[l] = (bcloud->fdisp)->dy[l]*ody/dy;
	(interp->fdisp)->dz[l] = (bcloud->fdisp)->dz[l]*odz/dz;

	(interp->disp)->dx[l] = ROUND((interp->fdisp)->dx[l]);
	(interp->disp)->dy[l] = ROUND((interp->fdisp)->dy[l]);
	(interp->disp)->dz[l] = ROUND((interp->fdisp)->dz[l]);

	interp->prob[l] = Gradient3::LinearInterpCentr(bcloud->prob[l], 
						       dx,dy,dz);
      }
      return interp;
    }
    
    
    BorderCloud3 *ChangeOrientationToLPS(BorderCloud3 *bcloud,
					 char *ori){
      BorderCloud3 *lps=NULL;
      int nobjs,l;
      nobjs = bcloud->nobjs;
      //Create BorderCloud3:
      lps = (BorderCloud3 *)calloc(1,sizeof(BorderCloud3));
      if(lps==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"BorderCloud3::ChangeOrientationToLPS");
      lps->prob = (Gradient3 **)calloc(nobjs+1,sizeof(Gradient3 *));
      if(lps->prob==NULL) 
	gft::Error((char *)MSG1,
		   (char *)"BorderCloud3::ChangeOrientationToLPS");
      lps->disp  = AdjRel3::Create(nobjs+1);
      lps->fdisp = AdjRel3f::ChangeOrientationToLPS(bcloud->fdisp, ori);
      lps->nobjs = bcloud->nobjs;
      
      for(l=1; l<=nobjs; l++){
	(lps->disp)->dx[l] = ROUND((lps->fdisp)->dx[l]);
	(lps->disp)->dy[l] = ROUND((lps->fdisp)->dy[l]);
	(lps->disp)->dz[l] = ROUND((lps->fdisp)->dz[l]);
	
	lps->prob[l] = Gradient3::ChangeOrientationToLPS(bcloud->prob[l], ori);
      }
      return lps;
    }
    

    void Normalize(BorderCloud3 *bcloud,
		   int omin,int omax,
		   int nmin,int nmax){
      int nobjs,l;
      nobjs = bcloud->nobjs;
      for(l=1; l<=nobjs; l++){
	Gradient3::Normalize(bcloud->prob[l],
			     omin,omax,
			     nmin,nmax);
      }
    }
    
    
    BorderCloud3 *Read(char *filename){
      BorderCloud3 *bcloud=NULL;
      char tmp[512],more[512];
      int nobjs,l;
      FILE *fp;
      
      fp = fopen(filename,"rb");
      if(fp == NULL)
	gft::Error((char *)MSG2,
		   (char *)"BorderCloud3::Read");
      
      fread(&nobjs,  sizeof(int),1,fp);
      
      //Create BorderCloud3:
      bcloud = (BorderCloud3 *)calloc(1,sizeof(BorderCloud3));
      if(bcloud==NULL)
	gft::Error((char *)MSG1,
		   (char *)"BorderCloud3::Read");
      bcloud->prob = (Gradient3 **)calloc(nobjs+1,sizeof(Gradient3 *));
      if(bcloud->prob==NULL)
	gft::Error((char *)MSG1,
		   (char *)"BorderCloud3::Read");
      bcloud->disp  = AdjRel3::Create(nobjs+1);
      bcloud->fdisp = AdjRel3f::Create(nobjs+1);
      bcloud->nobjs = nobjs;
      
      fread((bcloud->fdisp)->dx,sizeof(float),nobjs+1,fp);
      fread((bcloud->fdisp)->dy,sizeof(float),nobjs+1,fp);
      fread((bcloud->fdisp)->dz,sizeof(float),nobjs+1,fp);
      fclose(fp);
      
      for(l=1; l<=nobjs; l++){
	(bcloud->disp)->dx[l] = ROUND((bcloud->fdisp)->dx[l]);
	(bcloud->disp)->dy[l] = ROUND((bcloud->fdisp)->dy[l]);
	(bcloud->disp)->dz[l] = ROUND((bcloud->fdisp)->dz[l]);
	
	strcpy(tmp, filename);
	FileList::RemoveFileExtension(tmp);
	sprintf(more, "_%02d.sgr.bz2", l);
	strcat(tmp, more);
	bcloud->prob[l] = Gradient3::ReadCompressed(tmp);
      }
      return bcloud;
    }


    void Write(BorderCloud3 *bcloud, 
	       char *filename){
      char tmp[512],more[512];
      int l;
      FILE *fp;
      fp = fopen(filename,"wb");
      if(fp == NULL)
	gft::Error((char *)MSG2,
		   (char *)"BorderCloud3::Write");
      fwrite(&bcloud->nobjs,  sizeof(int),1,fp);
      fwrite((bcloud->fdisp)->dx,sizeof(float),bcloud->nobjs+1,fp);
      fwrite((bcloud->fdisp)->dy,sizeof(float),bcloud->nobjs+1,fp);
      fwrite((bcloud->fdisp)->dz,sizeof(float),bcloud->nobjs+1,fp);
      fclose(fp);

      for(l=1; l<=bcloud->nobjs; l++){
	strcpy(tmp, filename);
	FileList::RemoveFileExtension(tmp);
	sprintf(more, "_%02d.sgr.bz2", l);
	strcat(tmp, more);
	Gradient3::WriteCompressed(bcloud->prob[l], tmp);
      }
    }



  } //end BorderCloud3 namespace
} //end gft namespace
