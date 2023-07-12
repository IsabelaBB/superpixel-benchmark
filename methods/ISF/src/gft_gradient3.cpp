
#include "gft_gradient3.h"


namespace gft{
  namespace Gradient3{

    Gradient3 *Create(int xsize,int ysize,int zsize){
      Gradient3 *grad=NULL;
      grad = (Gradient3 *) calloc(1,sizeof(Gradient3));
      if(grad == NULL)
	gft::Error((char *)MSG1,
		   (char *)"Gradient3::Create");
      grad->Gx = Gradient3::Create(xsize,ysize,zsize);
      grad->Gy = Gradient3::Create(xsize,ysize,zsize);
      grad->Gz = Gradient3::Create(xsize,ysize,zsize);
      grad->mag = NULL;
      return grad;
    }


    void Destroy(Gradient3 **grad){
      Gradient3 *aux;
      aux = *grad;
      if(aux != NULL){
	if(aux->Gx != NULL) Scene32::Destroy(&aux->Gx);
	if(aux->Gy != NULL) Scene32::Destroy(&aux->Gy);
	if(aux->Gz != NULL) Scene32::Destroy(&aux->Gz);
	if(aux->mag != NULL) Scene32::Destroy(&aux->mag);
	free(aux);    
	*grad = NULL;
      }
    }


    Gradient3 *RemFrame(Gradient3 *fgrad, int sz){
      Gradient3 *grad=NULL;
      grad = (Gradient3 *) calloc(1,sizeof(Gradient3));
      if(grad == NULL)
	gft::Error((char *)MSG1,
		   (char *)"Gradient3::RemFrame");
      grad->Gx = Gradient3::RemFrame(fgrad->Gx, sz);
      grad->Gy = Gradient3::RemFrame(fgrad->Gy, sz);
      grad->Gz = Gradient3::RemFrame(fgrad->Gz, sz);
      if(fgrad->mag!=NULL){
	grad->mag = Gradient3::RemFrame(fgrad->mag, sz);
	Scene32::GetMaximumValue(grad->mag);
      }
      else
	grad->mag = NULL;
      
      return grad;
    }


    Gradient3 *LinearInterpCentr(Gradient3 *grad,
				 float dx,float dy,float dz){
      Gradient3 *interp=NULL;
      interp = (Gradient3 *) calloc(1,sizeof(Gradient3));
      if(interp == NULL)
	gft::Error((char *)MSG1,
		   (char *)"Gradient3::LinearInterpCentr");
      interp->Gx = FastLinearInterpCentr3(grad->Gx, dx,dy,dz);
      interp->Gy = FastLinearInterpCentr3(grad->Gy, dx,dy,dz);
      interp->Gz = FastLinearInterpCentr3(grad->Gz, dx,dy,dz);
      interp->mag = NULL;
      return interp;
    }


    Gradient3 *ChangeOrientationToLPS(Gradient3 *grad,
				      char *ori){
      Gradient3 *lps=NULL;
      Scene32::Scene32 *gx,*gy,*gz;
      
      lps = (Gradient3 *) calloc(1,sizeof(Gradient3));
      if(lps == NULL)
	gft::Error((char *)MSG1,
		   (char *)"Gradient3::ChangeOrientationToLPS");
      gx = Scene32::ChangeOrientationToLPS(grad->Gx, ori);
      gy = Scene32::ChangeOrientationToLPS(grad->Gy, ori);
      gz = Scene32::ChangeOrientationToLPS(grad->Gz, ori);
      lps->Gx  = NULL;
      lps->Gy  = NULL;
      lps->Gz  = NULL;
      lps->mag = NULL;
      
      if     (ori[0]=='L'){ lps->Gx = gx; }
      else if(ori[0]=='R'){ lps->Gx = gx; Scene32::Negate3inplace(gx); }
      else if(ori[0]=='P'){ lps->Gy = gx; }
      else if(ori[0]=='A'){ lps->Gy = gx; Scene32::Negate3inplace(gx); }
      else if(ori[0]=='S'){ lps->Gz = gx; }
      else if(ori[0]=='I'){ lps->Gz = gx; Scene32::Negate3inplace(gx); }
      else{ gft::Error((char *)"Invalid orientation",
		       (char *)"Gradient3::ChangeOrientationToLPS"); }
      
      if     (ori[1]=='L'){ lps->Gx = gy; }
      else if(ori[1]=='R'){ lps->Gx = gy; Scene32::Negate3inplace(gy); }
      else if(ori[1]=='P'){ lps->Gy = gy; }
      else if(ori[1]=='A'){ lps->Gy = gy; Scene32::Negate3inplace(gy); }
      else if(ori[1]=='S'){ lps->Gz = gy; }
      else if(ori[1]=='I'){ lps->Gz = gy; Scene32::Negate3inplace(gy); }
      else{ gft::Error((char *)"Invalid orientation",
		       (char *)"Gradient3::ChangeOrientationToLPS"); }
      
      if     (ori[2]=='L'){ lps->Gx = gz; }
      else if(ori[2]=='R'){ lps->Gx = gz; Scene32::Negate3inplace(gz); }
      else if(ori[2]=='P'){ lps->Gy = gz; }
      else if(ori[2]=='A'){ lps->Gy = gz; Scene32::Negate3inplace(gz); }
      else if(ori[2]=='S'){ lps->Gz = gz; }
      else if(ori[2]=='I'){ lps->Gz = gz; Scene32::Negate3inplace(gz); }
      else{ gft::Error((char *)"Invalid orientation",
		       (char *)"Gradient3::ChangeOrientationToLPS"); }
      
      if(lps->Gx==NULL || lps->Gy==NULL || lps->Gz==NULL)
	gft::Error((char *)"Invalid orientation",
		   (char *)"Gradient3::ChangeOrientationToLPS");
      
      return lps;
    }



    Gradient3 *Read(char *filename){
      Gradient3 *grad=NULL;
      int xsize,ysize,zsize,n;
      float dx,dy,dz;
      char msg[512];
      FILE *fp;
      fp = fopen(filename,"rb");
      if(fp == NULL){
	sprintf(msg,"Cannot open %s",filename);
	gft::Error((char *)msg,
		   (char *)"Gradient3::Read");
      }
      fread(&xsize, sizeof(int), 1, fp);
      fread(&ysize, sizeof(int), 1, fp);
      fread(&zsize, sizeof(int), 1, fp);
      grad = Gradient3::Create(xsize,ysize,zsize);
      fread(&dx, sizeof(float), 1, fp);
      fread(&dy, sizeof(float), 1, fp);
      fread(&dz, sizeof(float), 1, fp);
      grad->Gx->dx = dx;
      grad->Gx->dy = dy;
      grad->Gx->dz = dz;

      grad->Gy->dx = dx;
      grad->Gy->dy = dy;
      grad->Gy->dz = dz;

      grad->Gz->dx = dx;
      grad->Gz->dy = dy;
      grad->Gz->dz = dz;

      n = xsize*ysize*zsize;
      fread((grad->Gx)->data, sizeof(int), n, fp);
      fread((grad->Gy)->data, sizeof(int), n, fp);
      fread((grad->Gz)->data, sizeof(int), n, fp);
      fclose(fp);
      return grad;
    }
    

    void Write(Gradient3 *grad, char *filename){
      int xsize,ysize,zsize,n;
      float dx,dy,dz;
      char msg[512];
      FILE *fp;
      
      fp = fopen(filename,"wb");
      if(fp == NULL){
	sprintf(msg,"Cannot open %s",filename);
	gft::Error((char *)msg,
		   (char *)"Gradient3::Write");
      }
      xsize = (grad->Gx)->xsize;
      ysize = (grad->Gx)->ysize;
      zsize = (grad->Gx)->zsize;
      fwrite(&xsize, sizeof(int),  1, fp);
      fwrite(&ysize, sizeof(int),  1, fp);
      fwrite(&zsize, sizeof(int),  1, fp);
      dx = (grad->Gx)->dx;
      dy = (grad->Gx)->dy;
      dz = (grad->Gx)->dz;
      fwrite(&dx, sizeof(float),  1, fp);
      fwrite(&dy, sizeof(float),  1, fp);
      fwrite(&dz, sizeof(float),  1, fp);
      n = xsize*ysize*zsize;
      fwrite((grad->Gx)->data, sizeof(int), n, fp);
      fwrite((grad->Gy)->data, sizeof(int), n, fp);
      fwrite((grad->Gz)->data, sizeof(int), n, fp);
      fclose(fp);
    }



    Gradient3 *Spherical(Scene32::Scene32 *scn, float r){
      AdjRel3::AdjRel3 *A = AdjRel3::Spherical_mm(scn, r);
      AdjRel3::AdjVxl  *N;
      Gradient3 *grad=NULL;
      int p,q,i,vp,vq;
      float *mg=NULL;
      float gx,gy,gz,d;
      int fx,fy,fz,dp;
      
      N  = AdjRel3::AdjVoxels(scn, A);
      mg = AdjRel3::GetDistanceArray_mm(A, scn);
      AdjRel3::GetFrameSize(A, &fx, &fy, &fz);
      dp  = fx*1;
      dp += fy*scn->xsize;
      dp += fz*scn->xsize*scn->ysize;
      
      grad = Gradient3::Create(scn->xsize,
			       scn->ysize,
			       scn->zsize);
      grad->Gx->dx = scn->dx;
      grad->Gx->dy = scn->dy;
      grad->Gx->dz = scn->dz;

      grad->Gy->dx = scn->dx;
      grad->Gy->dy = scn->dy;
      grad->Gy->dz = scn->dz;
      
      grad->Gz->dx = scn->dx;
      grad->Gz->dy = scn->dy;
      grad->Gz->dz = scn->dz;
      
      for(p=dp; p<scn->n-dp; p++){
	vp = scn->data[p];
	gx = gy = gz = 0.0;
	for(i=1; i<N->n; i++){
	  q = p + N->dp[i];
	  vq = scn->data[q];
	  d = (float)(vq-vp);
	
	  gx  += (d*A->dx[i])/(mg[i]);
	  gy  += (d*A->dy[i])/(mg[i]);
	  gz  += (d*A->dz[i])/(mg[i]);
	}
	(grad->Gx)->data[p] = ROUND(10.0*gx*scn->dx/fx);
	(grad->Gy)->data[p] = ROUND(10.0*gy*scn->dy/fy);
	(grad->Gz)->data[p] = ROUND(10.0*gz*scn->dz/fz);
      }
      Scene32::ClearAdjFrame(grad->Gx, A);
      Scene32::ClearAdjFrame(grad->Gy, A);
      Scene32::ClearAdjFrame(grad->Gz, A);
      AdjRel3::Destroy(&A);
      AdjRel3::DestroyAdjVxl(&N);
      free(mg);
      return grad;
    }



    void ComputeMagnitude(Gradient3 *grad){
      Scene32::Scene32 *mag=NULL;
      float gx,gy,gz,d;
      int p,m,Gmax=0;
      
      if(grad->mag!=NULL) return;
      
      mag = Scene32::Create((grad->Gx)->xsize,
			    (grad->Gx)->ysize,
			    (grad->Gx)->zsize);
      mag->dx = grad->Gx->dx;
      mag->dy = grad->Gx->dy;
      mag->dz = grad->Gx->dz;
      for(p=0; p<mag->n; p++){
	gx = (float)(grad->Gx->data[p]);
	gy = (float)(grad->Gy->data[p]);
	gz = (float)(grad->Gz->data[p]);
	d = sqrtf(gx*gx + gy*gy + gz*gz);
	m = ROUND(d);
	if(m>Gmax) Gmax = m;
	mag->data[p] = m;
      }
      mag->maxval = Gmax;
      grad->mag = mag;
    }


    int MaximumMag(Gradient3 *grad){
      int Imax;
      if(grad->mag==NULL)
	ComputeMagnitude(grad);
      Imax = Scene32::GetMaximumValue(grad->mag);
      return Imax;
    }


    void Normalize(Gradient3 *grad,
		   int omin,int omax,
		   int nmin,int nmax){
      float gx,gy,gz,m,nm;
      int p;
      ComputeMagnitude(grad);
      for(p=0; p<(grad->Gx)->n; p++){
	m  = (grad->mag)->data[p];
	nm = (float)IntegerNormalize((grad->mag)->data[p],
				     omin, omax, 
				     nmin, nmax);
	gx = grad->Gx->data[p];
	gy = grad->Gy->data[p];
	gz = grad->Gz->data[p];
	if((grad->mag)->data[p]>0){
	  (grad->Gx)->data[p] = ROUND(gx*(nm/m));
	  (grad->Gy)->data[p] = ROUND(gy*(nm/m));
	  (grad->Gz)->data[p] = ROUND(gz*(nm/m));
	  (grad->mag)->data[p] = ROUND(nm);
	}
	else{
	  (grad->Gx)->data[p] = 0;
	  (grad->Gy)->data[p] = 0;
	  (grad->Gz)->data[p] = 0;
	  (grad->mag)->data[p] = 0;
	}
      }
      (grad->mag)->maxval = IntegerNormalize((grad->mag)->maxval,
					     omin, omax, 
					     nmin, nmax);
    }


    void PowerEnhancement(Gradient3 *grad){
      Scene32::Scene32 *pmag;
      float gx,gy,gz,m,pm;
      int p;
      ComputeMagnitude(grad);
      pmag = Scene32::Clone(grad->mag);
      Scene32::PowerEnhancement(pmag);
      for(p=0; p<pmag->n; p++){
	gx = (grad->Gx)->data[p];
	gy = (grad->Gy)->data[p];
	gz = (grad->Gz)->data[p];
	m  = (grad->mag)->data[p];
	pm = pmag->data[p];
	if((grad->mag)->data[p]>0){
	  grad->Gx->data[p] = ROUND((gx/m)*pm);
	  grad->Gy->data[p] = ROUND((gy/m)*pm);
	  grad->Gz->data[p] = ROUND((gz/m)*pm);
	  grad->mag->data[p] = ROUND(pm);
	}
	else{
	  grad->Gx->data[p] = 0;
	  grad->Gy->data[p] = 0;
	  grad->Gz->data[p] = 0;
	  grad->mag->data[p] = 0;
	}
      }
      Scene32::Destroy(&pmag);
    }
   
    

  } //end Gradient3 namespace
} //end gft namespace

