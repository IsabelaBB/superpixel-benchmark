
#include "gft_preproc3.h"


namespace gft{
  namespace Kernel3{

    Kernel3 *Create(AdjRel3::AdjRel3 *A){
      Kernel3 *K=NULL;
      int i;
      gft::Voxel max, min;
      
      max.x = max.y = max.z = INT_MIN;
      min.x = min.y = min.z = INT_MAX;
      
      K = (Kernel3 *) calloc(1,sizeof(Kernel3));
      if (K == NULL){
	gft::Error((char *)MSG1,
		   (char *)"Kernel3::Create");
      }
      K->val = gft::AllocFloatArray(A->n);
      K->adj = AdjRel3::Create(A->n);
      
      for (i=0;i<A->n;i++) {
	max.x = MAX(A->dx[i],max.x);
	max.y = MAX(A->dy[i],max.y);
	max.z = MAX(A->dz[i],max.z);
	min.x = MIN(A->dx[i],min.x);
	min.y = MIN(A->dy[i],min.y);
	min.z = MIN(A->dz[i],min.z);
	K->adj->dx[i] = A->dx[i];
	K->adj->dy[i] = A->dy[i];
	K->adj->dz[i] = A->dz[i];
      }
      K->xsize = max.x - min.x + 1;
      K->ysize = max.y - min.y + 1;
      K->zsize = max.z - min.z + 1;
      return(K);
    }


    Kernel3 *Clone(Kernel3 *K){
      Kernel3 *C;
      
      C = (Kernel3 *) calloc(1,sizeof(Kernel3));
      if(C == NULL)
	gft::Error((char *)MSG1,
		   (char *)"Kernel3::Clone");
      
      C->val = gft::AllocFloatArray(K->adj->n);
      memcpy(C->val, K->val,
	     sizeof(float)*K->adj->n);
      C->adj   = AdjRel3::Clone(K->adj);
      C->xsize = K->xsize;
      C->ysize = K->ysize;
      C->zsize = K->zsize;
      return C;
    }

    void Destroy(Kernel3 **K){
      Kernel3 *aux;
      aux = *K;
      if(aux != NULL){
	if (aux->val != NULL) gft::FreeFloatArray(&aux->val); 
	AdjRel3::Destroy(&(aux->adj));
	free(aux);    
	*K = NULL;
      }
    }


    Kernel3 *Normalize(Kernel3 *K){
      Kernel3 *C = Clone(K);
      float wt=0.0;
      int i;
      
      for(i=0; i<K->adj->n; i++)
	wt += K->val[i];
      for(i=0; i<C->adj->n; i++)
	C->val[i] /= wt;
      
      return C;
    }


    Kernel3 *SphericalGaussian(float R, float s, float f){
      float R2,r2;
      Kernel3 *K;
      AdjRel3::AdjRel3 *A;
      int i;
      R2 = R*R;
      A = AdjRel3::Spheric(R);
      K = Create(A);
      for (i=0;i<A->n;i++) {
	r2 = A->dx[i] * A->dx[i] + A->dy[i] * A->dy[i] + A->dz[i] * A->dz[i];
	K->val[i] = s*exp (-f*(r2/R2));
      }
      AdjRel3::Destroy(&A);
      return(K);
    }


  } //end Kernel3 namespace
} //end gft namespace


namespace gft{
  namespace Scene32{

    int    AreaPercentageHigherThreshold(Curve::Curve *hist,
					 float perc){
      Curve::Curve *nhist;
      double sum=0.0,ratio;
      int i;
      ratio = 1.0 - perc/100.0;
      nhist = Curve::Normalize(hist);
      i = nhist->n-1;
      while(sum<ratio && i>=0){
	sum += nhist->Y[i];
	i--;
      }
      i = MIN(nhist->n-1, i+1);
      Curve::Destroy(&nhist);
      return i;
    }


    void SuppressHighIntensities(Scene32 *scn){
      Curve::Curve *hist;
      int T,p,n;
      hist = Histogram(scn, 1);
      T = AreaPercentageHigherThreshold(hist, 99.95);
      n = scn->xsize*scn->ysize*scn->zsize;
      for(p=0; p<n; p++)
	if(scn->data[p]>T)
	  scn->data[p] = T;
      Curve::Destroy(&hist);
    }


    Scene32 *Convolution(Scene32 *scn, Kernel3::Kernel3 *K){
      Scene32 *cscn;
      gft::Voxel u,v;
      AdjRel3::AdjVxl *vxl;
      int p,q,i;
      float conv;
      
      cscn = Scene32::Create(scn->xsize,scn->ysize,scn->zsize);
      cscn->dx = scn->dx;
      cscn->dy = scn->dy;
      cscn->dz = scn->dz;
      vxl = AdjRel3::AdjVoxels(K->adj, scn);
      
      for(u.z=0; u.z<scn->zsize; u.z++)
	for(u.y=0; u.y<scn->ysize; u.y++)
	  for(u.x=0; u.x<scn->xsize; u.x++){
	    p = GetVoxelAddress(cscn,u.x,u.y,u.z);
	    conv=0.0;
	    for (i=0;i<K->adj->n;i++) {
	      v.x = u.x + K->adj->dx[i]; 
	      v.y = u.y + K->adj->dy[i];
	      v.z = u.z + K->adj->dz[i];
	      if (Scene32::IsValidVoxel(scn, v)){
		q = p + vxl->dp[i];
		conv += ((float)scn->data[q])*K->val[i];	   
	      }      
	    }
	    cscn->data[p] = ROUND(conv);
	  }
      
      AdjRel3::DestroyAdjVxl(&vxl);
      return(cscn);
    }


    Scene32 *OptConvolution3(Scene32 *scn, Kernel3::Kernel3 *K){
      Scene32 *cscn;
      AdjRel3::AdjVxl *vxl;
      int p,q,i,N,dx,dy,dz,dn;
      float conv;
      N = scn->xsize*scn->ysize*scn->zsize;
      cscn = Scene32::Create(scn->xsize,scn->ysize,scn->zsize);
      cscn->dx = scn->dx;
      cscn->dy = scn->dy;
      cscn->dz = scn->dz;
      vxl = AdjRel3::AdjVoxels(K->adj, scn);
      dx = dy = dz = 0;
      for(i=0;i<K->adj->n;i++){
	dx = MAX(dx, abs(K->adj->dx[i]));
	dy = MAX(dy, abs(K->adj->dy[i]));
	dz = MAX(dz, abs(K->adj->dz[i]));
      }
      dx *= 1;
      dy *= scn->xsize;
      dz *= scn->xsize * scn->ysize;
      dn = dx+dy+dz;
      
      for(p=dn; p<N-dn; p++){
	conv=0.0;
	for(i=0; i<K->adj->n; i++){
	  q = p + vxl->dp[i];
	  conv += ((float)scn->data[q])*K->val[i];	   
	}
	cscn->data[p] = ROUND(conv);
      }
      AdjRel3::DestroyAdjVxl(&vxl);
      return(cscn);
    }


    Scene32   *GaussianBlur(Scene32 *scn){
      Scene32 *blur;
      Kernel3::Kernel3 *K,*NK;
      K  = Kernel3::SphericalGaussian(1.0, 10.0, 1.0); //2.0, 10.0, 1.0
      NK = Kernel3::Normalize(K);
      blur = Convolution(scn, NK);
      blur->dx = scn->dx;
      blur->dy = scn->dy;
      blur->dz = scn->dz;
      Kernel3::Destroy(&NK);
      Kernel3::Destroy(&K);
      return blur;
    }

    Scene32   *FastGaussianBlur(Scene32 *scn){
      Scene32 *blur;
      Kernel3::Kernel3 *K,*NK;
      K  = Kernel3::SphericalGaussian(1.0, 10.0, 1.0); //2.0, 10.0, 1.0
      NK = Kernel3::Normalize(K);
      blur = FastConvolution(scn, NK);
      blur->dx = scn->dx;
      blur->dy = scn->dy;
      blur->dz = scn->dz;
      Kernel3::Destroy(&NK);
      Kernel3::Destroy(&K);
      return blur;
    }
    
    Scene32   *OptGaussianBlur(Scene32 *scn){
      Scene32 *blur;
      Kernel3::Kernel3 *K,*NK;
      K  = Kernel3::SphericalGaussian(1.0, 10.0, 1.0); //2.0, 10.0, 1.0
      NK = Kernel3::Normalize(K);
      blur = OptConvolution(scn, NK);
      blur->dx = scn->dx;
      blur->dy = scn->dy;
      blur->dz = scn->dz;
      Kernel3::Destroy(&NK);
      Kernel3::Destroy(&K);
      return blur;
    }


    Scene32   *FastOptGaussianBlur(Scene32 *scn){
      Scene32 *blur;
      Kernel3::Kernel3 *K,*NK;
      K  = Kernel3::SphericalGaussian(1.0, 10.0, 1.0); //2.0, 10.0, 1.0
      NK = Kernel3::Normalize(K);
      blur = FastOptConvolution(scn, NK);
      blur->dx = scn->dx;
      blur->dy = scn->dy;
      blur->dz = scn->dz;
      Kernel3::Destroy(&NK);
      Kernel3::Destroy(&K);
      return blur;
    }


    Scene32   *Subsampling(Scene32 *scn){
      Scene32 *scl;
      gft::Voxel u,v;
      int p,q;
      int xsize,ysize,zsize;
      xsize = (int)(scn->xsize/2) + (scn->xsize%2);
      ysize = (int)(scn->ysize/2) + (scn->ysize%2);
      zsize = (int)(scn->zsize/2) + (scn->zsize%2);
      if(xsize%2==0) xsize++;
      if(ysize%2==0) ysize++;
      if(zsize%2==0) zsize++;
      scl = Scene32::Create(xsize, ysize, zsize);
      scl->dx = 2.0*scn->dx;
      scl->dy = 2.0*scn->dy;
      scl->dz = 2.0*scn->dz;
      
      for(v.z=0; v.z<scl->zsize; v.z++)
	for(v.y=0; v.y<scl->ysize; v.y++)
	  for(v.x=0; v.x<scl->xsize; v.x++){
	    p = Scene32::VoxelAddress(scl, v);
	    u.x = (v.x-scl->xsize/2)*2 + scn->xsize/2;
	    u.y = (v.y-scl->ysize/2)*2 + scn->ysize/2;
	    u.z = (v.z-scl->zsize/2)*2 + scn->zsize/2;
	    
	    if(Scene32::IsValidVoxel(scn, u)){
	      q = Scene32::VoxelAddress(scn, u);
	      scl->data[p] = scn->data[q];
	    }
	  }
      
      return(scl);
    }


    Scene32 *LaplacianFilter(Scene32 *orig){
      Scene32 *lap;
      int dx,dy,dz,di,m,i,N;
      int acc;
      static int n0[9]={ -1,  0,  1, -1,  0,  1, -1,  0,  1};
      static int n1[9]={ -1, -1, -1,  0,  0,  0,  1,  1,  1};
      static int wt[9]={ -1, -1, -1, -1, -1, -1, -1, -1, -1};
      static int wc[9]={ -1, -1, -1, -1, 26, -1, -1, -1, -1};
      
      lap = Create(orig);
      N = orig->xsize*orig->ysize*orig->zsize;
      dx = 1;
      dy = orig->xsize;
      dz = orig->xsize * orig->ysize;
      di = dx+dy+dz;
      for(i=di; i<N-di; i++){
	// xy plane
	acc = 0;
	for(m=0;m<9;m++){
	  acc += wt[m]*orig->data[ i+dz+dx*n0[m]+dy*n1[m] ];
	  acc += wc[m]*orig->data[ i+ 0+dx*n0[m]+dy*n1[m] ];
	  acc += wt[m]*orig->data[ i-dz+dx*n0[m]+dy*n1[m] ];
	}
	lap->data[i] = acc;
      }
      return lap;
    }
    
    
    
    Scene32 *SobelFilter(Scene32 *scn){
      Scene32 *grad;
      int   dx,dy,dz,di,m,i,N;
      int   acc[3];
      float facc[3];
      static int n0[9]={-1, 0, 1,-1, 0, 1,-1, 0, 1};
      static int n1[9]={-1,-1,-1, 0, 0, 0, 1, 1, 1};
      static int we[9]={ 1, 2, 1, 2, 4, 2, 1, 2, 1};
      
      grad = Create(scn);
      N = scn->xsize*scn->ysize*scn->zsize;
      dx = 1;
      dy = scn->xsize;
      dz = scn->xsize * scn->ysize;
      di = dx+dy+dz;
      
      for(i=di; i<N-di; i++){
	// xy plane
	acc[0] = 0;
	for(m=0; m<9; m++)
	  acc[0] += we[m] * (scn->data[ i+dz+dx*n0[m]+dy*n1[m] ] -
			     scn->data[ i-dz+dx*n0[m]+dy*n1[m] ]);
	// yz plane
	acc[1] = 0;
	for(m=0; m<9; m++)
	  acc[1] += we[m] * (scn->data[ i+dx+dz*n0[m]+dy*n1[m] ] -
			     scn->data[ i-dx+dz*n0[m]+dy*n1[m] ]);
	// xz plane
	acc[2] = 0;
	for(m=0; m<9; m++)
	  acc[2] += we[m] * (scn->data[ i+dy+dx*n0[m]+dz*n1[m] ] -
			     scn->data[ i-dy+dx*n0[m]+dz*n1[m] ]);
	acc[0] >>= 4;
	acc[1] >>= 4;
	acc[2] >>= 4;
	facc[0] = (float)acc[0];
	facc[1] = (float)acc[1];
	facc[2] = (float)acc[2];
	grad->data[i] = ROUND(sqrtf((facc[0]*facc[0] + 
				     facc[1]*facc[1] +
				     facc[2]*facc[2])/3.0));
      }
      return grad;
    }
    

    Scene32 *SphericalGradient(Scene32 *scn, float r){
      AdjRel3::AdjRel3 *A = AdjRel3::Spherical_mm(scn, r);
      AdjRel3::AdjVxl  *N;
      Scene32 *grad=NULL;
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
      grad = Create(scn);
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
	gx = ROUND(10.0*gx*scn->dx/fx);
	gy = ROUND(10.0*gy*scn->dy/fy);
	gz = ROUND(10.0*gz*scn->dz/fz);
	grad->data[p] = ROUND(sqrtf(gx*gx + gy*gy + gz*gz));
      }
      Scene32::ClearAdjFrame(grad, A);
      AdjRel3::Destroy(&A);
      AdjRel3::DestroyAdjVxl(&N);
      free(mg);
      return grad;
    }



    
  } //end Scene32 namespace
} //end gft namespace

