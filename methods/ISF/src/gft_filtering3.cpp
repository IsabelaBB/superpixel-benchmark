
#include "gft_filtering3.h"

namespace gft{


  namespace Scene8{


    void ModeFilterLabel(Scene8 *label, float r){
      AdjRel3::AdjRel3 *A = AdjRel3::Spheric(r);
      AdjRel3::AdjVxl  *N;
      int fx,fy,fz,dp;
      int p,q,i,Lmax,l,lmax;
      int *frequency;
      Scene8 *sub=NULL,*mode=NULL,*tmp=NULL;
      struct{
	Voxel v1;
	Voxel v2;
      } box;
      
      MBB(label, &box.v1, &box.v2);
      sub = SubScene(label, box.v1, box.v2);

      AdjRel3::GetFrameSize(A, &fx, &fy, &fz);

      tmp = AddFrame(sub, fx, 0);
      Destroy(&sub);
      sub = tmp;
      mode = Clone(sub);
      
      Lmax = GetMaximumValue(sub);
      frequency = (int *)calloc(Lmax+1,sizeof(int));
      
      N = AdjRel3::AdjVoxels(A, sub);

      dp  = fx*1;
      dp += fy*sub->xsize;
      dp += fz*sub->xsize*sub->ysize;
      
      for(p=dp; p<sub->n-dp; p++){
	
	memset(frequency, 0, (Lmax+1)*sizeof(int));
	
	for(i=0; i<N->n; i++){
	  q = p + N->dp[i];
	  frequency[sub->data[q]]++;
	}
	lmax = sub->data[p];
	for(l=0; l<=Lmax; l++){
	  if(frequency[l]>frequency[lmax])
	    lmax = l;
	}
	mode->data[p] = lmax;
      }
      Destroy(&sub);

      sub = RemFrame(mode, fx);      
      Copy(label, sub, box.v1);
      
      Destroy(&sub);
      Destroy(&mode);
      AdjRel3::Destroy(&A);
      AdjRel3::DestroyAdjVxl(&N);
      free(frequency);
    }
    

  } //end Scene8 namespace



  
  namespace Scene32{


    void ModeFilterLabel(Scene32 *label, float r){
      AdjRel3::AdjRel3 *A = AdjRel3::Spheric(r);
      AdjRel3::AdjVxl  *N;
      int fx,fy,fz,dp;
      int p,q,i,Lmax,l,lmax;
      int *frequency;
      Scene32 *sub=NULL,*mode=NULL,*tmp=NULL;
      struct{
	Voxel v1;
	Voxel v2;
      } box;
      
      MBB(label, &box.v1, &box.v2);
      sub = SubScene(label, box.v1, box.v2);

      AdjRel3::GetFrameSize(A, &fx, &fy, &fz);

      tmp = AddFrame(sub, fx, 0);
      Destroy(&sub);
      sub = tmp;
      mode = Clone(sub);
      
      Lmax = GetMaximumValue(sub);
      frequency = (int *)calloc(Lmax+1,sizeof(int));
      
      N = AdjRel3::AdjVoxels(A, sub);

      dp  = fx*1;
      dp += fy*sub->xsize;
      dp += fz*sub->xsize*sub->ysize;
      
      for(p=dp; p<sub->n-dp; p++){
	
	memset(frequency, 0, (Lmax+1)*sizeof(int));
	
	for(i=0; i<N->n; i++){
	  q = p + N->dp[i];
	  frequency[sub->data[q]]++;
	}
	lmax = sub->data[p];
	for(l=0; l<=Lmax; l++){
	  if(frequency[l]>frequency[lmax])
	    lmax = l;
	}
	mode->data[p] = lmax;
      }
      Destroy(&sub);

      sub = RemFrame(mode, fx);      
      Copy(label, sub, box.v1);
      
      Destroy(&sub);
      Destroy(&mode);
      AdjRel3::Destroy(&A);
      AdjRel3::DestroyAdjVxl(&N);
      free(frequency);
    }
    

    
    Scene32  *AccAbsDiff(Scene32 *scn, float r){
      Scene32 *W;
      AdjRel3::AdjRel3 *A = AdjRel3::Spheric(r);
      int p,q,i,weight;
      Voxel u,v;
      
      W = Create(scn);
      for(p = 0; p < scn->n; p++){
	v.c.x = GetAddressX(scn, p);
	v.c.y = GetAddressY(scn, p);
	v.c.z = GetAddressZ(scn, p);
	
	for(i = 1; i < A->n; i++){
	  u.v = v.v + A->d[i].v;

	  if(gft::Scene32::IsValidVoxel(scn, u)){
	    q = gft::Scene32::GetVoxelAddress(scn, u);
	    weight = abs(scn->data[p] - scn->data[q]);
	    W->data[p] += weight;
	  }
	}
      }
      AdjRel3::Destroy(&A);
      return W;
    }


    
  } //end Scene32 namespace

  

} //end gft namespace





