
#include "gft_evaluation.h"

namespace gft{
  namespace Image32{
    
    float DiceSimilarity(Image32 *mask1, 
			 Image32 *mask2){
      float nelems_intersec,nelems_union;
      int  p,n;
      
      /* compute similarity between shapes */
      n = mask1->ncols * mask1->nrows;
      nelems_intersec = nelems_union = 0.0;
      for(p=0; p<n; p++){
	if(mask1->data[p]>0){
	  nelems_union++;
	  if(mask2->data[p]>0)
	    nelems_intersec++;
	}
	else{
	  if(mask2->data[p]>0)
	    nelems_union++;
	}
      }
      return((2.0*nelems_intersec)/(nelems_union+nelems_intersec));
    }


    float    JaccardSimilarity(Image32 *mask1,
			       Image32 *mask2){
      float nelems_intersec,nelems_union;
      int  p,n;
      
      /* compute similarity between shapes */
      n = mask1->ncols * mask1->nrows;
      nelems_intersec = nelems_union = 0.0;
      for(p=0; p<n; p++){
	if(mask1->data[p]>0){
	  nelems_union++;
	  if(mask2->data[p]>0)
	    nelems_intersec++;
	}
	else{
	  if(mask2->data[p]>0)
	    nelems_union++;
	}
      }
      return(nelems_intersec/nelems_union);
    }
    
    
    
    //mask1: Ground Truth
    //mask2: Segmentation Result
    int     AssessTP(Image32 *mask1, Image32 *mask2){
      int p,n,tp=0;
      
      n = mask1->ncols * mask1->nrows;
      for(p=0; p<n; p++){
	if(mask1->data[p]>0 && mask2->data[p]>0)
	  tp++;
      }
      return tp;
    }
    
    
    //mask1: Ground Truth
    //mask2: Segmentation Result
    int     AssessFN(Image32 *mask1, Image32 *mask2){
      int p,n,fn=0;
      
      n = mask1->ncols * mask1->nrows;
      for(p=0; p<n; p++){
	if(mask1->data[p]>0 && mask2->data[p]==0)
	  fn++;
      }
      return fn;
    }
    
    //mask1: Ground Truth
    //mask2: Segmentation Result
    int     AssessFP(Image32 *mask1, Image32 *mask2){
      int p,n,fp=0;
      
      n = mask1->ncols * mask1->nrows;
      for(p=0; p<n; p++){
	if(mask1->data[p]==0 && mask2->data[p]>0)
	  fp++;
      }
      return fp;
    }


    //mask1: Ground Truth
    //mask2: Segmentation Result
    int     AssessTN(Image32 *mask1, Image32 *mask2){
      int p,n,tn=0;
      
      n = mask1->ncols * mask1->nrows;
      for(p=0; p<n; p++){
	if(mask1->data[p]==0 && mask2->data[p]==0)
	  tn++;
      }
      return tn;
    }



    Image32 *GetObjError(Image32 *gtruth,
			 Image32 *mask){
      Image32 *err = NULL;
      int p,n;
      
      n = mask->ncols*mask->nrows;
      err = Create(mask->ncols, mask->nrows);
      
      for(p=0; p<n; p++){
	if(gtruth->data[p] > 0){
	  if(mask->data[p] == 0)
	    err->data[p] = 1;
	}
      }
      return err;
    }


    Image32 *GetBkgError(Image32 *gtruth,
			 Image32 *mask){
      Image32 *err = NULL;
      int p,n;
      
      n = mask->ncols*mask->nrows;
      err = Create(mask->ncols, mask->nrows);
      
      for(p=0; p<n; p++){
	if(gtruth->data[p] == 0){
	  if(mask->data[p] > 0)
	    err->data[p] = 1;
	}
      }
      return err;
    }
    

    //mask1: Ground Truth
    //mask2: Segmentation Result
    float    BoundaryError(Image32 *mask1,
			   Image32 *mask2){
      Image32 *border,*dist2 = NULL;
      gft::AdjRel::AdjRel *A = gft::AdjRel::Circular(1.5);
      float mean,d = 0.0;
      int n,p,nb = 0;
      
      //dist2  = DistTrans(mask1, A, BOTH);
      border = GetObjBorder(mask2);

      n = mask1->ncols*mask1->nrows;
      for(p=0; p<n; p++){
	if(border->data[p]>0){
	  nb++;
	  d += sqrtf((float)dist2->data[p]);
	}
      }
      Destroy(&border);
      Destroy(&dist2);
      gft::AdjRel::Destroy(&A);
      
      mean = (d/(float)nb);
      return mean;
    }


    //mask1: Ground Truth
    //mask2: Segmentation Result
    float    BoundaryFP(Image32 *mask1,
			Image32 *mask2){
      Image32 *border,*dist2 = NULL;
      gft::AdjRel::AdjRel *A = gft::AdjRel::Circular(1.5);
      float mean,d = 0.0;
      int n,p,nb = 0;

      //dist2  = SignedDistTrans(mask1, A, BOTH);
      border = GetObjBorder(mask2);

      n = mask1->ncols*mask1->nrows;
      for(p=0; p<n; p++){
	if(border->data[p]>0){
	  nb++;
	  if(dist2->data[p]<0)
	    d += sqrtf((float)(-dist2->data[p]));
	}
      }
      
      Destroy(&border);
      Destroy(&dist2);
      gft::AdjRel::Destroy(&A);
      
      mean = (d/(float)nb);
      return mean;
    }


    //mask1: Ground Truth
    //mask2: Segmentation Result
    float    BoundaryFN(Image32 *mask1,
			Image32 *mask2){
      Image32 *border,*dist2 = NULL;
      gft::AdjRel::AdjRel *A = gft::AdjRel::Circular(1.5);
      float mean,d = 0.0;
      int n,p,nb = 0;
      
      //dist2  = SignedDistTrans(mask1, A, BOTH);
      border = GetObjBorder(mask2);
      
      n = mask1->ncols*mask1->nrows;
      for(p=0; p<n; p++){
	if(border->data[p]>0){
	  nb++;
	  if(dist2->data[p]>0)
	    d += sqrtf((float)dist2->data[p]);
	}
      }
      
      Destroy(&border);
      Destroy(&dist2);
      gft::AdjRel::Destroy(&A);
      
      mean = (d/(float)nb);
      return mean;
    }

  } //end Image32 namespace


  namespace Scene32{
  
    float DiceSimilarity(Scene32 *mask1, 
			 Scene32 *mask2){
      float nelems_intersec,nelems_union;
      int  p,n;
      
      /* compute similarity between shapes */
      n = mask1->n;
      nelems_intersec = nelems_union = 0.0;
      for(p = 0; p < n; p++){
	if(mask1->data[p]>0){
	  nelems_union++;
	  if(mask2->data[p]>0)
	    nelems_intersec++;
	}
	else{
	  if(mask2->data[p]>0)
	    nelems_union++;
	}
      }
      return((2.0*nelems_intersec)/(nelems_union+nelems_intersec));
    }
    
    
  } //end Scene32 namespace

} //end gft namespace


