
#include "gft_marker.h"


namespace gft{
  namespace Image32{

    float MaxRadiusByErosion(Image32 *bin){
      Image32 *bin1=NULL,*edt=NULL;
      gft::AdjRel::AdjRel *A=NULL;
      int Imax,Emax;

      A = gft::AdjRel::Circular(1.0);
      bin1 = Threshold(bin, 0, 0);
      edt = Mask2EDT(bin1, A, INTERIOR, INT_MAX, 0);
      Emax = GetMaxVal(edt);
      Destroy(&bin1);
      Destroy(&edt);
      
      bin1 = Threshold(bin, 1, INT_MAX);
      edt = Mask2EDT(bin1, A, INTERIOR, INT_MAX, 0);
      Imax = GetMaxVal(edt);
      Destroy(&bin1);
      Destroy(&edt);
      gft::AdjRel::Destroy(&A);
      
      return MIN(sqrtf(Imax),sqrtf(Emax));
    }


    float MaxObjRadiusByErosion(Image32 *bin){
      Image32 *bin1=NULL,*edt=NULL;
      gft::AdjRel::AdjRel *A=NULL;
      int Imax;
      
      A = gft::AdjRel::Circular(1.0);
      bin1 = Threshold(bin, 1, INT_MAX);
      edt = Mask2EDT(bin1, A, INTERIOR, INT_MAX, 0);
      Imax = GetMaxVal(edt);
      Destroy(&bin1);
      Destroy(&edt);
      gft::AdjRel::Destroy(&A);
      
      return (sqrtf(Imax));
    }


    float MaxBkgRadiusByErosion(Image32 *bin){
      Image32 *bin1=NULL,*edt=NULL;
      gft::AdjRel::AdjRel *A=NULL;
      int Emax;
      
      A = gft::AdjRel::Circular(1.0);
      bin1 = Threshold(bin, 0, 0);
      edt = Mask2EDT(bin1, A, INTERIOR, INT_MAX, 0);
      Emax = GetMaxVal(edt);
      Destroy(&bin1);
      Destroy(&edt);
      gft::AdjRel::Destroy(&A);
      
      return (sqrtf(Emax));
    }



    Image32 *BkgMaskByErosion(Image32 *bin, float radius){
      gft::Set::Set *S=NULL;
      Image32 *bin1=NULL,*bin2=NULL;
      
      bin1 = Threshold(bin, 0, 0);
      bin2 = ErodeBin(bin1, &S, radius);
      Destroy(&bin1);
      gft::Set::Destroy(&S);
      return bin2;
    }
    

    Image32 *ObjMaskByErosion(Image32 *bin, float radius){
      gft::Set::Set *S=NULL;
      Image32 *bin1=NULL,*bin2=NULL;
      
      bin1 = Threshold(bin, 1, INT_MAX);
      bin2 = ErodeBin(bin1, &S, radius);
      Destroy(&bin1);
      gft::Set::Destroy(&S);
      return bin2;
    }



  } /*end Image32 namespace*/
} /*end gft namespace*/
