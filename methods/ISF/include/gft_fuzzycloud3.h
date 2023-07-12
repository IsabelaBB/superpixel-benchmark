#ifndef _GFT_FUZZYCLOUD3_H_
#define _GFT_FUZZYCLOUD3_H_

#include "gft_common.h"
#include "gft_adjrel3.h"
#include "gft_scene32.h"
#include "gft_filelist.h"
#include "gft_adjrel3f.h"
#include "gft_gradient3.h"

/*#include "cloud3.h"*/
/*#include "gradient.h"*/
/*#include "scene_addons.h"*/

#define MAX_PROB 100000

namespace gft{
  namespace RegionCloud3{

    typedef struct _regionCloud3 {
      AdjRel3f::AdjRel3f *disp;  //displacement.
      Scene32::Scene32  **prob;
      int nobjs;
      int nimages;
      //private:
      AdjRel3f::AdjRel3f *fdisp;  //float displacement.
    } RegionCloud3;

    RegionCloud3 *ByLabelList(FileList::FileList *L);
    void          Destroy(RegionCloud3 **rcloud);
    void          GetVoxelSize(RegionCloud3 *rcloud,
			       float *dx, 
			       float *dy, 
			       float *dz);
    RegionCloud3 *Subsampling(RegionCloud3 *rcloud);
    RegionCloud3 *LinearInterp(RegionCloud3 *rcloud,
			       float dx,float dy,float dz);
    RegionCloud3 *GaussianBlur(RegionCloud3 *rcloud);
    RegionCloud3 *ChangeOrientationToLPS(RegionCloud3 *rcloud,
					 char *ori);
    
    RegionCloud3 *Read(char *filename);
    void          Write(RegionCloud3 *rcloud, 
			char *filename);
    
    void RemoveElem(RegionCloud3 *rcloud,
		    Scene32::Scene32 *label);

  } //end RegionCloud3 namespace
} //end gft namespace
   
    //---------------------------------------------

namespace gft{
  namespace BorderCloud3{
    
    typedef struct _borderCloud3 {
      AdjRel3::AdjRel3 *disp;  //displacement.
      Gradient3::Gradient3 **prob;
      int nobjs;
      //private:
      AdjRel3f::AdjRel3f *fdisp;  //float displacement.
    } BorderCloud3;


    BorderCloud3 *ByRegionCloud(RegionCloud3::RegionCloud3 *rcloud);
    
    void          Destroy(BorderCloud3 **bcloud);
    void          GetVoxelSize(BorderCloud3 *bcloud,
			       float *dx, 
			       float *dy, 
			       float *dz);
    BorderCloud3 *Subsampling(BorderCloud3 *bcloud);
    BorderCloud3 *LinearInterp(BorderCloud3 *bcloud,
			       float dx,float dy,float dz);
    BorderCloud3 *ChangeOrientationToLPS(BorderCloud3 *bcloud,
					 char *ori);
    void          Normalize(BorderCloud3 *bcloud,
			    int omin,int omax,
			    int nmin,int nmax);
    
    BorderCloud3 *Read(char *filename);
    void          Write(BorderCloud3 *bcloud, 
			char *filename);
    
  } //end BorderCloud3 namespace
} //end gft namespace
   
#endif


