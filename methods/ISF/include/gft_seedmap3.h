
#ifndef _GFT_SEEDMAP3_H_
#define _GFT_SEEDMAP3_H_

#include "gft_common.h"
#include "gft_scene.h"
#include "gft_adjregion3.h"
#include "gft_segmentation3.h"
#include "gft_fuzzycloud3.h"
#include "gft_gradient3.h"


namespace gft{
  /**
   * \brief Data structures for fast accessing seed pixels of a fuzzy model (e.g., CSM).
   */
  namespace AdjSeedmap3{


    typedef struct _adjseedmap3 {
      AdjRel3::AdjRel3 *disp;  //displacement.
      AdjRegion3::AdjRegion3 **uncertainty;
      AdjRegion3::AdjRegion3 **object;
      AdjRegion3::AdjRegion3 **obj_border;
      AdjRegion3::AdjRegion3 **bkg_border;
      int nobjs;
    } AdjSeedmap3;


    AdjSeedmap3 *Create(int nobjs);
    AdjSeedmap3 *Create(RegionCloud3::RegionCloud3 *rcloud);
    void         Destroy(AdjSeedmap3 **asmap);

    
    void DrawObject(Scene8::Scene8 *scn,
		    AdjSeedmap3 *asmap,
		    Voxel u, int l, uchar val);
    void DrawObjBorder(Scene8::Scene8 *scn,
		       AdjSeedmap3 *asmap,
		       Voxel u, int l, uchar val);
    void DrawBkgBorder(Scene8::Scene8 *scn,
		       AdjSeedmap3 *asmap,
		       Voxel u, int l, uchar val);
    void DrawUncertainty(Scene8::Scene8 *scn,
			 AdjSeedmap3 *asmap,
			 Voxel u, int l, uchar val);
    
    void CopyUncertainty(Scene8::Scene8 *dest, 
			 Scene8::Scene8 *src,
			 AdjSeedmap3 *asmap,
			 Voxel u, int l);

    void AddUncertainty(Scene8::Scene8 *dest, 
			Scene8::Scene8 *src,
			AdjSeedmap3 *asmap,
			Voxel u, int l);

    void CloudArcWeight(Scene16::Scene16 *arcw,
			Scene16::Scene16 *wobj,
			Gradient3::Gradient3 *grad,
			Voxel u,
			BorderCloud3::BorderCloud3 *bcloud,
			AdjSeedmap3 *asmap,
			int l, float w);
    

  } //end AdjSeedmap3 namespace
} //end gft namespace

#endif


