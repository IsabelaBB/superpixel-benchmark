
#ifndef _GFT_ADJREGION3_H_
#define _GFT_ADJREGION3_H_

#include "gft_common.h"
#include "gft_scene.h"
#include "gft_adjrel3.h"

namespace gft{
  namespace AdjRegion3{

    typedef union _displacement3 {
      v4si v;
      int  data[4];
      struct{ int x,y,z; } axis;
    } Displacement3;

    typedef struct _adjregion3 {
      Displacement3 *d;
      int n;
      int *dp;

      //private:
      int xsize;
      int xysize;
      int max[3];
      int min[3];
    } AdjRegion3;


    AdjRegion3 *Create(int n);
    AdjRegion3 *Create(AdjRel3::AdjRel3 *A);
    AdjRegion3 *Create(Scene8::Scene8 *mask, Voxel Ref);
    void        Destroy(AdjRegion3 **adjreg);
    AdjRegion3 *Clone(AdjRegion3 *adjreg);
    AdjRegion3 *Merge(AdjRegion3 *r1, AdjRegion3 *r2);

    Scene8::Scene8 *Export2Mask(AdjRegion3 *adjreg);

    void    Draw(AdjRegion3 *adjreg,
		 Scene8::Scene8 *scn,
		 Voxel u,
		 uchar val);
    void    Draw(AdjRegion3 *adjreg,
		 Scene16::Scene16 *scn,
		 Voxel u,
		 ushort val);
    void    Draw(AdjRegion3 *adjreg,
		 Scene32::Scene32 *scn,
		 Voxel u,
		 int val);
    /*
    void    MT_Draw(AdjRegion3 *adjreg,
		    Scene8::Scene8 *scn,
		    Voxel u,
		    uchar val);
    */
    void    DrawOpt(AdjRegion3 *adjreg,
		    Scene8::Scene8 *scn,
		    int p, uchar val);
    void    DrawOpt(AdjRegion3 *adjreg,
		    Scene16::Scene16 *scn,
		    int p, ushort val);
    void    DrawOpt(AdjRegion3 *adjreg,
		    Scene32::Scene32 *scn,
		    int p, int val);
    /*
    void    MT_DrawOpt(AdjRegion3 *adjreg,
		       Scene8::Scene8 *scn,
		       int p, uchar val);
    */
    void    Optimize(AdjRegion3 *adjreg,
		     int xsize, int ysize);
    void    Optimize(AdjRegion3 *adjreg,
		     Scene8::Scene8 *scn);
    void    Optimize(AdjRegion3 *adjreg,
		     Scene16::Scene16 *scn);
    void    Optimize(AdjRegion3 *adjreg,
		     Scene32::Scene32 *scn);

    void    RefreshLimits(AdjRegion3 *adjreg);
    void    GetLimits(AdjRegion3 *adjreg,
		      int *dx_min, int *dy_min, int *dz_min,
		      int *dx_max, int *dy_max, int *dz_max);

    bool    FitInside(AdjRegion3 *adjreg, Voxel vx,
		      int xsize, int ysize, int zsize, int sz);
    bool    FitInside(AdjRegion3 *adjreg, Voxel vx,
		      Scene8::Scene8 *scn, int sz);
    bool    FitInside(AdjRegion3 *adjreg, Voxel vx,
		      Scene16::Scene16 *scn, int sz);
    bool    FitInside(AdjRegion3 *adjreg, Voxel vx,
		      Scene32::Scene32 *scn, int sz);

    float   InnerMean(AdjRegion3 *adjreg,
		      Voxel vx,
		      Scene16::Scene16 *scn);

    float   InnerSum(AdjRegion3 *adjreg,
		     Voxel vx,
		     Scene16::Scene16 *scn);


  } //end AdjRegion3 namespace
} //end gft namespace
    
#endif

