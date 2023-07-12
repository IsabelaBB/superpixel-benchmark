
#ifndef _GFT_ADJREL3_H_
#define _GFT_ADJREL3_H_

#include "gft_common.h"

namespace gft{
  namespace AdjRel3{

    typedef union _displacement3 {
      v4si v;
      int  data[4];
      struct{ int x,y,z; } axis;
    } Displacement3;

    typedef struct _adjrel3 {
      Displacement3 *d;
      int n;
    } AdjRel3;
    
    typedef struct _adjvxl {
      int *dp;
      int n;
    } AdjVxl;


    AdjRel3 *Create(int n);
    void     Destroy(AdjRel3 **A);
    AdjRel3 *Clone(AdjRel3 *A);

    AdjRel3 *Spheric(float r);
    AdjRel3 *Ellipsoid(float rx, float ry, float rz);
    AdjRel3 *SphericalShell(float inner_radius,
			    float outer_radius);
    AdjRel3 *Box(int xsize, int ysize, int zsize);


    void     Scale(AdjRel3 *A,
		   float Sx, float Sy, float Sz);
    void     ClipX(AdjRel3 *A, int lower, int higher);
    void     ClipY(AdjRel3 *A, int lower, int higher);
    void     ClipZ(AdjRel3 *A, int lower, int higher);

    int      GetFrameSize(AdjRel3 *A);
    void     GetFrameSize(AdjRel3 *A, 
			  int *sz_x, 
			  int *sz_y, 
			  int *sz_z);
    float   *GetDistanceArray(AdjRel3 *A);

    void     DestroyAdjVxl(AdjVxl **N);
    
  } //end AdjRel3 namespace
} //end gft namespace


#include "gft_scene.h"

namespace gft{
  namespace AdjRel3{

    AdjVxl  *AdjVoxels(AdjRel3 *A, int xsize, int ysize);
    AdjVxl  *AdjVoxels(AdjRel3 *A, Scene32::Scene32 *scn);
    AdjVxl  *AdjVoxels(AdjRel3 *A, Scene16::Scene16 *scn);
    AdjVxl  *AdjVoxels(AdjRel3 *A, Scene8::Scene8 *scn);
    AdjVxl  *AdjVoxels(AdjRel3 *A, Scene::Scene *scn);

    AdjRel3 *Spherical_mm(Scene32::Scene32 *scn, float r_mm);
    AdjRel3 *SphericalGrid_mm(float dx, float dy, float dz,
			      float r_mm, float spacement_mm);
    AdjRel3 *SphericalGrid_mm(Scene32::Scene32 *scn, 
			      float r_mm, float spacement_mm);

    float   *GetDistanceArray_mm(AdjRel3 *A, Scene32::Scene32 *scn);
   
  } //end AdjRel3 namespace
} //end gft namespace

#endif

