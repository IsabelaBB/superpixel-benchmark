
#ifndef _GFT_SCENE32_H_
#define _GFT_SCENE32_H_

#include "gft_common.h"
#include "gft_image32.h"

extern "C" {
#include "nifti1_io.h"
}

namespace gft{
  namespace Scene32{

  /**
   * It supports both linear and three-dimensional access 
   * (i.e., scn->data[p] or scn->array[z][y][x] for a voxel
   * (x,y,z) at address p=x+y*xsize+z*xsize*ysize).
   */
  typedef struct _scene32 {
    int *data;
    int ***array;
    int xsize,ysize,zsize;
    float dx,dy,dz;
    int maxval, n;
    nifti_image *nii_hdr;
  } Scene32;

  /**
   * \brief A constructor.
   */
  Scene32 *Create(int xsize,int ysize,int zsize);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  Scene32 *Create(Scene32 *scn);
  /**
   * \brief A destructor.
   */
  void     Destroy(Scene32 **scn);
  void     Copy(Scene32 *dest, Scene32 *src);
  void     Copy(Scene32 *dest, Scene32 *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  Scene32 *Clone(Scene32 *scn);
  Scene32 *SubScene(Scene32 *scn, Voxel l, Voxel h);
  Scene32 *SubScene(Scene32 *scn,
		    int xl, int yl, int zl,
		    int xh, int yh, int zh);
  void     Fill(Scene32 *scn, int value);

  Scene32 *Read(char *filename);
  Scene32 *ReadCompressed(char *filename);
  Scene32 *ReadNifti1(char *filename);
  void     Write(Scene32 *scn, char *filename);
  void     WriteCompressed(Scene32 *scn, char *filename);
  void     WriteNifti1(Scene32 *scn, char *filename);

  inline int   GetValue(Scene32 *scn, Voxel v);
  inline int   GetValue(Scene32 *scn, int p);
  inline int   GetValue(Scene32 *scn, int x, int y, int z);
  float GetValue_trilinear(Scene32 *scn, float x, float y, float z);
  int   GetValue_nn(Scene32 *scn, float x, float y, float z);

  inline int GetAddressX(Scene32 *scn, int p);
  inline int GetAddressY(Scene32 *scn, int p);
  inline int GetAddressZ(Scene32 *scn, int p);
  inline int GetVoxelAddress(Scene32 *scn, Voxel v);
  inline int GetVoxelAddress(Scene32 *scn, int x, int y, int z);

  inline bool  IsValidVoxel(Scene32 *scn, int x, int y, int z);
  inline bool  IsValidVoxel(Scene32 *scn, Voxel v);
  
  int      GetMaximumValue(Scene32 *scn);
  int      GetMinimumValue(Scene32 *scn);

  gft::Image32::Image32 *GetSliceX(Scene32 *scn, int x);
  gft::Image32::Image32 *GetSliceY(Scene32 *scn, int y);
  gft::Image32::Image32 *GetSliceZ(Scene32 *scn, int z);

  void PutSliceX(Scene32 *scn, Image32::Image32 *img, int x);
  void PutSliceY(Scene32 *scn, Image32::Image32 *img, int y);
  void PutSliceZ(Scene32 *scn, Image32::Image32 *img, int z);
  
  Scene32 *MBB(Scene32 *scn);
  void     MBB(Scene32 *scn, Voxel *l, Voxel *h);

  Scene32 *AddFrame(Scene32 *scn,  int sz, int value);
  Scene32 *RemFrame(Scene32 *fscn, int sz);


  //---------inline definitions------------------
  inline bool  IsValidVoxel(Scene32 *scn, int x, int y, int z){
    if((x >= 0)&&(x < scn->xsize)&&
       (y >= 0)&&(y < scn->ysize)&&
       (z >= 0)&&(z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline bool  IsValidVoxel(Scene32 *scn, Voxel v){
    if((v.c.x >= 0)&&(v.c.x < scn->xsize)&&
       (v.c.y >= 0)&&(v.c.y < scn->ysize)&&
       (v.c.z >= 0)&&(v.c.z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline int   GetValue(Scene32 *scn, Voxel v){
    return (scn->array[v.c.z][v.c.y][v.c.x]);
  }
  inline int   GetValue(Scene32 *scn, int p){
    return (scn->data[p]);
  }
  inline int   GetValue(Scene32 *scn, int x, int y, int z){
    return (scn->array[z][y][x]);
  }

  inline int GetAddressX(Scene32 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))%scn->xsize);
  }
  inline int GetAddressY(Scene32 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))/scn->xsize);
  }
  inline int GetAddressZ(Scene32 *scn, int p){
    return (p/(scn->xsize*scn->ysize));
  }
  inline int GetVoxelAddress(Scene32 *scn, Voxel v){
    return (v.c.x + v.c.y*scn->xsize + 
	    v.c.z*scn->xsize*scn->ysize);
  }
  inline int GetVoxelAddress(Scene32 *scn, int x, int y, int z){
    return (x + y*scn->xsize + z*scn->xsize*scn->ysize);
  }

  
  } //end Scene32 namespace
} //end gft namespace


#include "gft_scene16.h"
#include "gft_scene8.h"

namespace gft{
  namespace Scene32{

    Scene16::Scene16 *ConvertTo16(Scene32 *scn);
    Scene8::Scene8   *ConvertTo8(Scene32 *scn);

  } //end Scene32 namespace
} //end gft namespace


#include "gft_scene64.h"

namespace gft{
  namespace Scene32{

    Scene64::Scene64 *ComputeIntegralScene(Scene32 *scn);

    //The dimensions of the window (i.e., xsize,ysize,zsize) 
    //should be given in millimeters.
    void DrawWindow(Scene32 *scn,
		    int val,
		    float xsize, 
		    float ysize, 
		    float zsize,
		    Voxel u);

    //The dimensions of the window (i.e., xsize,ysize,zsize) 
    //should be given in millimeters.
    int WindowNvoxels(Scene32 *scn,
		      float xsize,
		      float ysize,
		      float zsize);
    
  } //end Scene32 namespace
} //end gft namespace


#endif
