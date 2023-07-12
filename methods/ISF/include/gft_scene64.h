
#ifndef _GFT_SCENE64_H_
#define _GFT_SCENE64_H_

#include "gft_common.h"

extern "C" {
#include "nifti1_io.h"
}

namespace gft{
  namespace Scene64{

  /**
   * It supports both linear and three-dimensional access 
   * (i.e., scn->data[p] or scn->array[z][y][x] for a voxel
   * (x,y,z) at address p=x+y*xsize+z*xsize*ysize).
   */
  typedef struct _scene64 {
    long long *data;
    long long ***array;
    int xsize,ysize,zsize;
    float dx,dy,dz;
    long long maxval;
    int n;
    nifti_image *nii_hdr;
  } Scene64;

  /**
   * \brief A constructor.
   */
  Scene64 *Create(int xsize,int ysize,int zsize);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  Scene64 *Create(Scene64 *scn);
  /**
   * \brief A destructor.
   */
  void     Destroy(Scene64 **scn);
  void     Copy(Scene64 *dest, Scene64 *src);
  void     Copy(Scene64 *dest, Scene64 *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  Scene64 *Clone(Scene64 *scn);
  Scene64 *SubScene(Scene64 *scn, Voxel l, Voxel h);
  Scene64 *SubScene(Scene64 *scn,
		    int xl, int yl, int zl,
		    int xh, int yh, int zh);
  void     Fill(Scene64 *scn, long long value);


  inline int GetAddressX(Scene64 *scn, int p);
  inline int GetAddressY(Scene64 *scn, int p);
  inline int GetAddressZ(Scene64 *scn, int p);
  inline int GetVoxelAddress(Scene64 *scn, Voxel v);
  inline int GetVoxelAddress(Scene64 *scn, int x, int y, int z);

  inline bool  IsValidVoxel(Scene64 *scn, int x, int y, int z);
  inline bool  IsValidVoxel(Scene64 *scn, Voxel v);
  
  long long    GetMaximumValue(Scene64 *scn);
  long long    GetMinimumValue(Scene64 *scn);


  //---------inline definitions------------------
  inline bool  IsValidVoxel(Scene64 *scn, int x, int y, int z){
    if((x >= 0)&&(x < scn->xsize)&&
       (y >= 0)&&(y < scn->ysize)&&
       (z >= 0)&&(z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline bool  IsValidVoxel(Scene64 *scn, Voxel v){
    if((v.c.x >= 0)&&(v.c.x < scn->xsize)&&
       (v.c.y >= 0)&&(v.c.y < scn->ysize)&&
       (v.c.z >= 0)&&(v.c.z < scn->zsize))
      return(true);
    else
      return(false);
  }

  
  inline int GetAddressX(Scene64 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))%scn->xsize);
  }
  inline int GetAddressY(Scene64 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))/scn->xsize);
  }
  inline int GetAddressZ(Scene64 *scn, int p){
    return (p/(scn->xsize*scn->ysize));
  }
  inline int GetVoxelAddress(Scene64 *scn, Voxel v){
    return (v.c.x + v.c.y*scn->xsize + 
	    v.c.z*scn->xsize*scn->ysize);
  }
  inline int GetVoxelAddress(Scene64 *scn, int x, int y, int z){
    return (x + y*scn->xsize + z*scn->xsize*scn->ysize);
  }


  //The dimensions of the window (i.e., xsize,ysize,zsize) 
  //should be given in millimeters.
  float ComputeWindowSum(Scene64 *Iscn,
			 float xsize, 
			 float ysize, 
			 float zsize,
			 Voxel u);

  //The dimensions of the window (i.e., xsize,ysize,zsize) 
  //should be given in millimeters.
  float ComputeWindowDensity(Scene64 *Iscn,
			     float xsize, 
			     float ysize, 
			     float zsize,
			     Voxel u);
  
  //The dimensions of the window (i.e., xsize,ysize,zsize) 
  //should be given in millimeters.
  Voxel  FindWindowOfMaximumSum(Scene64 *Iscn, 
				float xsize, 
				float ysize, 
				float zsize);
  
  
  } //end Scene64 namespace
} //end gft namespace


#endif
