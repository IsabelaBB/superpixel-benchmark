
#ifndef _GFT_SCENE8_H_
#define _GFT_SCENE8_H_

#include "gft_common.h"

namespace gft{
  namespace Scene8{

  /**
   * It supports both linear and three-dimensional access 
   * (i.e., scn->data[p] or scn->array[z][y][x] for a voxel
   * (x,y,z) at address p=x+y*xsize+z*xsize*ysize).
   */
  typedef struct _scene8 {
    uchar *data;
    uchar ***array;
    int xsize,ysize,zsize;
    float dx,dy,dz;
    uchar maxval;
    int n;
  } Scene8;

  /**
   * \brief A constructor.
   */
  Scene8 *Create(int xsize,int ysize,int zsize);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  Scene8 *Create(Scene8 *scn);
  /**
   * \brief A destructor.
   */
  void    Destroy(Scene8 **scn);
  void    Copy(Scene8 *dest, Scene8 *src);
  void    Copy(Scene8 *dest, Scene8 *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  Scene8 *Clone(Scene8 *scn);
  Scene8 *SubScene(Scene8 *scn, Voxel l, Voxel h);
  Scene8 *SubScene(Scene8 *scn,
		   int xl, int yl, int zl,
		   int xh, int yh, int zh);
  void    Fill(Scene8 *scn, uchar value);

  void    Write(Scene8 *scn, char *filename);

  inline uchar GetValue(Scene8 *scn, Voxel v);
  inline uchar GetValue(Scene8 *scn, int p);
  inline uchar GetValue(Scene8 *scn, int x, int y, int z);
  uchar  GetValue_nn(Scene8 *scn, float x, float y, float z);

  inline int GetAddressX(Scene8 *scn, int p);
  inline int GetAddressY(Scene8 *scn, int p);
  inline int GetAddressZ(Scene8 *scn, int p);
  inline int GetVoxelAddress(Scene8 *scn, Voxel v);
  inline int GetVoxelAddress(Scene8 *scn, int x, int y, int z);

  inline bool IsValidVoxel(Scene8 *scn, int x, int y, int z);
  inline bool IsValidVoxel(Scene8 *scn, Voxel v);
  
  uchar   GetMaximumValue(Scene8 *scn);
  uchar   GetMinimumValue(Scene8 *scn);

  Scene8 *MBB(Scene8 *scn);
  void    MBB(Scene8 *scn, Voxel *l, Voxel *h);

  Voxel   Centroid(Scene8 *scn);

  Scene8 *AddFrame(Scene8 *scn,  int sz, uchar value);
  Scene8 *RemFrame(Scene8 *fscn, int sz);


  //---------inline definitions------------------
  inline bool  IsValidVoxel(Scene8 *scn, int x, int y, int z){
    if((x >= 0)&&(x < scn->xsize)&&
       (y >= 0)&&(y < scn->ysize)&&
       (z >= 0)&&(z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline bool IsValidVoxel(Scene8 *scn, Voxel v){
    if((v.c.x >= 0)&&(v.c.x < scn->xsize)&&
       (v.c.y >= 0)&&(v.c.y < scn->ysize)&&
       (v.c.z >= 0)&&(v.c.z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline uchar GetValue(Scene8 *scn, Voxel v){
    return (scn->array[v.c.z][v.c.y][v.c.x]);
  }
  inline uchar GetValue(Scene8 *scn, int p){
    return (scn->data[p]);
  }
  inline uchar GetValue(Scene8 *scn, int x, int y, int z){
    return (scn->array[z][y][x]);
  }

  inline int GetAddressX(Scene8 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))%scn->xsize);
  }
  inline int GetAddressY(Scene8 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))/scn->xsize);
  }
  inline int GetAddressZ(Scene8 *scn, int p){
    return (p/(scn->xsize*scn->ysize));
  }
  inline int GetVoxelAddress(Scene8 *scn, Voxel v){
    return (v.c.x + v.c.y*scn->xsize + 
	    v.c.z*scn->xsize*scn->ysize);
  }
  inline int GetVoxelAddress(Scene8 *scn, int x, int y, int z){
    return (x + y*scn->xsize + z*scn->xsize*scn->ysize);
  }

  
  } //end Scene8 namespace
} //end gft namespace


#include "gft_scene32.h"
#include "gft_scene16.h"

namespace gft{
  namespace Scene8{

    Scene32::Scene32 *ConvertTo32(Scene8 *scn);
    Scene16::Scene16 *ConvertTo16(Scene8 *scn);

  } //end Scene8 namespace
} //end gft namespace


#endif
