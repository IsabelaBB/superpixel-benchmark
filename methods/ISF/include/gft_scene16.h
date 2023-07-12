
#ifndef _GFT_SCENE16_H_
#define _GFT_SCENE16_H_

#include "gft_common.h"

namespace gft{
  namespace Scene16{

  /**
   * It supports both linear and three-dimensional access 
   * (i.e., scn->data[p] or scn->array[z][y][x] for a voxel
   * (x,y,z) at address p=x+y*xsize+z*xsize*ysize).
   */
  typedef struct _scene16 {
    ushort *data;
    ushort ***array;
    int xsize,ysize,zsize;
    float dx,dy,dz;
    ushort maxval;
    int n;
  } Scene16;


  /**
   * \brief A constructor.
   */
  Scene16 *Create(int xsize,int ysize,int zsize);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  Scene16 *Create(Scene16 *scn);
  /**
   * \brief A destructor.
   */
  void     Destroy(Scene16 **scn);
  void     Copy(Scene16 *dest, Scene16 *src);
  void     Copy(Scene16 *dest, Scene16 *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  Scene16 *Clone(Scene16 *scn);
  Scene16 *SubScene(Scene16 *scn, Voxel l, Voxel h);
  Scene16 *SubScene(Scene16 *scn,
		    int xl, int yl, int zl,
		    int xh, int yh, int zh);
  void     Fill(Scene16 *scn, ushort value);

  void     Write(Scene16 *scn, char *filename);

  inline ushort GetValue(Scene16 *scn, Voxel v);
  inline ushort GetValue(Scene16 *scn, int p);
  inline ushort GetValue(Scene16 *scn, int x, int y, int z);
  ushort GetValue_nn(Scene16 *scn, float x, float y, float z);

  inline int GetAddressX(Scene16 *scn, int p);
  inline int GetAddressY(Scene16 *scn, int p);
  inline int GetAddressZ(Scene16 *scn, int p);
  inline int GetVoxelAddress(Scene16 *scn, Voxel v);
  inline int GetVoxelAddress(Scene16 *scn, int x, int y, int z);

  inline bool IsValidVoxel(Scene16 *scn, int x, int y, int z);
  inline bool IsValidVoxel(Scene16 *scn, Voxel v);
  
  ushort GetMaximumValue(Scene16 *scn);
  ushort GetMinimumValue(Scene16 *scn);

  Scene16 *MBB(Scene16 *scn);
  void     MBB(Scene16 *scn, Voxel *l, Voxel *h);

  Scene16 *AddFrame(Scene16 *scn,  int sz, ushort value);
  Scene16 *RemFrame(Scene16 *fscn, int sz);


  //---------inline definitions------------------
  inline bool  IsValidVoxel(Scene16 *scn, int x, int y, int z){
    if((x >= 0)&&(x < scn->xsize)&&
       (y >= 0)&&(y < scn->ysize)&&
       (z >= 0)&&(z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline bool IsValidVoxel(Scene16 *scn, Voxel v){
    if((v.c.x >= 0)&&(v.c.x < scn->xsize)&&
       (v.c.y >= 0)&&(v.c.y < scn->ysize)&&
       (v.c.z >= 0)&&(v.c.z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline ushort  GetValue(Scene16 *scn, Voxel v){
    return (scn->array[v.c.z][v.c.y][v.c.x]);
  }
  inline ushort  GetValue(Scene16 *scn, int p){
    return (scn->data[p]);
  }
  inline ushort  GetValue(Scene16 *scn, int x, int y, int z){
    return (scn->array[z][y][x]);
  }

  inline int GetAddressX(Scene16 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))%scn->xsize);
  }
  inline int GetAddressY(Scene16 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))/scn->xsize);
  }
  inline int GetAddressZ(Scene16 *scn, int p){
    return (p/(scn->xsize*scn->ysize));
  }
  inline int GetVoxelAddress(Scene16 *scn, Voxel v){
    return (v.c.x + v.c.y*scn->xsize + 
	    v.c.z*scn->xsize*scn->ysize);
  }
  inline int GetVoxelAddress(Scene16 *scn, int x, int y, int z){
    return (x + y*scn->xsize + z*scn->xsize*scn->ysize);
  }

  
  } //end Scene16 namespace
} //end gft namespace


#include "gft_scene32.h"
#include "gft_scene8.h"

namespace gft{
  namespace Scene16{

    Scene32::Scene32 *ConvertTo32(Scene16 *scn);
    Scene8::Scene8   *ConvertTo8(Scene16 *scn);

  } //end Scene16 namespace
} //end gft namespace


#endif
