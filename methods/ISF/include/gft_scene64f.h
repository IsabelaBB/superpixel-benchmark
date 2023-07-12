
#ifndef _GFT_SCENE64F_H_
#define _GFT_SCENE64F_H_

#include "gft_common.h"

namespace gft{
  namespace Scene64f{

  /**
   * It supports both linear and three-dimensional access 
   * (i.e., scn->data[p] or scn->array[z][y][x] for a voxel
   * (x,y,z) at address p=x+y*xsize+z*xsize*ysize).
   */
  typedef struct _scene64f {
    double *data;
    double ***array;
    int xsize,ysize,zsize;
    float dx,dy,dz;
    double maxval;
    int n;
  } Scene64f;

  /**
   * \brief A constructor.
   */
  Scene64f *Create(int xsize,int ysize,int zsize);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  Scene64f *Create(Scene64f *scn);
  /**
   * \brief A destructor.
   */
  void     Destroy(Scene64f **scn);
  void     Copy(Scene64f *dest, Scene64f *src);
  void     Copy(Scene64f *dest, Scene64f *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  Scene64f *Clone(Scene64f *scn);
  Scene64f *SubScene(Scene64f *scn, Voxel l, Voxel h);
  Scene64f *SubScene(Scene64f *scn,
		    int xl, int yl, int zl,
		    int xh, int yh, int zh);
  void     Fill(Scene64f *scn, double value);

  //Scene64f *Read(char *filename);
  void     Write(Scene64f *scn, char *filename);

  inline double   GetValue(Scene64f *scn, Voxel v);
  inline double   GetValue(Scene64f *scn, int p);
  inline double   GetValue(Scene64f *scn, int x, int y, int z);
  double GetValue_trilinear(Scene64f *scn, float x, float y, float z);
  double   GetValue_nn(Scene64f *scn, float x, float y, float z);

  inline int GetAddressX(Scene64f *scn, int p);
  inline int GetAddressY(Scene64f *scn, int p);
  inline int GetAddressZ(Scene64f *scn, int p);
  inline int GetVoxelAddress(Scene64f *scn, Voxel v);
  inline int GetVoxelAddress(Scene64f *scn, int x, int y, int z);

  inline bool  IsValidVoxel(Scene64f *scn, int x, int y, int z);
  inline bool  IsValidVoxel(Scene64f *scn, Voxel v);
  
  double      GetMaximumValue(Scene64f *scn);
  double      GetMinimumValue(Scene64f *scn);

  Scene64f *MBB(Scene64f *scn);
  void     MBB(Scene64f *scn, Voxel *l, Voxel *h);

  Scene64f *AddFrame(Scene64f *scn,  int sz, double value);
  Scene64f *RemFrame(Scene64f *fscn, int sz);


  //---------inline definitions------------------
  inline bool  IsValidVoxel(Scene64f *scn, int x, int y, int z){
    if((x >= 0)&&(x < scn->xsize)&&
       (y >= 0)&&(y < scn->ysize)&&
       (z >= 0)&&(z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline bool  IsValidVoxel(Scene64f *scn, Voxel v){
    if((v.c.x >= 0)&&(v.c.x < scn->xsize)&&
       (v.c.y >= 0)&&(v.c.y < scn->ysize)&&
       (v.c.z >= 0)&&(v.c.z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline double   GetValue(Scene64f *scn, Voxel v){
    return (scn->array[v.c.z][v.c.y][v.c.x]);
  }
  inline double   GetValue(Scene64f *scn, int p){
    return (scn->data[p]);
  }
  inline double   GetValue(Scene64f *scn, int x, int y, int z){
    return (scn->array[z][y][x]);
  }

  inline int GetAddressX(Scene64f *scn, int p){
    return ((p%(scn->xsize*scn->ysize))%scn->xsize);
  }
  inline int GetAddressY(Scene64f *scn, int p){
    return ((p%(scn->xsize*scn->ysize))/scn->xsize);
  }
  inline int GetAddressZ(Scene64f *scn, int p){
    return (p/(scn->xsize*scn->ysize));
  }
  inline int GetVoxelAddress(Scene64f *scn, Voxel v){
    return (v.c.x + v.c.y*scn->xsize + 
	    v.c.z*scn->xsize*scn->ysize);
  }
  inline int GetVoxelAddress(Scene64f *scn, int x, int y, int z){
    return (x + y*scn->xsize + z*scn->xsize*scn->ysize);
  }

  
  } //end Scene64f namespace
} //end gft namespace


#include "gft_scene16.h"
#include "gft_scene8.h"

namespace gft{
  namespace Scene64f{

    Scene16::Scene16 *ConvertTo16(Scene64f *scn);
    Scene8::Scene8   *ConvertTo8(Scene64f *scn);

  } //end Scene64f namespace
} //end gft namespace


#endif
