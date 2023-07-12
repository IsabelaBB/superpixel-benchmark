
#ifndef _GFT_SCENE32F_H_
#define _GFT_SCENE32F_H_

#include "gft_common.h"

namespace gft{
  namespace Scene32f{

  /**
   * It supports both linear and three-dimensional access 
   * (i.e., scn->data[p] or scn->array[z][y][x] for a voxel
   * (x,y,z) at address p=x+y*xsize+z*xsize*ysize).
   */
  typedef struct _scene32f {
    float *data;
    float ***array;
    int xsize,ysize,zsize;
    float dx,dy,dz;
    float maxval;
    int n;
  } Scene32f;

  /**
   * \brief A constructor.
   */
  Scene32f *Create(int xsize,int ysize,int zsize);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  Scene32f *Create(Scene32f *scn);
  /**
   * \brief A destructor.
   */
  void     Destroy(Scene32f **scn);
  void     Copy(Scene32f *dest, Scene32f *src);
  void     Copy(Scene32f *dest, Scene32f *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  Scene32f *Clone(Scene32f *scn);
  Scene32f *SubScene(Scene32f *scn, Voxel l, Voxel h);
  Scene32f *SubScene(Scene32f *scn,
		    int xl, int yl, int zl,
		    int xh, int yh, int zh);
  void     Fill(Scene32f *scn, float value);

  //Scene32f *Read(char *filename);
  void     Write(Scene32f *scn, char *filename);

  inline float   GetValue(Scene32f *scn, Voxel v);
  inline float   GetValue(Scene32f *scn, int p);
  inline float   GetValue(Scene32f *scn, int x, int y, int z);
  float GetValue_trilinear(Scene32f *scn, float x, float y, float z);
  float   GetValue_nn(Scene32f *scn, float x, float y, float z);

  inline int GetAddressX(Scene32f *scn, int p);
  inline int GetAddressY(Scene32f *scn, int p);
  inline int GetAddressZ(Scene32f *scn, int p);
  inline int GetVoxelAddress(Scene32f *scn, Voxel v);
  inline int GetVoxelAddress(Scene32f *scn, int x, int y, int z);

  inline bool  IsValidVoxel(Scene32f *scn, int x, int y, int z);
  inline bool  IsValidVoxel(Scene32f *scn, Voxel v);
  
  float      GetMaximumValue(Scene32f *scn);
  float      GetMinimumValue(Scene32f *scn);

  Scene32f *MBB(Scene32f *scn);
  void     MBB(Scene32f *scn, Voxel *l, Voxel *h);

  Scene32f *AddFrame(Scene32f *scn,  int sz, float value);
  Scene32f *RemFrame(Scene32f *fscn, int sz);


  //---------inline definitions------------------
  inline bool  IsValidVoxel(Scene32f *scn, int x, int y, int z){
    if((x >= 0)&&(x < scn->xsize)&&
       (y >= 0)&&(y < scn->ysize)&&
       (z >= 0)&&(z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline bool  IsValidVoxel(Scene32f *scn, Voxel v){
    if((v.c.x >= 0)&&(v.c.x < scn->xsize)&&
       (v.c.y >= 0)&&(v.c.y < scn->ysize)&&
       (v.c.z >= 0)&&(v.c.z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline float   GetValue(Scene32f *scn, Voxel v){
    return (scn->array[v.c.z][v.c.y][v.c.x]);
  }
  inline float   GetValue(Scene32f *scn, int p){
    return (scn->data[p]);
  }
  inline float   GetValue(Scene32f *scn, int x, int y, int z){
    return (scn->array[z][y][x]);
  }

  inline int GetAddressX(Scene32f *scn, int p){
    return ((p%(scn->xsize*scn->ysize))%scn->xsize);
  }
  inline int GetAddressY(Scene32f *scn, int p){
    return ((p%(scn->xsize*scn->ysize))/scn->xsize);
  }
  inline int GetAddressZ(Scene32f *scn, int p){
    return (p/(scn->xsize*scn->ysize));
  }
  inline int GetVoxelAddress(Scene32f *scn, Voxel v){
    return (v.c.x + v.c.y*scn->xsize + 
	    v.c.z*scn->xsize*scn->ysize);
  }
  inline int GetVoxelAddress(Scene32f *scn, int x, int y, int z){
    return (x + y*scn->xsize + z*scn->xsize*scn->ysize);
  }

  
  } //end Scene32f namespace
} //end gft namespace


#include "gft_scene16.h"
#include "gft_scene8.h"

namespace gft{
  namespace Scene32f{

    Scene16::Scene16 *ConvertTo16(Scene32f *scn);
    Scene8::Scene8   *ConvertTo8(Scene32f *scn);

  } //end Scene32f namespace
} //end gft namespace


#endif
