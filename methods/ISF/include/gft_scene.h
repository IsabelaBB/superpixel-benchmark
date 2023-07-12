
#ifndef _GFT_SCENE_H_
#define _GFT_SCENE_H_

#include "gft_common.h"
#include "gft_scene8.h"
#include "gft_scene16.h"
#include "gft_scene32.h"
#include "gft_scene32f.h"
#include "gft_scene64f.h"

typedef enum {integer, reals} gft_SceneType;

namespace gft{
  namespace Scene{

  typedef struct _scene {
    uchar nbits;
    gft_SceneType type;
    union{
      Scene8::Scene8   *scn8;
      Scene16::Scene16 *scn16;
      Scene32::Scene32 *scn32;
      Scene32f::Scene32f *scn32f;
      Scene64f::Scene64f *scn64f;
    } ptr;
  } Scene;


  /**
   * \brief A constructor.
   */
  Scene *Create(int xsize,int ysize,int zsize, int nbits, gft_SceneType type);
  /**
   * \brief A constructor for integer data.
   */
  Scene *Create(int xsize,int ysize,int zsize, int nbits);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  Scene *Create(Scene *scn);
  /**
   * \brief A destructor.
   */
  void   Destroy(Scene **scn);
  void   Copy(Scene *dest, Scene *src);
  void   Copy(Scene *dest, Scene *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  Scene *Clone(Scene *scn);
  Scene *SubScene(Scene *scn, Voxel l, Voxel h);
  Scene *SubScene(Scene *scn,
		  int xl, int yl, int zl,
		  int xh, int yh, int zh);
  void   Fill(Scene *scn, int value);

  Scene *Read(char *filename);
  void   Write(Scene *scn, char *filename);

  void   SetValue(Scene *scn, int p, int value);

  int    GetValue(Scene *scn, Voxel v);
  int    GetValue(Scene *scn, int p);
  int    GetValue(Scene *scn, int x, int y, int z);
  int    GetValue_nn(Scene *scn, float x, float y, float z);

  double    GetValue64f(Scene *scn, Voxel v);
  double    GetValue64f(Scene *scn, int p);
  double    GetValue64f(Scene *scn, int x, int y, int z);
  double   GetValue64f_nn(Scene *scn, float x, float y, float z);

  int    GetNumberOfVoxels(Scene *scn);

  int    GetAddressX(Scene *scn, int p);
  int    GetAddressY(Scene *scn, int p);
  int    GetAddressZ(Scene *scn, int p);
  int    GetVoxelAddress(Scene *scn, Voxel v);
  int    GetVoxelAddress(Scene *scn, int x, int y, int z);

  bool   IsValidVoxel(Scene *scn, int x, int y, int z);
  bool   IsValidVoxel(Scene *scn, Voxel v);
  
  int    GetMaximumValue(Scene *scn);
  int    GetMinimumValue(Scene *scn);

  Scene *MBB(Scene *scn);
  void   MBB(Scene *scn, Voxel *l, Voxel *h);

  Scene *AddFrame(Scene *scn,  int sz, int value);
  Scene *RemFrame(Scene *fscn, int sz);

  } //end Scene namespace
} //end gft namespace

#endif

