
#ifndef _GFT_FILTERING3_H_
#define _GFT_FILTERING3_H_

#include "gft_common.h"
#include "gft_scene.h"
#include "gft_adjrel3.h"

namespace gft{


  namespace Scene8{

    void ModeFilterLabel(Scene8 *label, float r);

  } //end Scene8 namespace


  namespace Scene32{
    
    void ModeFilterLabel(Scene32 *label, float r);

    Scene32  *AccAbsDiff(Scene32 *scn, float r);
    
  } //end Scene32 namespace


} //end gft namespace


#endif










