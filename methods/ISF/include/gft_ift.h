
#ifndef _GFT_IFT_H_
#define _GFT_IFT_H_

#include "gft_common.h"
#include "gft_sparsegraph.h"
#include "gft_pqueue32.h"
#include "gft_queue.h"
#include "gft_heap.h"

//----------------------------------------
//  Method wrappers:

//The 'label' scene should be pre-initialized as follows:
//  label->val[p]=NIL, unlabeled voxel.
//  label->val[p]=0,   background voxel.
//  label->val[p]=1,   object#1 voxel.
//  label->val[p]=2,   object#2 voxel.
//  ...
//---------------------------------------

namespace gft{
  namespace ift{

    void method_IFTW_FIFO_InnerCut(SparseGraph::SparseGraph *sg,
				   int *S,
				   Image32::Image32 *label);
    void method_IFTW_FIFO_OuterCut(SparseGraph::SparseGraph *sg,
				   int *S,
				   Image32::Image32 *label);

    Image32::Image32 *pred_IFTSUM(SparseGraph::SparseGraph *sg,
				  int *S,
				  Image32::Image32 *label,
				  float power, int type);
    
    //---------------------------------------
    /*    
    Image32::Image32 *EDT(Image32::Image32 *bin, float r);
    */
  } //end ift namespace
} //end gft namespace

#endif

