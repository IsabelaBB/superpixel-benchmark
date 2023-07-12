#ifndef _GFT_HEAP64F_H_
#define _GFT_HEAP64F_H_

#include "gft_common.h"
#include "gft_gpqueue_by_Falcao.h"

#include "gft_heap.h"

namespace gft{
  namespace Heap64f{

    typedef struct _heap64f {
      double *cost;
      char *color;
      int *pixel;
      int *pos;
      int last;
      int n;
    } Heap64f;


    /* Auxiliary Functions */

    char IsFull(Heap64f *H);
    char IsEmpty(Heap64f *H);
    Heap64f *Create(int n, double *cost);
    void Destroy(Heap64f **H);

    void Insert_MaxPolicy(Heap64f *H, int pixel);
    void Insert_MinPolicy(Heap64f *H, int pixel);
    
    void Remove_MaxPolicy(Heap64f *H, int *pixel);
    void Remove_MinPolicy(Heap64f *H, int *pixel);

    void Update_MaxPolicy(Heap64f *H, int p, double value);
    void Update_MinPolicy(Heap64f *H, int p, double value);

    void GoUp_MaxPolicy(Heap64f *H, int i);
    void GoUp_MinPolicy(Heap64f *H, int i);

    void GoDown_MaxPolicy(Heap64f *H, int i);
    void GoDown_MinPolicy(Heap64f *H, int i);

    void Reset(Heap64f *H);

    void Delete_MaxPolicy(Heap64f *H, int pixel);
    void Delete_MinPolicy(Heap64f *H, int pixel);
    

  } //end Heap64f namespace
} //end gft namespace

#endif

