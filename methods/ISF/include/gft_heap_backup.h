#ifndef _GFT_HEAP_H_
#define _GFT_HEAP_H_

#include "gft_common.h"
#include "gft_gpqueue_by_Falcao.h"


namespace gft{
  namespace Heap{

    typedef struct _heap {
      float *cost;
      char *color;
      int *pixel;
      int *pos;
      int last;
      int n;
      char removal_policy; /* 0 is MINVALUE and 1 is MAXVALUE */
    } Heap;


    /* Auxiliary Functions */

#define HEAP_DAD(i) ((i - 1) / 2)
#define HEAP_LEFTSON(i) (2 * i + 1)
#define HEAP_RIGHTSON(i) (2 * i + 2)

    void SetRemovalPolicy(Heap *H, char policy);
    char IsFull(Heap *H);
    char IsEmpty(Heap *H);
    Heap *Create(int n, float *cost);
    void Destroy(Heap **H);
    char Insert(Heap *H, int pixel);
    char Remove(Heap *H, int *pixel);
    void Update(Heap *H, int p, float value);
    void GoUp(Heap *H, int i);
    void GoDown(Heap *H, int i);
    void Reset(Heap *H);

    char Delete(Heap *H, int pixel);
    
  } //end Heap namespace
} //end gft namespace

#endif

