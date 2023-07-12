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
    } Heap;


    /* Auxiliary Functions */

    //#define HEAP_DAD(i) ((i - 1) / 2)
    //#define HEAP_LEFTSON(i) (2 * i + 1)
    //#define HEAP_RIGHTSON(i) (2 * i + 2)

#define HEAP_DAD(i) (i/2)
#define HEAP_LEFTSON(i) (2 * i)
#define HEAP_RIGHTSON(i) (2 * i + 1)
    
    char IsFull(Heap *H);
    char IsEmpty(Heap *H);
    Heap *Create_MaxPolicy(int n, float **cost);
    Heap *Create_MinPolicy(int n, float **cost);
    void Destroy(Heap **H);

    char Insert_MaxPolicy(Heap *H, int pixel);
    char Insert_MinPolicy(Heap *H, int pixel);
    
    char Remove_MaxPolicy(Heap *H, int *pixel);
    char Remove_MinPolicy(Heap *H, int *pixel);

    void Update_MaxPolicy(Heap *H, int p, float value);
    void Update_MinPolicy(Heap *H, int p, float value);

    void GoUp_MaxPolicy(Heap *H, int i);
    void GoUp_MinPolicy(Heap *H, int i);

    void GoDown_MaxPolicy(Heap *H, int i);
    void GoDown_MinPolicy(Heap *H, int i);

    void Reset(Heap *H);

    char Delete_MaxPolicy(Heap *H, int pixel);
    char Delete_MinPolicy(Heap *H, int pixel);
    
  } //end Heap namespace
} //end gft namespace

#endif

