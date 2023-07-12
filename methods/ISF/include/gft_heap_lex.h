#ifndef _GFT_HEAP_LEX_H_
#define _GFT_HEAP_LEX_H_

#include "gft_common.h"
#include "gft_heap.h"
#include "gft_gpqueue_by_Falcao.h"

namespace gft{
  namespace Heap_lex{

    typedef struct _heap_lex {
      float *cost1,*cost2;
      char *color;
      int *pixel;
      int *pos;
      int last;
      int n;
      char removal_policy; /* 0 is MINVALUE and 1 is MAXVALUE */
    } Heap_lex;
   

    void SetRemovalPolicy(Heap_lex *H, char policy);
    char IsFull(Heap_lex *H);
    char IsEmpty(Heap_lex *H);
    Heap_lex *Create(int n, float *cost1, float *cost2);
    void Destroy(Heap_lex **H);
    char Insert(Heap_lex *H, int pixel);
    char Remove(Heap_lex *H, int *pixel);
    void Update(Heap_lex *H, int p, float value1, float value2);
    void GoUp(Heap_lex *H, int i);
    void GoDown(Heap_lex *H, int i);
    void Reset(Heap_lex *H);

  } //end Heap_lex namespace
} //end gft namespace

#endif



