
#ifndef _GFT_STACK_H
#define _GFT_STACK_H 1

#include "gft_common.h"

namespace gft{
  namespace Stack{

    typedef struct _stack {
      int *data;
      int top;
      int n;
    } Stack;
    
    Stack *Create(int n);
    void   Destroy(Stack **S);
    void   Push(Stack *S, int p);
    /**
     * @return Returns NIL if empty.
     */
    int    Pop(Stack *S);

    void Clear(Stack *S);
    
    inline int IsEmpty(Stack *S){ return(S->top == -1); }

  } //end Stack namespace
} //end gft namespace

#endif
