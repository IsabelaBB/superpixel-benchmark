
#include "gft_stack.h"

namespace gft{
  namespace Stack{

    Stack *Create(int n) {
      Stack *S;
      S = (Stack *) malloc(sizeof(Stack));
      if(S==NULL) gft::Error((char *)MSG1,
			     (char *)"Stack::Create");
      S->n   = n;
      S->top = -1;
      S->data = gft::AllocIntArray(n);
      return S;
    }

    void    Destroy(Stack **S) {
      Stack *aux = *S;
      if(aux) {
	if(aux->data) gft::FreeIntArray(&aux->data);
	free(aux);
	*S = NULL;
      }
    }

    void Clear(Stack *S){
      S->top = -1;
    }
    
    void    Push(Stack *S, int p) {
      S->data[++(S->top)] = p;
    }

    int     Pop(Stack *S) {
      if(S->top == -1) return -1;
      return(S->data[S->top--]);
    }

    
  } //end Stack namespace
} //end gft namespace

