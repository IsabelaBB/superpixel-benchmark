
#include "gft_arraylist.h"

namespace gft{
  namespace ArrayList{

    void myfree(void **mem);
    
    void myfree(void **mem){
      if(*mem!=NULL){
	free(*mem);
	*mem=NULL;
      }
    }


    ArrayList *Create(int cap){
      ArrayList *A;
      int i;
      
      A = (ArrayList *) calloc(1,sizeof(ArrayList));
      if(A == NULL)
	gft::Error((char *)MSG1,(char *)"ArrayList::Create");
      
      A->array = (void **) calloc(cap, sizeof(void *));
      if(A->array == NULL)
	gft::Error((char *)MSG1,(char *)"ArrayList::Create");
      
      for(i=0; i<cap; i++)
	A->array[i] = NULL;
      A->cap = cap;
      A->n = 0;
      A->clean = myfree;
      return A;
    }

    
    void Destroy(ArrayList **A){
      void (*clean)(void**);
      ArrayList *aux;
      int i;
      
      aux = *A;
      clean = aux->clean;
      if(aux != NULL){
	if(aux->array != NULL){
	  if(clean!=NULL){
	    for(i=0; i<aux->n; i++)
	      if(aux->array[i]!=NULL)
		(*clean)(&aux->array[i]);
	  }
	  free(aux->array);
	}
	free(aux);
	*A = NULL;
      }
    }


    void SetCleanFunc(ArrayList *A,
		      void (*clean)(void**)){
      A->clean = clean;
    }
    
    
    void AddElement(ArrayList *A, 
		    void *elem){
      int i;
      
      if(A->n < A->cap){
	A->array[A->n] = elem;
	A->n++;
      }
      else{
	A->cap = ROUND(A->cap*1.25)+1;
	A->array = (void **)realloc(A->array,
				    A->cap*sizeof(void *));
	if(A->array == NULL)
	  gft::Error((char *)MSG1,(char *)"ArrayList::AddElement");
	for(i=A->n; i<A->cap; i++)
	  A->array[i] = NULL;
	AddElement(A, elem);
      }
    }


    void *GetElement(ArrayList *A, 
		     int index){
      if(index<0 || index>=A->n)
	return NULL;
      return A->array[index];
    }
    

    void  DelElement(ArrayList *A, 
		     int index){
      void (*clean)(void**);
      int i;
      
      if(index<0 || index>=A->n)
	return;
      
      clean = A->clean;
      (*clean)(&A->array[index]);
      
      for(i=index+1; i<A->n; i++)
	A->array[i-1] = A->array[i];
      A->n--;
    }


    void  DelElement(ArrayList *A,
		     void **elem){
      int i;
      for(i=0; i<A->n; i++){
	if(*elem == GetElement(A, i)){
	  DelElement(A, i);
	  *elem=NULL;
	}
      }
    }
    

    void  Resize(ArrayList *A, int n){
      void (*clean)(void**);
      int n_old,i;
      
      n_old = A->n;
      
      A->cap = n;
      A->n = MIN(n,A->n);
      
      clean = A->clean;
      for(i=A->n; i<n_old; i++)
	if(A->array[i]!=NULL)
	  (*clean)(&A->array[i]);
      
      A->array = (void **)realloc(A->array,
				  n*sizeof(void *));
      if(A->array == NULL)
	gft::Error((char *)MSG1,(char *)"ArrayList::Resize");
      for(i=A->n; i<A->cap; i++)
	A->array[i] = NULL;
    }
    
    
    void  Trim2Size(ArrayList *A){
      Resize(A, A->n);
    }


  } /*end ArrayList namespace*/
} /*end gft namespace*/

