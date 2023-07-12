
#include "gft_queue.h"

namespace gft{
  namespace Queue{


    Queue *Create(int nbuckets){
      Queue *Q;
      
      Q = (Queue *) malloc(sizeof(Queue));
      if(Q==NULL) gft::Error((char *)MSG1,
			     (char *)"Queue::Create");
      Q->nbuckets = nbuckets;
      Q->nadded = 0;
      Q->put = Q->get = 0;
      Q->data = gft::AllocIntArray(nbuckets);
      return Q;
    }


    void    Destroy(Queue **Q){
      Queue *aux;
      aux = *Q;
      if(aux!=NULL){
	if(aux->data!=NULL) 
	  free(aux->data);
	free(aux);
	*Q = NULL;
      }  
    }
    
    
    void    Push(Queue *Q, int p){
      if(!IsFull(Q)){
	Q->nadded++;
	Q->data[Q->put] = p;
	Q->put = (Q->put + 1) % Q->nbuckets;
      }
    }
    

    /* returns NIL if empty. */
    int     Pop(Queue *Q){ 
      int v;
      if(IsEmpty(Q)) return NIL;
      v = Q->data[Q->get];
      Q->get = (Q->get + 1) % Q->nbuckets;
      Q->nadded--;
      return v;
    }
    

    void    Reset(Queue *Q){
      Q->put = Q->get = 0;
      Q->nadded = 0;
    }
    

    bool IsEmpty(Queue *Q){ 
      return (Q->nadded==0); 
    }
    
    
    bool IsFull(Queue *Q){ 
      return (Q->nadded==Q->nbuckets); 
    }



  } /*end Queue namespace*/
} /*end gft namespace*/

