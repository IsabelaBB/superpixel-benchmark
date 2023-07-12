
#include "gft_pqueue32.h"

namespace gft{
  namespace PQueue32{


    PQueue32 *Create(int nbuckets, int nelems, int *value){
      PQueue32 *Q=NULL;

      Q = (PQueue32 *) malloc(1*sizeof(PQueue32));
  
      if (Q != NULL) {
	Q->C.first = (int *)malloc((nbuckets+1) * sizeof(int));
	Q->C.last  = (int *)malloc((nbuckets+1) * sizeof(int));
	Q->C.nbuckets = nbuckets;
	if ( (Q->C.first != NULL) && (Q->C.last != NULL) ){
	  Q->L.elem = (PQNode *)malloc(nelems*sizeof(PQNode));
	  Q->L.nelems = nelems;
	  Q->L.value  = value;
	  if (Q->L.elem != NULL){
	    Reset(Q);
	  } else
	    gft::Error((char *)MSG1,
		       (char *)"PQueue32::Create");
	} else
	  gft::Error((char *)MSG1,
		     (char *)"PQueue32::Create");
      } else 
	gft::Error((char *)MSG1,
		   (char *)"PQueue32::Create");
      return(Q);
    }


    void Destroy(PQueue32 **Q){
      PQueue32 *aux;
      aux = *Q;
      if (aux != NULL) {
	if (aux->C.first != NULL) free(aux->C.first);
	if (aux->C.last  != NULL) free(aux->C.last);
	if (aux->L.elem  != NULL) free(aux->L.elem);
	free(aux);
	*Q = NULL;
      }
    }


    PQueue32 *Grow(PQueue32 **Q, int nbuckets){
      PQueue32 *Q1;
      int i,bucket;
      
      Q1 = Create(nbuckets,(*Q)->L.nelems,(*Q)->L.value);
      Q1->nadded = (*Q)->nadded;
      Q1->C.minvalue  = (*Q)->C.minvalue;
      Q1->C.maxvalue  = (*Q)->C.maxvalue;
      for (i=0; i<(*Q)->C.nbuckets; i++){
	if ((*Q)->C.first[i]!=NIL){
	  bucket = (*Q)->L.value[(*Q)->C.first[i]]%Q1->C.nbuckets;
	  Q1->C.first[bucket] = (*Q)->C.first[i];
	  Q1->C.last[bucket]  = (*Q)->C.last[i];
	}
      }
      if ((*Q)->C.first[(*Q)->C.nbuckets]!=NIL){
	bucket = Q1->C.nbuckets;
	Q1->C.first[bucket] = (*Q)->C.first[(*Q)->C.nbuckets];
	Q1->C.last[bucket]  = (*Q)->C.last[(*Q)->C.nbuckets];
      }
      for (i=0; i < (*Q)->L.nelems; i++) 
	Q1->L.elem[i] = (*Q)->L.elem[i];
      
      Destroy(Q);
      return(Q1);
    }
    

    void Reset(PQueue32 *Q){
      int i;
      Q->nadded = 0;
      Q->C.minvalue = INT_MAX;
      Q->C.maxvalue = INT_MIN;
      for (i=0; i < Q->C.nbuckets+1; i++)
	Q->C.first[i]=Q->C.last[i]=NIL;
    
      for (i=0; i < Q->L.nelems; i++) {
	Q->L.elem[i].next =  Q->L.elem[i].prev = NIL;
	Q->L.elem[i].color = WHITE;
      }
    }



    // Generic version with circular and growing features
    void   InsertElem(PQueue32 **Q, int elem){
      int bucket,value,dvalue;
      
      (*Q)->nadded++;
      value = (*Q)->L.value[elem];
      if(value==INT_MAX || value==INT_MIN)
	bucket=(*Q)->C.nbuckets;
      else{
	if(value < (*Q)->C.minvalue)
	  (*Q)->C.minvalue = value;
	if(value > (*Q)->C.maxvalue)
	  (*Q)->C.maxvalue = value;
	
	dvalue = (*Q)->C.maxvalue - (*Q)->C.minvalue;
	if (dvalue > ((*Q)->C.nbuckets-1)){
	  (*Q) = Grow(Q, 2*(dvalue)+1);
	  gft::Warning((char *)"Doubling queue size",
		       (char *)"PQueue32::InsertElem");
	}
	bucket = value%(*Q)->C.nbuckets;
      }
      if ((*Q)->C.first[bucket] == NIL){ 
	(*Q)->C.first[bucket]   = elem;  
	(*Q)->L.elem[elem].prev = NIL;
      }else {
	(*Q)->L.elem[(*Q)->C.last[bucket]].next = elem;
	(*Q)->L.elem[elem].prev = (*Q)->C.last[bucket];
      }
      
      (*Q)->C.last[bucket]     = elem;
      (*Q)->L.elem[elem].next  = NIL;
      (*Q)->L.elem[elem].color = GRAY;
    }


    void   RemoveElem(PQueue32 *Q, int elem){
      int prev,next,bucket;
  
      if(Q->L.elem[elem].color!=GRAY) return;
      
      Q->nadded--;
      if ((Q->L.value[elem]==INT_MAX)||(Q->L.value[elem]==INT_MIN))
	bucket = Q->C.nbuckets;
      else
	bucket = Q->L.value[elem]%Q->C.nbuckets;
      
      prev = Q->L.elem[elem].prev;
      next = Q->L.elem[elem].next;
  
      // if elem is the first element
      if (Q->C.first[bucket] == elem) {
	Q->C.first[bucket] = next;
	if (next == NIL) // elem is also the last one
	  Q->C.last[bucket] = NIL;
	else
	  Q->L.elem[next].prev = NIL;
      }
      else{   // elem is in the middle or it is the last
	Q->L.elem[prev].next = next;
	if (next == NIL) // if it is the last
	  Q->C.last[bucket] = prev;
	else 
	  Q->L.elem[next].prev = prev;
      }
      Q->L.elem[elem].color = BLACK;
    }


    void   UpdateElem(PQueue32 **Q, int elem, int newvalue){
      RemoveElem(*Q, elem);
      (*Q)->L.value[elem] = newvalue;
      InsertElem(Q, elem);
    }
    

    /// \cond
    inline int _FindMinBucket(PQueue32 *Q){
      int current,last;
      
      current = Q->C.minvalue%Q->C.nbuckets;
      // moves to next element
      if(Q->C.first[current] == NIL){
	last = current;
	
	do{
	  current = (current + 1) % (Q->C.nbuckets);
	}while((Q->C.first[current] == NIL) && (current != last));
	
	if(Q->C.first[current] != NIL)
	  Q->C.minvalue = Q->L.value[Q->C.first[current]];
	else{
	  if(Q->C.first[Q->C.nbuckets] != NIL){
	    current = Q->C.nbuckets;
	    Q->C.minvalue = Q->L.value[Q->C.first[current]];
	  }
	  else
	    gft::Error((char *)"PQueue is empty",
		       (char *)"_FindMinBucket");
	}
      }
      return current;
    }
    /// \endcond
    

    /// \cond    
    inline int _FindMaxBucket(PQueue32 *Q){
      int current,last;
      
      current = Q->C.maxvalue%Q->C.nbuckets;
      // moves to next element
      if(Q->C.first[current] == NIL){
	last = current;
	
	do{
	  current--;
	  if(current<0) current = Q->C.nbuckets-1;
	}while((Q->C.first[current] == NIL) && (current != last));
	
	if(Q->C.first[current] != NIL)
	  Q->C.maxvalue = Q->L.value[Q->C.first[current]];
	else{
	  if(Q->C.first[Q->C.nbuckets] != NIL){
	    current = Q->C.nbuckets;
	    Q->C.maxvalue = Q->L.value[Q->C.first[current]];
	  }
	  else
	    gft::Error((char *)"PQueue is empty",
		       (char *)"_FindMaxBucket");
	}
      }
      return current;
    }
    /// \endcond


    /// \cond
    inline int _BucketFIFO(PQueue32 *Q, int bucket){
      int elem=NIL, next;
      
      elem = Q->C.first[bucket];
      next = Q->L.elem[elem].next;
      if(next == NIL) { // there was a single element in the list
	Q->C.first[bucket] = Q->C.last[bucket] = NIL;
      }
      else {
	Q->C.first[bucket] = next;
	Q->L.elem[next].prev = NIL;
      }
      Q->L.elem[elem].color = BLACK;
      return elem;
    }
    /// \endcond


    /// \cond
    inline int _BucketLIFO(PQueue32 *Q, int bucket){
      int elem=NIL, prev;
      
      elem = Q->C.last[bucket];
      prev = Q->L.elem[elem].prev;
      if(prev == NIL){ // there was a single element in the list
	Q->C.last[bucket] = Q->C.first[bucket] = NIL;
      }
      else {
	Q->C.last[bucket] = prev;
	Q->L.elem[prev].next = NIL;
      }
      Q->L.elem[elem].color = BLACK;
      return elem;
    }
    /// \endcond


    int    RemoveMinFIFO(PQueue32 *Q){
      int bucket;
      Q->nadded--;
      bucket = _FindMinBucket(Q);
      return _BucketFIFO(Q, bucket);
    }
    

    int    RemoveMinLIFO(PQueue32 *Q){
      int bucket;
      Q->nadded--;
      bucket = _FindMinBucket(Q);
      return _BucketLIFO(Q, bucket);
    }


    int    RemoveMaxFIFO(PQueue32 *Q){
      int bucket;
      Q->nadded--;
      bucket = _FindMaxBucket(Q);
      return _BucketFIFO(Q, bucket);
    }


    int    RemoveMaxLIFO(PQueue32 *Q){
      int bucket;
      Q->nadded--;
      bucket = _FindMaxBucket(Q);
      return _BucketLIFO(Q, bucket);
    }


    void FastInsertElem(PQueue32 *Q, int elem){
      int bucket;
      Q->nadded++;
      bucket = Q->L.value[elem];
      
      if(bucket < Q->C.minvalue)
	Q->C.minvalue = bucket;
      if(bucket > Q->C.maxvalue)
	Q->C.maxvalue = bucket;

      if(Q->C.first[bucket] == NIL){ 
	Q->C.first[bucket]   = elem;  
	Q->L.elem[elem].prev = NIL;
      }else {
	Q->L.elem[Q->C.last[bucket]].next = elem;
	Q->L.elem[elem].prev = Q->C.last[bucket];
      }
      Q->C.last[bucket]     = elem;
      Q->L.elem[elem].next  = NIL;
      Q->L.elem[elem].color = GRAY;
    }

    void FastInsertElemAsFirst(PQueue32 *Q, int elem){
      int bucket;
      
      Q->nadded++;
      bucket = Q->L.value[elem];
      if(bucket < Q->C.minvalue)
	Q->C.minvalue = bucket;
      if(bucket > Q->C.maxvalue)
	Q->C.maxvalue = bucket;
      
      if(Q->C.first[bucket] == NIL){ 
	Q->C.last[bucket]    = elem;
	Q->L.elem[elem].next = NIL;
      }
      else{
	Q->L.elem[Q->C.first[bucket]].prev = elem;
	Q->L.elem[elem].next = Q->C.first[bucket];
      }
      Q->C.first[bucket]   = elem;  
      Q->L.elem[elem].prev = NIL;
      Q->L.elem[elem].color = GRAY;
    }
    

    void FastRemoveElem(PQueue32 *Q, int elem){
      int prev,next,bucket;
      Q->nadded--;
      bucket = Q->L.value[elem];
      prev = Q->L.elem[elem].prev;
      next = Q->L.elem[elem].next;
      // if elem is the first element
      if (Q->C.first[bucket] == elem) {
	Q->C.first[bucket] = next;
	if (next == NIL) // elem is also the last one
	  Q->C.last[bucket] = NIL;
	else
	  Q->L.elem[next].prev = NIL;
      }
      else{   // elem is in the middle or it is the last
	Q->L.elem[prev].next = next;
	if (next == NIL) // if it is the last
	  Q->C.last[bucket] = prev;
	else 
	  Q->L.elem[next].prev = prev;
      }
      Q->L.elem[elem].color = BLACK;
    }


    void FastUpdateElem(PQueue32 *Q, int elem, int newvalue){
      FastRemoveElem(Q, elem);
      Q->L.value[elem] = newvalue;
      FastInsertElem(Q, elem);
    }
    

    /// \cond
    inline int _FastFindMinBucket(PQueue32 *Q){
      int current;
      current = Q->C.minvalue;
      // moves to next element
      if(Q->C.first[current] == NIL){
	do{
	  current++;
	}while((current<Q->C.nbuckets)&&(Q->C.first[current] == NIL));
	
	if(current < Q->C.nbuckets)
	  Q->C.minvalue = current;
	else
	  gft::Error((char *)"PQueue is empty",
		     (char *)"_FastFindMinBucket");
      }
      return current;
    }
    /// \endcond


    /// \cond
    inline int _FastFindMaxBucket(PQueue32 *Q){
      int current;
      
      current = Q->C.maxvalue;
      // moves to next element
      if(Q->C.first[current] == NIL){
	do{
	  current--;
	}while((current>=0)&&(Q->C.first[current] == NIL));
	
	if(current >= 0)
	  Q->C.maxvalue = current;
	else
	  gft::Error((char *)"PQueue is empty",
		     (char *)"_FastFindMaxBucket");
      }
      return current;
    }
    /// \endcond


    int    FastRemoveMinFIFO(PQueue32 *Q){
      int bucket;
      Q->nadded--;
      bucket = _FastFindMinBucket(Q);
      return _BucketFIFO(Q, bucket);
    }
    
    int    FastRemoveMinLIFO(PQueue32 *Q){
      int bucket;
      Q->nadded--;
      bucket = _FastFindMinBucket(Q);
      return _BucketLIFO(Q, bucket);
    }

    int    FastRemoveMaxFIFO(PQueue32 *Q){
      int bucket;
      Q->nadded--;
      bucket = _FastFindMaxBucket(Q);
      return _BucketFIFO(Q, bucket);
    }

    int    FastRemoveMaxLIFO(PQueue32 *Q){
      int bucket;
      Q->nadded--;
      bucket = _FastFindMaxBucket(Q);
      return _BucketLIFO(Q, bucket);
    }


  } //end PQueue32 namespace
} //end gft namespace

