#include "gft_heap.h"

namespace gft{
  namespace Heap{

    void SetRemovalPolicy(Heap *H, char policy){
      if(H->removal_policy != policy){
	H->removal_policy = policy;
      }
    }

    void GoUp(Heap *H, int i) {
      int j = HEAP_DAD(i);
      
      if(H->removal_policy == MINVALUE){
	
	while ((i > 0) && (H->cost[H->pixel[j]] > H->cost[H->pixel[i]])) {
	  gft::SwapInt(&H->pixel[j], &H->pixel[i]);
	  H->pos[H->pixel[i]] = i;
	  H->pos[H->pixel[j]] = j;
	  i = j;
	  j = HEAP_DAD(i);
	}
      }
      else{ /* removal_policy == MAXVALUE */
	
	while ((i > 0) && (H->cost[H->pixel[j]] < H->cost[H->pixel[i]])) {
	  gft::SwapInt(&H->pixel[j], &H->pixel[i]);
	  H->pos[H->pixel[i]] = i;
	  H->pos[H->pixel[j]] = j;
	  i = j;
	  j = HEAP_DAD(i);
	}
      }
    }

    void GoDown(Heap *H, int i) {
      int j, left = HEAP_LEFTSON(i), right = HEAP_RIGHTSON(i);
      
      j = i;
      if(H->removal_policy == MINVALUE){
	
	if ((left <= H->last) && 
	    (H->cost[H->pixel[left]] < H->cost[H->pixel[i]]))
	  j = left;
	if ((right <= H->last) && 
	    (H->cost[H->pixel[right]] < H->cost[H->pixel[j]]))
	  j = right;
      }
      else{ /* removal_policy == MAXVALUE */
	
	if ((left <= H->last) && 
	    (H->cost[H->pixel[left]] > H->cost[H->pixel[i]]))
	  j = left;
	if ((right <= H->last) && 
	    (H->cost[H->pixel[right]] > H->cost[H->pixel[j]]))
	  j = right;
      }
      
      if(j != i) {
	gft::SwapInt(&H->pixel[j], &H->pixel[i]);
	H->pos[H->pixel[i]] = i;
	H->pos[H->pixel[j]] = j;
	GoDown(H, j);
      }
    }

    char IsFull(Heap *H) {
      if (H->last == (H->n - 1))
	return 1;
      else
	return 0;
    }
    
    char IsEmpty(Heap *H) {
      if (H->last == -1){
	//Reset(H);
	return 1;
      }else
	return 0;
    }
    
    Heap *Create(int n, float *cost) {
      Heap *H = NULL;
      int i;
      
      if (cost == NULL) {
	fprintf(stdout,"Cannot create heap without cost map in Heap::Create");
	return NULL;
      }
      
      H = (Heap *) malloc(sizeof(Heap));
      if (H != NULL) {
	H->n       = n;
	H->cost    = cost;
	H->color   = (char *) malloc(sizeof(char) * n);
	H->pixel   = (int *) malloc(sizeof(int) * n);
	H->pos     = (int *) malloc(sizeof(int) * n);
	H->last    = -1;
	H->removal_policy = MINVALUE;
	if (H->color == NULL || H->pos == NULL || H->pixel == NULL)
	  gft::Error((char *)MSG1,(char *)"Heap::Create");
	for (i = 0; i < H->n; i++) {
	  H->color[i] = WHITE;
	  H->pos[i]   = -1;
	  H->pixel[i] = -1;
	}    
      } 
      else
	gft::Error((char *)MSG1,(char *)"Heap::Create");

      return H;
    }

    void Destroy(Heap **H) {
      Heap *aux = *H;
      if (aux != NULL) {
	if (aux->pixel != NULL) free(aux->pixel);
	if (aux->color != NULL) free(aux->color);
	if (aux->pos != NULL)   free(aux->pos);
	free(aux);
	*H = NULL;
      }
    }
    
    char Insert(Heap *H, int pixel) {
      if (!IsFull(H)) {
	H->last++;
	H->pixel[H->last] = pixel;
	H->color[pixel]   = GRAY;
	H->pos[pixel]     = H->last;
	GoUp(H, H->last); 
	return 1;
      } else 
	return 0;
    }

    char Remove(Heap *H, int *pixel) {
      if (!IsEmpty(H)) {
	*pixel = H->pixel[0];
	H->pos[*pixel]   = -1;
	H->color[*pixel] = BLACK;
	H->pixel[0]      = H->pixel[H->last];
	H->pos[H->pixel[0]] = 0;
	H->pixel[H->last] = -1;
	H->last--;
	GoDown(H, 0);
	return 1;
      } else 
	return 0;
    }

    /* Quando o novo custo eh pior.
    void Update(Heap *H, int p, float value) {
	bool up = value < H->cost[p] ? true : false ;
	
	H->cost[p] = value;

	if (H->color[p] == BLACK) printf("ferrou\n");

	if (H->color[p] == WHITE)
		Insert(H, p);
	else {
		if (up) {
			GoUp(H, H->pos[p]);
			
		}
		else {
			GoDown(H, H->pos[p]);
			
		}
	}
   }
    */

    void Update(Heap *H, int p, float value){
      H->cost[p] = value;
      
      if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert(H, p);
      else
	GoUp(H, H->pos[p]);
    }
    
    void Reset(Heap *H){
      int i;
      
      for (i=0; i < H->n; i++) {
	H->color[i] = WHITE;
	H->pos[i]   = -1;
	H->pixel[i] = -1;
      }
      H->last = -1;
    }


char Delete	(Heap *H, int pixel) {
	int i, t = 0;
	if (!IsEmpty(H)) {
		Update(H,pixel,-2);
		Remove(H, &i);
		H->color[pixel] = WHITE;
		return 1;
	}
	else
		return 0;
}


    

  } /*end Heap namespace*/
} /*end gft namespace*/


