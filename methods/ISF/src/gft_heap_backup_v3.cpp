#include "gft_heap.h"

namespace gft{
  namespace Heap{


    void GoUp_MaxPolicy(Heap *H, int i){
      int j = HEAP_DAD(i);
      int p = H->pixel[i];
      float c = H->cost[p];
      while(H->cost[H->pixel[j]] < c){
	H->pixel[i] = H->pixel[j];
	H->pos[H->pixel[j]] = i;
	i = j;
	j = HEAP_DAD(i);
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }


    void GoUp_MinPolicy(Heap *H, int i){
      int j = HEAP_DAD(i);
      int p = H->pixel[i];
      float c = H->cost[p];
      while(H->cost[H->pixel[j]] > c){
	H->pixel[i] = H->pixel[j];
	H->pos[H->pixel[j]] = i;
	i = j;
	j = HEAP_DAD(i);
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }    


    void GoDown_MaxPolicy(Heap *H, int i) {
      int j, left, right, p = H->pixel[i];
      float c = H->cost[p];
      j = i;
      while(1){
	i = j;
	left = HEAP_LEFTSON(i);
	right = left + 1;
	if(H->cost[H->pixel[right]] > H->cost[H->pixel[left]])
	  j = right;
	else
	  j = left;
	
	if(H->cost[H->pixel[j]] > c) {
	  H->pixel[i] = H->pixel[j];
	  H->pos[H->pixel[j]] = i;
	}
	else break;
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }
    

    void GoDown_MinPolicy(Heap *H, int i) {
      int j, left, right, p = H->pixel[i];
      float c = H->cost[p];
      j = i;
      while(1){
	i = j;
	left = HEAP_LEFTSON(i);
	right = left + 1;
	if(H->cost[H->pixel[right]] < H->cost[H->pixel[left]])
	  j = right;
	else
	  j = left;
	
	if(H->cost[H->pixel[j]] < c) {
	  H->pixel[i] = H->pixel[j];
	  H->pos[H->pixel[j]] = i;
	}
	else break;
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }

    
    char IsFull(Heap *H) {
      if (H->last == H->n)
	return 1;
      else
	return 0;
    }
    
    char IsEmpty(Heap *H) {
      if (H->last == 0)
	return 1;
      else
	return 0;
    }


    Heap *Create_MaxPolicy(int n, float **cost) {
      Heap *H = Create_MinPolicy(n, cost);
      H->cost[n] = -FLT_MAX;
      H->cost[n+1] = FLT_MAX;
      return H;
    }

    
    Heap *Create_MinPolicy(int n, float **cost) {
      Heap *H = NULL;
      int i;
      
      if (cost == NULL) {
	fprintf(stdout,"Cannot create heap without cost map in Heap::Create");
	return NULL;
      }
      
      H = (Heap *) malloc(sizeof(Heap));
      if (H != NULL) {
	H->n       = n;
	H->cost    = gft::AllocFloatArray(n+2); 
	H->color   = (char *) malloc(sizeof(char) * n);
	H->pixel   = (int *) malloc(sizeof(int) * (n*2+2));
	H->pos     = (int *) malloc(sizeof(int) * n);
	H->last    = 0;
	if (H->color == NULL || H->pos == NULL || H->pixel == NULL)
	  gft::Error((char *)MSG1,(char *)"Heap::Create");
	for (i = 0; i < H->n; i++) {
	  H->color[i] = WHITE;
	  H->pos[i]   = -1;
	}
	H->pixel[0] = n+1;
	for (i = 1; i < H->n*2+2; i++) {
	  H->pixel[i] = n;
	}
	H->cost[n] = FLT_MAX;
	H->cost[n+1] = -FLT_MAX;
	*cost = H->cost;
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
    
    char Insert_MaxPolicy(Heap *H, int pixel) {
      //if (!IsFull(H)) {
      if(H->last < H->n){
	H->last++;
	H->pixel[H->last] = pixel;
	H->color[pixel]   = GRAY;
	H->pos[pixel]     = H->last;
	GoUp_MaxPolicy(H, H->last);
	return 1;
      } else 
	return 0;
    }

    
    char Insert_MinPolicy(Heap *H, int pixel) {
      //if (!IsFull(H)) {
      if(H->last < H->n){
	H->last++;
	H->pixel[H->last] = pixel;
	H->color[pixel]   = GRAY;
	H->pos[pixel]     = H->last;
	GoUp_MinPolicy(H, H->last);
	return 1;
      } else 
	return 0;
    }

    
    char Remove_MaxPolicy(Heap *H, int *pixel) {
      if(H->last == 1){
	*pixel = H->pixel[1];
	H->pos[*pixel]   = -1;
	H->color[*pixel] = BLACK;
	H->pixel[1] = H->n; 
	H->last = 0;
	return 1;
      }
      else if(H->last > 1){
	*pixel = H->pixel[1];
	H->pos[*pixel]   = -1;
	H->color[*pixel] = BLACK;
	H->pixel[1]      = H->pixel[H->last];
	H->pos[H->pixel[1]] = 1;
	H->pixel[H->last] = H->n;
	H->last--;
	GoDown_MaxPolicy(H, 1);
	return 1;
      } else 
	return 0;
    }


    char Remove_MinPolicy(Heap *H, int *pixel) {
      if(H->last == 1){
	*pixel = H->pixel[1];
	H->pos[*pixel]   = -1;
	H->color[*pixel] = BLACK;
	H->pixel[1] = H->n;
	H->last = 0;
	return 1;
      }
      else if(H->last > 1){
	*pixel = H->pixel[1];
	H->pos[*pixel]   = -1;
	H->color[*pixel] = BLACK;
	H->pixel[1]      = H->pixel[H->last];
	H->pos[H->pixel[1]] = 1;
	H->pixel[H->last] = H->n;
	H->last--;
	GoDown_MinPolicy(H, 1);
	return 1;
      } else 
	return 0;
    }


    void Update_MaxPolicy(Heap *H, int p, float value){
      H->cost[p] = value;
      
      if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert_MaxPolicy(H, p);
      else
	GoUp_MaxPolicy(H, H->pos[p]);
    }


    void Update_MinPolicy(Heap *H, int p, float value){
      H->cost[p] = value;
      
      if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert_MinPolicy(H, p);
      else
	GoUp_MinPolicy(H, H->pos[p]);
    }
    
    
    void Reset(Heap *H){
      int i;
      
      for (i=0; i < H->n; i++) {
	H->color[i] = WHITE;
	H->pos[i]   = -1;
      }
      H->pixel[0] = H->n+1;
      for (i = 1; i < H->n*2+2; i++) {
	H->pixel[i] = H->n;
      }
      H->last = 0;
    }


    char Delete_MaxPolicy(Heap *H, int pixel) {
      int i, t = 0;
      if (!IsEmpty(H)) {
	Update_MaxPolicy(H,pixel,FLT_MAX);
	Remove_MaxPolicy(H, &i);
	H->color[pixel] = WHITE;
	return 1;
      }
      else
	return 0;
    }
    

    char Delete_MinPolicy(Heap *H, int pixel) {
      int i, t = 0;
      if (!IsEmpty(H)) {
	Update_MinPolicy(H,pixel,-FLT_MAX);
	Remove_MinPolicy(H, &i);
	H->color[pixel] = WHITE;
	return 1;
      }
      else
	return 0;
    }
  

  } /*end Heap namespace*/
} /*end gft namespace*/


