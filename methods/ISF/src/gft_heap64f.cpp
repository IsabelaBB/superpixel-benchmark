#include "gft_heap64f.h"

namespace gft{
  namespace Heap64f{

    void GoUp_MaxPolicy(Heap64f *H, int i){
      int j = HEAP_DAD(i);
      int p = H->pixel[i];
      double c = H->cost[p];
      while((i > 1)&&(H->cost[H->pixel[j]] < c)){
	H->pixel[i] = H->pixel[j];
	H->pos[H->pixel[j]] = i;
	i = j;
	j = HEAP_DAD(i);
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }


    void GoUp_MinPolicy(Heap64f *H, int i){
      int j = HEAP_DAD(i);
      int p = H->pixel[i];
      double c = H->cost[p];
      while((i > 1)&&(H->cost[H->pixel[j]] > c)){
	H->pixel[i] = H->pixel[j];
	H->pos[H->pixel[j]] = i;
	i = j;
	j = HEAP_DAD(i);
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }    

    

    void GoDown_MaxPolicy(Heap64f *H, int i) {
      int j, left, right, p = H->pixel[i];
      double c = H->cost[p];
      j = i;
      while(1){
	i = j;
	left = HEAP_LEFTSON(i);
	right = left + 1;
	if(right <= H->last){
	  if(H->cost[H->pixel[right]] > H->cost[H->pixel[left]])
	    j = right;
	  else
	    j = left;
	}
	else if(left <= H->last)
	  j = left;
	else
	  break;
	
	if(H->cost[H->pixel[j]] > c) {
	  H->pixel[i] = H->pixel[j];
	  H->pos[H->pixel[j]] = i;
	}
	else break;
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }
    

    void GoDown_MinPolicy(Heap64f *H, int i) {
      int j, left, right, p = H->pixel[i];
      double c = H->cost[p];
      j = i;
      while(1){
	i = j;
	left = HEAP_LEFTSON(i);
	right = left + 1;
	if(right <= H->last){
	  if(H->cost[H->pixel[right]] < H->cost[H->pixel[left]])
	    j = right;
	  else
	    j = left;
	}
	else if(left <= H->last)
	  j = left;
	else
	  break;
	
	if(H->cost[H->pixel[j]] < c) {
	  H->pixel[i] = H->pixel[j];
	  H->pos[H->pixel[j]] = i;
	}
	else break;
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }

    
    char IsFull(Heap64f *H) {
      if (H->last == H->n)
	return 1;
      else
	return 0;
    }
    
    char IsEmpty(Heap64f *H) {
      if (H->last == 0)
	return 1;
      else
	return 0;
    }

    
    Heap64f *Create(int n, double *cost) {
      Heap64f *H = NULL;
      int i;
      
      if (cost == NULL) {
	fprintf(stdout,"Cannot create heap without cost map in Heap64f::Create");
	return NULL;
      }
      
      H = (Heap64f *) malloc(sizeof(Heap64f));
      if (H != NULL) {
	H->n       = n;
	H->cost    = cost;
	H->color   = (char *) malloc(sizeof(char) * n);
	H->pixel   = (int *) malloc(sizeof(int) * (n+1));
	H->pos     = (int *) malloc(sizeof(int) * n);
	H->last    = 0;
	if (H->color == NULL || H->pos == NULL || H->pixel == NULL)
	  gft::Error((char *)MSG1,(char *)"Heap64f::Create");
	for (i = 0; i < H->n; i++) {
	  H->color[i] = WHITE;
	  H->pos[i]   = -1;
	  H->pixel[i] = -1;
	}
	H->pixel[n] = -1;
      }
      else
	gft::Error((char *)MSG1,(char *)"Heap64f::Create");

      return H;
    }

    void Destroy(Heap64f **H) {
      Heap64f *aux = *H;
      if (aux != NULL) {
	if (aux->pixel != NULL) free(aux->pixel);
	if (aux->color != NULL) free(aux->color);
	if (aux->pos != NULL)   free(aux->pos);
	free(aux);
	*H = NULL;
      }
    }
    
    void Insert_MaxPolicy(Heap64f *H, int pixel) {
      H->last++;
      H->pixel[H->last] = pixel;
      H->color[pixel]   = GRAY;
      H->pos[pixel]     = H->last;
      GoUp_MaxPolicy(H, H->last);
    }

    
    void Insert_MinPolicy(Heap64f *H, int pixel) {
      H->last++;
      H->pixel[H->last] = pixel;
      H->color[pixel]   = GRAY;
      H->pos[pixel]     = H->last;
      GoUp_MinPolicy(H, H->last);
    }

    
    void Remove_MaxPolicy(Heap64f *H, int *pixel) {
      *pixel = H->pixel[1];
      H->pos[*pixel]   = -1;
      H->color[*pixel] = BLACK;
      if(H->last == 1){
	H->pixel[1] = -1; //Linha nao obrigatoria.
	H->last = 0;
      }
      else if(H->last > 1){
	H->pixel[1]      = H->pixel[H->last];
	H->pos[H->pixel[1]] = 1;
	H->pixel[H->last] = -1;
	H->last--;
	GoDown_MaxPolicy(H, 1);
      }
    }


    void Remove_MinPolicy(Heap64f *H, int *pixel) {
      *pixel = H->pixel[1];
      H->pos[*pixel]   = -1;
      H->color[*pixel] = BLACK;
      if(H->last == 1){
	H->pixel[1] = -1; //Linha nao obrigatoria.
	H->last = 0;
      }
      else if(H->last > 1){
	H->pixel[1]      = H->pixel[H->last];
	H->pos[H->pixel[1]] = 1;
	H->pixel[H->last] = -1;
	H->last--;
	GoDown_MinPolicy(H, 1);
      }
    }


    void Update_MaxPolicy(Heap64f *H, int p, double value){
      H->cost[p] = value;
      
      //if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert_MaxPolicy(H, p);
      else
	GoUp_MaxPolicy(H, H->pos[p]);
    }


    void Update_MinPolicy(Heap64f *H, int p, double value){
      H->cost[p] = value;
      
      //if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert_MinPolicy(H, p);
      else
	GoUp_MinPolicy(H, H->pos[p]);
    }
    
    
    void Reset(Heap64f *H){
      int i;
      
      for (i=0; i < H->n; i++) {
	H->color[i] = WHITE;
	H->pos[i]   = -1;
	H->pixel[i] = -1;
      }
      H->pixel[H->n] = -1;
      H->last = 0;
    }


    void Delete_MaxPolicy(Heap64f *H, int pixel) {
      int q;
      //if (IsEmpty(H)) printf("error: delete pixel%d\n", pixel);
      q = H->pixel[H->last];
      H->pixel[H->pos[pixel]] = q;
      H->pos[q] = H->pos[pixel];
      H->pixel[H->last] = -1;
      H->last--;
      H->pos[pixel] = -1;
      H->color[pixel] = WHITE;
      if(H->cost[pixel] > H->cost[q])
	GoDown_MaxPolicy(H, H->pos[q]);
      else 
	GoUp_MaxPolicy(H, H->pos[q]);
    }
    

    void Delete_MinPolicy(Heap64f *H, int pixel) {
      int q;
      //if (IsEmpty(H)) printf("error: delete pixel%d\n", pixel);
      q = H->pixel[H->last];
      H->pixel[H->pos[pixel]] = q;
      H->pos[q] = H->pos[pixel];
      H->pixel[H->last] = -1;
      H->last--;
      H->pos[pixel] = -1;
      H->color[pixel] = WHITE;
      if(H->cost[pixel] < H->cost[q])
	GoDown_MinPolicy(H, H->pos[q]);
      else 
	GoUp_MinPolicy(H, H->pos[q]);
    }

    

  } /*end Heap64f namespace*/
} /*end gft namespace*/


