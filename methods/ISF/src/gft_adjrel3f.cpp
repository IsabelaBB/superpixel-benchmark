
#include "adjrel3f.h"

namespace gft{
  namespace AdjRel3f{

    AdjRel3f *Create(int n){
      AdjRel3f *A=NULL;
      
      A = (AdjRel3f *) calloc(1,sizeof(AdjRel3f));
      if (A != NULL){
	A->dx = gft::AllocFloatArray(n);
	A->dy = gft::AllocFloatArray(n);
	A->dz = gft::AllocFloatArray(n);
	A->n  = n;
      } else {
	gft::Error((char *)MSG1,(char *)"AdjRel3f::Create");
      }
      return(A);
    }


    void Destroy(AdjRel3f **A){
      AdjRel3f *aux;
      
      aux = *A;
      if (aux != NULL){
	if (aux->dx != NULL) gft::FreeFloatArray(&aux->dx);
	if (aux->dy != NULL) gft::FreeFloatArray(&aux->dy);
	if (aux->dz != NULL) gft::FreeFloatArray(&aux->dz);
	free(aux);
	*A = NULL;
      }   
    }

    
    AdjRel3f *Clone(AdjRel3f *A){
      AdjRel3f *C;
      int i;
      
      C = Create(A->n);
      for(i=0; i < A->n; i++){
	C->dx[i] = A->dx[i];
	C->dy[i] = A->dy[i];
	C->dz[i] = A->dz[i];
      }
      return C;
    }
    

    AdjRel3f *ChangeOrientationToLPS(AdjRel3f *A,
				     char *ori){
      AdjRel3f *lps=NULL;
      int i;
      
      lps = Create(A->n);
      for(i=0; i < A->n; i++){
	if     (ori[0]=='L'){ lps->dx[i] =  A->dx[i]; }
	else if(ori[0]=='R'){ lps->dx[i] = -A->dx[i]; }
	else if(ori[0]=='P'){ lps->dy[i] =  A->dx[i]; }
	else if(ori[0]=='A'){ lps->dy[i] = -A->dx[i]; }
	else if(ori[0]=='S'){ lps->dz[i] =  A->dx[i]; }
	else if(ori[0]=='I'){ lps->dz[i] = -A->dx[i]; }
	else{ gft::Error((char *)"Invalid orientation",
			 (char *)"AdjRel3f::ChangeOrientationToLPS"); }
	
	if     (ori[1]=='L'){ lps->dx[i] =  A->dy[i]; }
	else if(ori[1]=='R'){ lps->dx[i] = -A->dy[i]; }
	else if(ori[1]=='P'){ lps->dy[i] =  A->dy[i]; }
	else if(ori[1]=='A'){ lps->dy[i] = -A->dy[i]; }
	else if(ori[1]=='S'){ lps->dz[i] =  A->dy[i]; }
	else if(ori[1]=='I'){ lps->dz[i] = -A->dy[i]; }
	else{ gft::Error((char *)"Invalid orientation",
			 (char *)"AdjRel3f::ChangeOrientationToLPS"); }
	
	if     (ori[2]=='L'){ lps->dx[i] =  A->dz[i]; }
	else if(ori[2]=='R'){ lps->dx[i] = -A->dz[i]; }
	else if(ori[2]=='P'){ lps->dy[i] =  A->dz[i]; }
	else if(ori[2]=='A'){ lps->dy[i] = -A->dz[i]; }
	else if(ori[2]=='S'){ lps->dz[i] =  A->dz[i]; }
	else if(ori[2]=='I'){ lps->dz[i] = -A->dz[i]; }
	else{ gft::Error((char *)"Invalid orientation",
			 (char *)"AdjRel3f::ChangeOrientationToLPS"); }
      }
      
      return lps;
    }
    

  } /*end AdjRel3f namespace*/
} /*end gft namespace*/
