
#include "gft_common.h"

namespace gft{
  
  char *AllocCharArray(int n){
    char *v=NULL;
    v = (char *) calloc(n,sizeof(char));
    if (v == NULL)
      Error((char *)MSG1,(char *)"gft::AllocCharArray");
    return(v);
  }

  /*  
  uchar *AllocUCharArray(int n){
    uchar *v=NULL;
    v = (uchar *) calloc(n,sizeof(uchar));
    if (v == NULL)
      Error(MSG1,"gft::AllocUCharArray");
    return(v);
  }
  */
  uchar  *AllocUCharArray(int n){
    uchar *v=NULL;
    
    v = (uchar *)_mm_malloc(sizeof(uchar)*n, 16);
    if(v==NULL)
      Error((char *)MSG1,(char *)"gft::AllocUCharArray");
    memset((void *)v, 0, sizeof(uchar)*n);
    return(v);
  }

  /*
  ushort *AllocUShortArray(int n){
    ushort *v=NULL;
    v = (ushort *) calloc(n,sizeof(ushort));
    if (v == NULL)
      Error(MSG1,"gft::AllocUShortArray");
    return(v);
  }
  */
  ushort *AllocUShortArray(int n){
    ushort *v=NULL;
    
    v = (ushort *)_mm_malloc(sizeof(ushort)*n, 16);
    if(v==NULL)
      Error((char *)MSG1,(char *)"gft::AllocUShortArray");
    memset((void *)v, 0, sizeof(ushort)*n);
    return(v);
  }


  uint *AllocUIntArray(int n){
    uint *v=NULL;
    v = (uint *) calloc(n,sizeof(uint));
    if (v == NULL)
      Error((char *)MSG1,(char *)"gft::AllocUIntArray");
    return(v);
  }
  
  /*
  int *AllocIntArray(int n){
    int *v=NULL;
    v = (int *) calloc(n,sizeof(int));
    if (v == NULL)
      Error(MSG1,"gft::AllocIntArray");
    return(v);
  }
  */  
  int *AllocIntArray(int n){
    int *v=NULL;

    v = (int *)_mm_malloc(sizeof(int)*n, sizeof(int)*4);
    if(v==NULL)
      Error((char *)MSG1,(char *)"gft::AllocIntArray");
    memset((void *)v, 0, sizeof(int)*n);
    return(v);
  }


  long long *AllocLongLongArray(int n){
    long long *v=NULL;
    v = (long long *) calloc(n,sizeof(long long));
    if (v == NULL)
      Error((char *)MSG1,(char *)"gft::AllocLongLongArray");
    return(v);
  }

  
  /*
  float *AllocFloatArray(int n){
    float *v=NULL;
    v = (float *) calloc(n,sizeof(float));
    if (v == NULL)
      Error(MSG1,"gft::AllocFloatArray");
    return(v);
  }
  */
  float *AllocFloatArray(int n){
    float *v=NULL;

    //if(posix_memalign((void **)&v, 16, sizeof(float)*n)!=0){
    //  printf("memalign error\n");
    //  exit(1);
    //}
    v = (float *)_mm_malloc(sizeof(float)*n, sizeof(float)*4);
    if(v==NULL) 
      Error((char *)MSG1,(char *)"gft::AllocFloatArray");
    memset((void *)v, 0, sizeof(float)*n);
    return(v);
  }


  double *AllocDoubleArray(int n){
    double *v=NULL;
    v = (double *) calloc(n,sizeof(double));
    if (v == NULL)
      Error((char *)MSG1,(char *)"gft::AllocDoubleArray");
    return(v);
  }

  void    FreeIntArray(int **a){
    _mm_free(*a);
    *a = NULL;
  }

  void    FreeFloatArray(float **a){
    _mm_free(*a);
    *a = NULL;
  }

  void    FreeDoubleArray(double **a){
    free(*a);
    *a = NULL;
  }
  
  void    FreeUCharArray(uchar **a){
    _mm_free(*a);
    *a = NULL;
  }

  void    FreeUShortArray(ushort **a){
    _mm_free(*a);
    *a = NULL;
  }
  
  void Error(char *msg,char *func){ /* It prints error message and exits
				       the program. */
    fprintf(stderr,"Error:%s in %s\n",msg,func);
    exit(-1);
  }
  
  void Warning(char *msg,char *func){ /* It prints warning message and
					 leaves the routine. */
    fprintf(stdout,"Warning:%s in %s\n",msg,func);
  }


  int NCFgets(char *s, int m, FILE *f) {
    while(fgets(s,m,f)!=NULL)
      if (s[0]!='#') return 1;
    return 0;
  }
  

  /**
   * Gera um número inteiro aleatório no intervalo [low,high].
   http://www.ime.usp.br/~pf/algoritmos/aulas/random.html
   Execute RandomSeed antes, para atualizar a semente.
  **/
  int RandomInteger (int low, int high){
    int k;
    double d;
    
    d = (double) rand () / ((double) RAND_MAX + 1);
    k = d * (high - low + 1);
    return low + k;
  }
  
  void RandomSeed() { 
    unsigned int t=0;
    timer *T;
    
    T = (timer *)malloc(sizeof(timer));
    gettimeofday(T,NULL); 

    if (T!=NULL){
      t = (unsigned int)(T->tv_usec*0.001);
      srand(t);
      free(T);
    }
  }



} /*end gft namespace*/
