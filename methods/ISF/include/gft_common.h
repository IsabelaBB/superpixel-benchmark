/**
 * \file gft_common.h
 * \brief Header file for common definitions and function prototypes.
 */

#ifndef _GFT_COMMON_H_
#define _GFT_COMMON_H_

extern "C" {
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
  
#include <xmmintrin.h>
}

/**
 * \brief Base namespace for common definitions and prototypes.
 */
namespace gft{

/* Error messages */
#define MSG1  "Cannot allocate memory space"
#define MSG2  "Cannot open file"
#define MSG3  "Invalid option"
#define MSG4  "Could not locate nifti header"
#define MSG5  "Nifti-1 data type not supported"

/* Common data types to all programs */
#ifndef __cplusplus
typedef enum boolean {false,true} bool;
#endif
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned char uchar;

typedef struct timeval timer;


/* Common definitions */
#define PI          3.1415926536
#define INTERIOR    0
#define EXTERIOR    1
#define BOTH        2
#define WHITE       0
#define GRAY        1
#define BLACK       2
#define NIL        -1
#define INCREASING  1
#define DECREASING  0
#define Epsilon     1E-05

/* Common operations */

  /**
   * \def MAX(x,y)
   * \brief A macro that returns the maximum of \a x and \a y.
   */
#ifndef MAX
#define MAX(x,y) (((x) > (y))?(x):(y))
#endif
  /**
   * \def MIN(x,y)
   * \brief A macro that returns the minimum of \a x and \a y.
   */
#ifndef MIN
#define MIN(x,y) (((x) < (y))?(x):(y))
#endif

#define ROUND(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))

#define SIGN(x) ((x >= 0)?1:-1)

#define SQUARE(x) ((x)*(x))

  /**
   * \brief Vector of four single floats.
   */
  typedef float  v4sf  __attribute__ ((vector_size(16),aligned(16)));
  
  /**
   * \brief Vector of four single integers.
   */
  typedef int    v4si  __attribute__ ((vector_size(16),aligned(16)));

  /**
   * \brief Vector of eight unsigned 8-bit integers.
   */
  typedef uchar  v8qi  __attribute__ ((vector_size(8),aligned(16)));

  /**
   * \brief Vector of sixteen unsigned 8-bit integers.
   */
  typedef uchar  v16qi __attribute__ ((vector_size(16),aligned(16)));

  /**
   * \brief Vector of eight unsigned short integers.
   */
  typedef ushort v8hi  __attribute__ ((vector_size(16),aligned(16)));
  
  
  typedef union _voxel {
    v4si v;
    int  data[4];
    struct{ int x,y,z; } c;
  } Voxel;


  typedef struct _pixel {
    int x,y;
  } Pixel;
  

  char   *AllocCharArray(int n);  /* It allocates 1D array of n characters */


  /**
   * \brief It allocates 1D array of n characters.
   */
  uchar  *AllocUCharArray(int n);

  /**
   * \brief It allocates 1D array of n ushorts.
   */
  ushort *AllocUShortArray(int n);

  uint   *AllocUIntArray(int n); /* It allocates 1D array of n uints */
  
  /**
   * \brief It allocates 1D array of n integers.
   */
  int    *AllocIntArray(int n);

  /**
   * \brief It allocates 1D array of n long long integers.
   */
  long long  *AllocLongLongArray(int n);  
  
  /**
   * \brief It allocates 1D array of n floats.
   */
  float  *AllocFloatArray(int n);

  double *AllocDoubleArray(int n);/* It allocates 1D array of n doubles */

  void    FreeIntArray(int **a);
  void    FreeFloatArray(float **a);
  void    FreeDoubleArray(double **a);
  void    FreeUCharArray(uchar **a);
  void    FreeUShortArray(ushort **a);

  /**
   * \brief It prints error message and exits the program.
   */
  void Error(char *msg,char *func);

  /**
   * \brief It prints warning message and leaves the routine.
   */
  void Warning(char *msg,char *func);

  /**
   * \brief It changes content between a and b.
   */
  inline void SwapInt(int *a, int *b){
    int c;
    c  = *a;  
    *a = *b;  
    *b = c;
  }

  /**
   * \brief It changes content between a and b.
   */
  inline void SwapFloat(float *a, float *b){
    float c;
    c  = *a;
    *a = *b;
    *b = c;
  }

  
  int NCFgets(char *s, int m, FILE *f); /* It skips # comments */
  
  /**
   * Gera um número inteiro aleatório no intervalo [low,high].
   http://www.ime.usp.br/~pf/algoritmos/aulas/random.html
   Execute RandomSeed antes, para atualizar a semente.
  */
  int  RandomInteger (int low, int high);
  void RandomSeed();

} //end gft namespace

#endif

