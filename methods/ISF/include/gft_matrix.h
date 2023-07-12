
#ifndef _GFT_MATRIX_H_
#define _GFT_MATRIX_H_

#include "gft_common.h"
#include "gft_image32.h"

namespace gft{
  namespace Matrix{

    /**
     * It supports both linear and two-dimensional access 
     * (i.e., M->data[p] or M->array[i][j] for an entry
     * (i,j) at address p=j+i*ncols).
     */
    typedef struct _matrix {
      float *data;
      float **array;
      int ncols,nrows;
    } Matrix;


    Matrix *Create(int ncols,int nrows);
    void    Destroy(Matrix **mat);
    Matrix *Clone(Matrix *mat);
    void    Copy(Matrix *dest, 
		 Matrix *src);
    
    Matrix *Invert(Matrix *A);
    Matrix *Transpose(Matrix *A);
    
    Matrix *Mult(Matrix *A, 
		 Matrix *B);
    Matrix *MultByScalar(Matrix *A, float k);

    Matrix *Sub(Matrix *A, 
		Matrix *B);
    Matrix *Add(Matrix *A, 
		Matrix *B);

    float   GetTrace(Matrix *M);
    
    void    Print(Matrix *M);
    void    PrintDimension(Matrix *M);
    
    float   ComputeDistanceL2(Matrix *Y, 
			      Matrix *X);
    
    void    Fill(Matrix *M, float value);
    void    ChangeValue(Matrix *M, 
			float old_value,
			float new_value);
    
    bool    IsValidEntry(Matrix *M,
			 int i, int j);
    
    Image32::Image32 *Convert2Image(Matrix *M);
    
    Matrix *Read(char *filename);
    void    Write(Matrix *M,
		  char *filename);
    
    float   GetMinimumValue(Matrix *M);
    float   GetMaximumValue(Matrix *M);

    Image32::Image32 *Threshold(Matrix *M,
			        float lower, 
			        float higher);
    
    /**
     * @param axis an option (0->x / 1->y / 2->z).
     */
    Matrix* RotationMatrix3(int axis, 
			    float th); 
    
    Matrix* TranslationMatrix3(float dx, float dy, float dz);
    
    Matrix* TransformVoxel(Matrix *m, gft::Voxel v);

  } //end Matrix namespace
} //end gft namespace


#endif


