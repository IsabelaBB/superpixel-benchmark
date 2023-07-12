#ifndef _GFT_IMAGE64_H_
#define _GFT_IMAGE64_H_

#include "gft_common.h"
#include "gft_image32.h"

namespace gft{
  namespace Image64{

    typedef struct _Image64 {
      long long *data;
      long long **array; /* long long signed integer type.
			    At least 64 bits in size.
			    Specified since the C99 version
			    of the standard. */
      int nrows; /* numero de linhas (altura) */
      int ncols; /* numero de colunas (largura) */
      int n;     /* numero de pixels */
    } Image64;


    Image64 *Create(int ncols, int nrows);
    void     Destroy(Image64 **img);
    
    Image64 *ConvertToImage64(Image32::Image32 *img);
    
    Image64 *Read(char *filename);
    void     Write(Image64 *img, char *filename);
    Image64 *ComputeIntegralImage(Image64 *img);
    
    /* recebe uma imagem integral e calcula
       a soma dentro do retangulo com x em [x1, x2] e y em [y1, y2] */
    long long ComputeSum(Image64 *iimg,
			 int x1, int y1,
			 int x2, int y2);
    
    /* Calcula a variancia do retangulo
       com x no intervalo fechado [x1, x2], e y em [y1, y2] */
    double ComputeVariance(Image64 *img,
			   int x1, int y1,
			   int x2, int y2);
    
    
    /* Calcula a variancia do retangulo
       com x no intervalo fechado [x1, x2], e y em [y1, y2] */
    double ComputeVariance(Image64 *iimg,
			   Image64 *iimg2,
			   int x1, int y1,
			   int x2, int y2);
    
    Image64 *Squared(Image64 *img);
    
    
    /* Tamanho do lado do quadrado.*/
    int GetPosWindowMaxSum(Image32::Image32 *img, int tam_x, int tam_y);

  } /*end Image64 namespace*/
} /*end gft namespace*/

#endif

