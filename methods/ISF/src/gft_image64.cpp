
#include "gft_image64.h"

namespace gft{
  namespace Image64{

    Image64 *Create(int ncols, int nrows){
      Image64 *img = NULL;
      int i;
      
      img = (Image64 *)malloc(sizeof(Image64));
      if(img == NULL){
	printf("Error: Insufficient memory.\n");
	exit(1);
      }
      img->nrows = nrows;
      img->ncols = ncols;
      img->n = nrows*ncols;
      
      img->data = (long long*)malloc(ncols*nrows*sizeof(long long));
      if(img->data == NULL){
	printf("Error: Insufficient memory.\n");
	exit(1);
      }

      img->array = (long long**)malloc(nrows*sizeof(long long*));
      if(img->array == NULL){
	printf("Error: Insufficient memory.\n");
	exit(1);
      }
      for(i = 0; i < nrows; i++){
	img->array[i] = (img->data + i*ncols);
      }
      return img;
    }
    

    void Destroy(Image64 **img){
      Image64 *tmp;
      if(img!=NULL){
	tmp = *img;
	if(tmp != NULL){
	  if(tmp->data  != NULL) free(tmp->data);
	  if(tmp->array != NULL) free(tmp->array);
	  free(tmp);
	  *img = NULL;
	}
      }
    }

    
    Image64 *ConvertToImage64(Image32::Image32 *img){
      Image64 *img64;
      int p;
      img64 = Create(img->ncols, img->nrows);
      for(p = 0; p < img->ncols*img->nrows; p++){
	img64->data[p] = img->data[p];
      }
      return img64;
    }


    Image64 *Read(char *filename){
      Image64 *img;
      int ncols,nrows,Imax,i,j,val;
      char buf[512];
      FILE *fp;
      
      fp = fopen(filename, "r");
      if(fp == NULL){
	printf("Error: Can't open file.\n");
	exit(1);
      }
      
      fscanf(fp, "%s", buf);
      if(strcmp(buf, "P2") != 0){
	printf("Error: Invalid file type.\n");
	exit(1);
      }
      
      fscanf(fp, "%s", buf);
      if(buf[0] == '#'){
	printf("Error: programa nao le comentarios; remova com um editor de texto primeiro.\n");
	exit(1);
      }
      ncols = atoi(buf);
      
      fscanf(fp, "%d", &nrows);
      fscanf(fp, "%d", &Imax); /* nao e usado */
      
      img = Create(ncols, nrows);
      
      for(i = 0; i < nrows; i++){
	for(j = 0; j < ncols; j++){
	  fscanf(fp, "%d", &val);
	  img->array[i][j] = (long long)val;
	}
      }
      fclose(fp);
      
      return img;
    }


    void Write(Image64 *img, char *filename){
      FILE *fp;
      int i,j;
      fp = fopen(filename, "w");
      if(fp == NULL){
	printf("Error: Can't open file.\n");
	exit(1);
      }
      fprintf(fp, "P2\n");
      fprintf(fp, "%d\n%d\n%d\n", img->ncols, img->nrows, 255);
      for(i = 0; i < img->nrows; i++){
	for(j = 0; j < img->ncols; j++){
	  fprintf(fp, " %3d", (int)img->array[i][j]);
	}
	fprintf(fp,"\n");
      }
      fclose(fp);
    }


    Image64 *ComputeIntegralImage(Image64 *img){
      Image64 *iimg;
      long long sum;
      int i,j;
      
      iimg = Create(img->ncols, img->nrows);
      
      for(i = 0; i < img->nrows; i++){
	for(j = 0; j < img->ncols; j++){
	  sum = img->array[i][j];
	  if(j-1 >= 0)
	    sum += iimg->array[i][j-1];
	  if(i-1 >= 0)
	    sum += iimg->array[i-1][j];
	  if(j-1 >= 0 && i-1 >= 0)
	    sum -= iimg->array[i-1][j-1];
	  
	  iimg->array[i][j] = sum;
	}
      }
      return iimg;
    }
    

    /* recebe uma imagem integral e calcula
       a soma dentro do retangulo com x em [x1, x2] e y em [y1, y2] */
    long long ComputeSum(Image64 *iimg,
			 int x1, int y1,
			 int x2, int y2){
      long long sum;
      
      /*
	if(x1 < 0) x1 = 0;
	if(y1 < 0) y1 = 0;
	if(x2 >= iimg->ncols) x2 = iimg->ncols - 1;
	if(y2 >= iimg->nrows) y2 = iimg->nrows - 1;
      */
      
      sum = iimg->array[y2][x2];
      if(x1 > 0) sum -= iimg->array[y2][x1-1];
      if(y1 > 0) sum -= iimg->array[y1-1][x2];
      if(x1 > 0 && y1 > 0) sum += iimg->array[y1-1][x1-1];
      return sum;
    }


    /* Calcula a variancia do retangulo
       com x no intervalo fechado [x1, x2], e y em [y1, y2] */
    double ComputeVariance(Image64 *img,
			   int x1, int y1,
			   int x2, int y2){
      int i,j,n;
      double var, sum = 0.0, sum2 = 0.0;
      
      if(x1 < 0) x1 = 0;
      if(y1 < 0) y1 = 0;
      if(x2 >= img->ncols) x2 = img->ncols - 1;
      if(y2 >= img->nrows) y2 = img->nrows - 1;
      
      n = (x2-x1+1)*(y2-y1+1);
      for(i = y1; i <= y2; i++){
	for(j = x1; j <= x2; j++){
	  sum  += img->array[i][j];
	  sum2 += (img->array[i][j] * img->array[i][j]);
	}
      }
      var = sum2/n - (sum/n)*(sum/n);
      return var;
    }
    
    
    /* Calcula a variancia do retangulo
       com x no intervalo fechado [x1, x2], e y em [y1, y2] */
    double ComputeVariance(Image64 *iimg,
			   Image64 *iimg2,
			   int x1, int y1,
			   int x2, int y2){
      double sum, sum2;
      double var;
      int n;
      n = (x2-x1+1)*(y2-y1+1);
      sum  = (double)ComputeSum(iimg,  x1, y1, x2, y2);
      sum2 = (double)ComputeSum(iimg2, x1, y1, x2, y2);
      var = sum2/n - (sum/n)*(sum/n);
      return var;
    }
    
    Image64 *Squared(Image64 *img){
      Image64 *img2;
      int i,j;
      img2 = Create(img->ncols, img->nrows);
      for(i = 0; i < img->nrows; i++){
	for(j = 0; j < img->ncols; j++){
	  img2->array[i][j] = img->array[i][j]*img->array[i][j];
	}
      }
      return img2;
    }
    
    
    
    /* Tamanho do lado do quadrado.*/
    int GetPosWindowMaxSum(Image32::Image32 *img, int tam_x, int tam_y){
      Image64 *img64,*iimg;
      long long max,val;
      int i,j,pmax = 0;
      
      img64 = ConvertToImage64(img);
      iimg  = ComputeIntegralImage(img64);
      max = INT_MIN;
      for(i = 0; i+tam_y-1 < img->nrows; i++){
	for(j = 0; j+tam_x-1 < img->ncols; j++){
	  val = ComputeSum(iimg, j, i, j+tam_x-1, i+tam_y-1);
	  if(val > max){
	    max = val;
	    pmax = j + i*img->ncols;
	  }
	}
      }
      Destroy(&img64);
      Destroy(&iimg);
      return pmax;
    }



  } /*end Image64 namespace*/
} /*end gft namespace*/

