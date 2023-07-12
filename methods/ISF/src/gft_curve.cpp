
#include "gft_curve.h"

namespace gft{
  namespace Curve{

    Curve *Create(int n){
      Curve *curve=NULL;

      curve = (Curve *) calloc(1,sizeof(Curve));
      if(curve != NULL){
	curve->X = gft::AllocFloatArray(n);
	curve->Y = gft::AllocFloatArray(n);
	curve->n = n;
      } 
      else gft::Error((char *)MSG1,
		      (char *)"Curve::Create");
      return(curve);
    }


    void   Destroy(Curve **curve){
      Curve *aux;

      aux = *curve;
      if(aux != NULL){
	if(aux->X != NULL) gft::FreeFloatArray(&aux->X);
	if(aux->Y != NULL) gft::FreeFloatArray(&aux->Y);
	free(aux);
	*curve = NULL;
      }
    }

    Curve *Clone(Curve *curve){
      Curve *curvec;
      
      curvec = Create(curve->n);
      memcpy(curvec->X,curve->X,curve->n*sizeof(float));
      memcpy(curvec->Y,curve->Y,curve->n*sizeof(float));
      return(curvec);
    }


    //Compatible with Gnuplot.
    Curve *Read(char *filename){
      FILE *fp;
      int i,n;
      Curve *curve;
      float x,y;
      
      fp = fopen(filename,"r");
      if(fp == NULL)
	gft::Error((char *)"Cannot open file",
		   (char *)"Curve::Read");

      fscanf(fp,"#%d\n",&n);
      curve= Create(n);
      for(i=0; i<curve->n; i++) {
	fscanf(fp,"%f %f\n",&x,&y);
	curve->X[i] = x;
	curve->Y[i] = y;
      }
      fclose(fp);
      return (curve);
    }


    void   Write(Curve *curve,char *filename){
      FILE *fp;
      int i;
      
      fp = fopen(filename,"w");
      if(fp == NULL)
	gft::Error((char *)"Cannot open file",
		   (char *)"Curve::Write");

      fprintf(fp,"#%d\n",curve->n);
      for(i=0; i<curve->n; i++)
	fprintf(fp,"%f %f\n",curve->X[i],curve->Y[i]);
      
      fclose(fp);
    }


    int    LowerPercentage(Curve *curve,
			   float perc){
      Curve *ncurve;
      float sum=0.0,ratio;
      int i=0;
      
      ratio = perc/100.0;
      ncurve = Normalize(curve);
      while(sum<ratio && i<ncurve->n){
	sum += ncurve->Y[i];
	i++;
      }
      i = MAX(0, i-1);
      Destroy(&ncurve);
      return i;
    }

    
    int    HigherPercentage(Curve *curve,
			    float perc){
      Curve *ncurve;
      float sum=0.0,ratio;
      int i;
      
      ratio = 1.0 - perc/100.0;
      ncurve = Normalize(curve);
      i = ncurve->n-1;
      while(sum<ratio && i>=0){
	sum += ncurve->Y[i];
	i--;
      }
      i = MIN(ncurve->n-1, i+1);
      Destroy(&ncurve);
      return i;
    }


    int    Median(Curve *curve){
      float count,n;
      int i;
      n = 0.0;
      for(i=0; i<curve->n; i++)
	n += curve->Y[i];

      i = -1;
      count = 0.0;
      while(count<n/2.0){
	i++;
	count += curve->Y[i];
      }
      return i;
    }

    
    int    Otsu(Curve *curve){
      Curve *ncurve=Normalize(curve);
      double p1,p2,m1,m2,s1,s2,J,Jmax=-1.0;
      int i,T,Topt=0,Imax=ncurve->n-1;
      
      for(T=1; T<Imax; T++){
	p1 = 0.0;
	for(i=0; i<=T; i++) 
	  p1 += ncurve->Y[i];
	p2 = 1.0 - p1;
	if((p1 > 0.0)&&(p2 > 0.0)){
	  m1 = 0.0;
	  for(i=0; i<=T; i++) 
	    m1 += ncurve->Y[i]*i;
	  m1 /= p1;
	  m2 = 0.0;
	  for(i=T+1; i<=Imax; i++) 
	    m2 += ncurve->Y[i]*i;
	  m2 /= p2;
	  s1 = 0.0;
	  for(i=0; i<=T; i++) 
	    s1 += ncurve->Y[i]*(i-m1)*(i-m1);
	  s1 /= p1;
	  s2 = 0.0;
	  for (i=T+1; i<=Imax; i++) 
	    s2 += ncurve->Y[i]*(i-m2)*(i-m2);
	  s2 /= p2;
	  J = (p1*p2*(m1-m2)*(m1-m2))/(p1*s1+p2*s2);
	}else{
	  J = 0.0;      
	}
	if (J > Jmax){
	  Jmax = J;
	  Topt = T;
	}
      }
      Destroy(&ncurve);
      return(Topt);
    }


    Curve *Normalize(Curve *curve){
      Curve *ncurve;
      ncurve = Clone(curve);
      Normalizeinplace(ncurve);
      return (ncurve);
    }


    void Normalizeinplace(Curve *curve){
      float sum;
      int i;
      sum = 0.0;
      for(i=0; i<curve->n; i++)
	sum += curve->Y[i];
      for(i=0; i<curve->n; i++)
	curve->Y[i] /= sum;
    }

    void Sort(Curve *curve){
      int i,j,m;
      float tmp;      
      for(i = curve->n-1; i > 0; i--){
	/* Coloca em dados[i] o maximo de [0,i]. */
	m = 0;
	for(j = 1; j <= i; j++){
	  if(curve->Y[j] > curve->Y[m])
	    m = j;
	}
	tmp = curve->Y[i];
	curve->Y[i] = curve->Y[m];
	curve->Y[m] = tmp;

	tmp = curve->X[i];
	curve->X[i] = curve->X[m];
	curve->X[m] = tmp;
      }
    }


    void InvertXY(Curve *curve){
      double tmp;
      int i;
      for (i=0; i<curve->n; i++){
	tmp = curve->X[i];
	curve->X[i] = curve->Y[i];
	curve->Y[i] = tmp;
      }
    }


  } //end Curve namespace
} //end gft namespace

