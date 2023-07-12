
#include "gft_color.h"

namespace gft{
  namespace Color{

    int Triplet(int a,int b,int c) {
      return(((a&0xff)<<16)|((b&0xff)<<8)|(c&0xff));
    }
    
    
    #define GMAX(a,b,c) (((a>=b)&&(a>=c))?a:((b>c)?b:c))
    #define GMIN(a,b,c) (((a<=b)&&(a<=c))?a:((b<c)?b:c))

    double f(const double &t){
      if(t > 0.008856)
	return cbrt(t); //pow(t, 1.0/3.0);
      else
	return (903.3*t + 16.0)/116.0;
    }

    float f(const float &t){
      if(t > 0.008856)
	return cbrtf(t); //pow(t, 1.0/3.0);
      else
	return (903.3*t + 16.0)/116.0;
    }

    
    void RGB2Lab(const int& R, const int& G, const int& B,
		 double &l, double &a, double &b){
      const double tresh = 0.008856;
      double X,Y,Z;
      const double Xn = 0.950456;  //reference white
      const double Yn = 1.0;       //reference white
      const double Zn = 1.088754;  //reference white

      double R0 = (R)/255.0;
      double G0 = (G)/255.0;
      double B0 = (B)/255.0;

      double fx,fy,fz;
      
      if(R0 <= 0.04045)	R0 = R0/12.92;
      else	        R0 = pow((R0+0.055)/1.055,2.4);

      if(G0 <= 0.04045)	G0 = G0/12.92;
      else		G0 = pow((G0+0.055)/1.055,2.4);

      if(B0 <= 0.04045)	B0 = B0/12.92;
      else		B0 = pow((B0+0.055)/1.055,2.4);

      /*
      // http://www.cs.rit.edu/~ncs/color/t_convert.html
      X =  0.412453 * R0 + 0.357580 * G0 + 0.180423 * B0;
      Y =  0.212671 * R0 + 0.715160 * G0 + 0.072169 * B0;
      Z =  0.019334 * R0 + 0.119193 * G0 + 0.950227 * B0;
      */

      X = R0 * 0.4124564 + G0 * 0.3575761 + B0 * 0.1804375;
      Y = R0 * 0.2126729 + G0 * 0.7151522 + B0 * 0.0721750;
      Z = R0 * 0.0193339 + G0 * 0.1191920 + B0 * 0.9503041;

      //printf("XYZ: (%f, %f, %f)\n",X,Y,Z);

      //-------------------------------------
      /*
      fx = f(X/Xn);
      fy = f(Y/Yn);
      fz = f(Z/Zn);
      */
      //-------------------------------------
      X = X/Xn;
      Y = Y/Yn;
      Z = Z/Zn;
      
      if(X > tresh)
	fx = cbrt(X); //pow(X, 1.0/3.0);
      else
	fx = (903.3*X + 16.0)/116.0;

      if(Y > tresh)
	fy = cbrt(Y); //pow(Y, 1.0/3.0);
      else
	fy = (903.3*Y + 16.0)/116.0;

      if(Z > tresh)
	fz = cbrt(Z); //pow(Z, 1.0/3.0);
      else
	fz = (903.3*Z + 16.0)/116.0;
      
      //-------------------------------------
      
      /*
      if(Y/Yn > tresh)
	lab.l = 116.0 * cbrt(Y/Yn) - 16.0;
      else
	lab.l = 903.3 * Y/Yn;
      */
      l = 116.0*fy - 16.0;
      a = 500.0 * (fx - fy);
      b = 200.0 * (fy - fz);
    }
    

    void RGB2Lab(const int& R, const int& G, const int& B,
		 float &l, float &a, float &b){
      float X,Y,Z;
      const float Xn = 0.950456;  //reference white
      const float Yn = 1.0;       //reference white
      const float Zn = 1.088754;  //reference white

      float R0 = (R)/255.0;
      float G0 = (G)/255.0;
      float B0 = (B)/255.0;

      float fx,fy,fz;
      
      if(R0 <= 0.04045)	R0 = R0/12.92;
      else	        R0 = powf((R0+0.055)/1.055,2.4);

      if(G0 <= 0.04045)	G0 = G0/12.92;
      else		G0 = powf((G0+0.055)/1.055,2.4);

      if(B0 <= 0.04045)	B0 = B0/12.92;
      else		B0 = powf((B0+0.055)/1.055,2.4);

      /*
      // http://www.cs.rit.edu/~ncs/color/t_convert.html
      X =  0.412453 * R0 + 0.357580 * G0 + 0.180423 * B0;
      Y =  0.212671 * R0 + 0.715160 * G0 + 0.072169 * B0;
      Z =  0.019334 * R0 + 0.119193 * G0 + 0.950227 * B0;
      */

      X = R0 * 0.4124564 + G0 * 0.3575761 + B0 * 0.1804375;
      Y = R0 * 0.2126729 + G0 * 0.7151522 + B0 * 0.0721750;
      Z = R0 * 0.0193339 + G0 * 0.1191920 + B0 * 0.9503041;

      //printf("XYZ: (%f, %f, %f)\n",X,Y,Z);

      fx = f(X/Xn);
      fy = f(Y/Yn);
      fz = f(Z/Zn);

      /*
      if(Y/Yn > tresh)
	lab.l = 116.0 * cbrt(Y/Yn) - 16.0;
      else
	lab.l = 903.3 * Y/Yn;
      */
      l = 116.0*fy - 16.0;
      a = 500.0 * (fx - fy);
      b = 200.0 * (fy - fz);
    }
    

    /*
    ColorLab RGB2Lab(ColorRGB rgb){
      ColorLab lab;
      float r,g,b;
      float tresh = 0.008856;
      //  Bergo:    
      float Xn  = 0.950456;
      float Yn  = 1.0;
      float Zn  = 1.088854;
      
      //http://www.color-image.com/2011/10/the-reference-white-in-adobe-photoshop-lab-mode/
      //float Xn = 0.9642;
      //float Yn = 1.0;
      //float Zn = 0.8249;

      float r13    = 1.0/3.0;
      float r16116 = 16.0/116.0;
      float cnst   = 7.787;

      r = rgb.r;
      g = rgb.g;
      b = rgb.b;      

      float R0 = (r)/255.0;
      float G0 = (g)/255.0;
      float B0 = (b)/255.0;

      //R0 = powf(R0,2.2); /// gamma
      //G0 = powf(G0,2.2); /// gamma
      //B0 = powf(B0,2.2); /// gamma
      
      // sRGB (d65) Matrix Cxr
      float X =  0.412453 * R0 + 0.357580 * G0 + 0.180423 * B0;
      float Y =  0.212671 * R0 + 0.715160 * G0 + 0.072169 * B0;
      float Z =  0.019334 * R0 + 0.119193 * G0 + 0.950227 * B0;
      
      float xr = X/Xn;
      float yr = Y/Yn;
      float zr = Z/Zn;
      
      float L, A, B;
      float xr0, yr0, zr0;


      if (xr > tresh){
        xr0 = pow(xr, r13);
      }
      else{
        xr0 = cnst * xr + r16116;
      }

      if (yr > tresh){
        yr0 = pow(yr, r13);
      }
      else{
        yr0 = cnst * yr + r16116;
      }

      if (zr > tresh){
        zr0 = pow(zr, r13);
      }
      else{
        zr0 = cnst * zr + r16116;
      }

      if (yr > tresh){
        L = 116.0 * pow(yr, r13) - 16.0;
      }
      else{
        L = 903.3 * yr;
      }

      A = 500.0 * (xr0 - yr0);
      B = 200.0 * (yr0 - zr0);

      lab.l = L;
      lab.a = A;
      lab.b = B;
      return lab;
    }
    */

    /*
    ColorRGB Lab2RGB(ColorLab lab){
      ColorRGB rgb;
      float thresh1 = 0.008856;
      float thresh2 = 0.206893;
      
      float Xn = 0.950456;
      float Yn = 1.000;
      float Zn = 1.088854;

      float r13    = 1.0/3.0;
      float r16116 = 16.0/116.0;
      float cnst   = 7.787;

      float X, Y, Z;
      float L = lab.l;
      float A = lab.a;
      float B = lab.b;

      float P = (L+16.0)/116.0;
      float yr,fy;

      float fx;
      float fz;
      float R1,G1,B1;

      if (L > 7.9996){
	Y = Yn * P * P * P;
      }
      else{
        Y = L / 903.3;
      }

      yr = Y/Yn;
      if ( yr > thresh1 ){
        fy = pow(yr, r13);
      }
      else{
        fy = cnst * yr + r16116;
      }

      fx = A / 500.0 + fy;
      fz = fy - B / 200.0;

      if (fx > thresh2 ){
        X = Xn * fx * fx * fx;
      }
      else{
        X = Xn/cnst * ( fx - r16116 );
      }

      if (fz > thresh2 ){
        Z = Zn * fz * fz * fz;
      }
      else{
        Z = Zn/cnst * ( fz - r16116 );
      }

      // sRGB (d65) Matrix Crx
      R1 =   3.240479 * X - 1.537150 * Y - 0.498535 * Z;
      G1 = - 0.969256 * X + 1.875992 * Y + 0.041556 * Z;
      B1 =   0.055648 * X - 0.204043 * Y + 1.057311 * Z;

      R1 = powf(R1,1/2.2); /// gamma
      G1 = powf(G1,1/2.2); /// gamma
      B1 = powf(B1,1/2.2); /// gamma

      rgb.r = ROUND(MAX(0.0, MIN(255.0, 255.0*R1)));
      rgb.g = ROUND(MAX(0.0, MIN(255.0, 255.0*G1)));
      rgb.b = ROUND(MAX(0.0, MIN(255.0, 255.0*B1)));
      return rgb;
    }
    */

    ColorHSV RGB2HSV(ColorRGB rgb){
      ColorHSV hsv;
      float r,g,b,max,min;
      float h,s,v;
      
      r = (float)rgb.r/255.0;
      g = (float)rgb.g/255.0;
      b = (float)rgb.b/255.0;
      
      max = MAX(r,MAX(g,b));
      min = MIN(r,MIN(g,b));
      
      if(max==min)
	h = 0;
      else if(max == r && g >= b)
	h = 60.0*(g-b)/(max-min) + 0.0;
      else if(max == r && g < b)
	h = 60.0*(g-b)/(max-min) + 360.0;
      else if(max == g)
	h = 60.0*(b-r)/(max-min) + 120.0;
      else /* max == b*/
	h = 60.0*(r-g)/(max-min) + 240.0;
      
      s = (max-min)/max;
      v = max;
      
      hsv.h = ROUND(h);         // [0..360]
      hsv.s = ROUND(s * 255.0); // [0..255]
      hsv.v = ROUND(v * 255.0); // [0..255]
      return hsv;
    }


    ColorRGB HSV2RGB(ColorHSV hsv){
      ColorRGB rgb;
      float s,v;
      float r,g,b;
      float p,q,t;
      int Hi;
      float f;
      
      s = (float)hsv.s/255.0;
      v = (float)hsv.v/255.0;
      
      if(hsv.s == 0){ /* resultado cinza */
	r = v;
	g = v;
	b = v;
      }
      else{
	Hi = (hsv.h/60) % 6;
	f = hsv.h/60.0 - Hi;
	p = v*(1.0-s);
	q = v*(1.0-f*s);
	t = v*(1.0-(1.0-f)*s);
	switch(Hi){
	case 0: r = v; g = t; b = p; break;
	case 1: r = q; g = v; b = p; break;
	case 2: r = p; g = v; b = t; break;
	case 3: r = p; g = q; b = v; break;
	case 4: r = t; g = p; b = v; break;
	case 5: r = v; g = p; b = q; break;
	}
      }
      rgb.r = ROUND(r*255.0);
      rgb.g = ROUND(g*255.0);
      rgb.b = ROUND(b*255.0);
      return rgb;
    }


    /* returns a random color */
    int RandomColor(){
      static int seeded = 0;
      int a,b,c;
      
      if (!seeded) {
	srand(time(0));
	seeded = 1;
      }
      
      a = 1+(int) (255.0*rand()/(RAND_MAX+1.0));
      b = 1+(int) (255.0*rand()/(RAND_MAX+1.0));
      c = 1+(int) (255.0*rand()/(RAND_MAX+1.0));
      return(Triplet(a,b,c));
    }

    
    /*
      MergeRGB
      Merges 2 colors, (1-ratio) of a with (ratio) of b.
      Works in RGB space.
      
      input:  2 RGB triplets, merge ratio
      output: RGB triplet
      author: bergo
    */
    int MergeRGB(int a,int b,float ratio){
      float c[6];
      int r[3];
      
      c[0] = (float) (a >> 16);
      c[1] = (float) (0xff & (a >> 8));
      c[2] = (float) (0xff & a);
      
      c[3] = (float) (b >> 16);
      c[4] = (float) (0xff & (b >> 8));
      c[5] = (float) (0xff & b);
      
      c[0] = (1.0-ratio) * c[0] + ratio * c[3];
      c[1] = (1.0-ratio) * c[1] + ratio * c[4];
      c[2] = (1.0-ratio) * c[2] + ratio * c[5];
      
      r[0] = (int) c[0];
      r[1] = (int) c[1];
      r[2] = (int) c[2];
      
      return( (r[0]<<16) | (r[1]<<8) | r[2] );
    }
    
    
    /*    
    int RGB2HSV(int vi) {
      float r = (float) Channel0(vi),
	    g = (float) Channel1(vi),
    	    b = (float) Channel2(vi), v, x, f;
      float a[3];
      
      int i;
      
      r/=255.0;
      g/=255.0;
      b/=255.0;
      
      // RGB are each on [0, 1]. S and V are returned on [0, 1] and H is returned on [0, 6]. 
      
      x = GMIN(r, g, b);
      v = GMAX(r, g, b);
      if (v == x) {
	a[0]=0.0;
	a[1]=0.0;
	a[2]=v;
      } else {
	f = (r == x) ? g - b : ((g == x) ? b - r : r - g);
	i = (r == x) ? 3 : ((g == x) ? 5 : 1);
	a[0]=((float)i)-f/(v-x);
	a[1]=(v-x)/v;
	//a[2]=v;
	a[2]=0.299*r+0.587*g+0.114*b;
      }
      
      // (un)normalize
      a[0]*=60.0;
      a[1]*=255.0;
      a[2]*=255.0;
      
      return(Triplet(ROUND(a[0]),ROUND(a[1]),ROUND(a[2])));
    }
    */

    /*    
    int HSV2RGB(int vi) {
      // H is given on [0, 6]. S and V are given on [0, 1]. 
      // RGB are each returned on [0, 1].
      float h = (float)Channel0(vi),
            s = (float)Channel1(vi),
            v = (float)Channel2(vi), m, n, f;
      float a[3]={0,0,0};
      int i;
      
      h/=60.0;
      s/=255.0;
      v/=255.0;
      
      if (s==0.0) {
	a[0]=a[1]=a[2]=v;
      } else {
	i = (int) floor(h);
	f = h - (float)i;
	if(!(i & 1)) f = 1 - f; // if i is even
	m = v * (1 - s);
	n = v * (1 - s * f);
	switch (i) {
	case 6:
	case 0: a[0]=v; a[1]=n; a[2]=m; break;
	case 1: a[0]=n; a[1]=v; a[2]=m; break;
	case 2: a[0]=m; a[1]=v; a[2]=n; break;
	case 3: a[0]=m; a[1]=n; a[2]=v; break;
	case 4: a[0]=n; a[1]=m; a[2]=v; break;
	case 5: a[0]=v; a[1]=m; a[2]=n; break;
	}
      }
      
      // (un)normalize
      for(i=0;i<3;i++)
	a[i]*=255;
      
      return(Triplet(ROUND(a[0]),ROUND(a[1]),ROUND(a[2])));
    }
    */


  } /*end Color namespace*/
} /*end gft namespace*/

