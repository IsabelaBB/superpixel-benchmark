#ifndef MYRGB2LSB
#define MYRGB2LAB

#include<cmath>

// Change from RGB colour space to LAB colour space



void lookup_table(float * g_s, float * g_r, int max_hw, float sigma_s, float sigma_c) {
    
    //pre-compute spatial gaussian and compute neighborhood half-size
    //FIX-ME: can be more efficient with a pre-computation of the size of the
    //        neighhborhood (less memory allocation)
    
//     int hs     = max_hw;
    float  s2  = sigma_s*sigma_s;
    for(int i=0; i< max_hw; i++){
        float v = exp(-0.5*(i*i)/(s2));
        if(v<0.1){
//             hs = i-1;
            break;
        }
        g_s[i]=v;
    }
    //pre-compute range gaussian
    
    float r2 = sigma_c*sigma_c;
    for(int i=0; i<256; i++){
        g_r[i] = exp(-0.5*(i*i)/r2);
    }
    
}



void RGB2XYZ(unsigned char sR,unsigned char sG,unsigned char sB,float&	X,float& Y,float& Z)
{
    float R = sR/255.0;
    float G = sG/255.0;
    float B = sB/255.0;
    
    float r, g, b;
    
    if(R <= 0.04045)	r = R/12.92;
    else				r = pow((R+0.055)/1.055,2.4);
    if(G <= 0.04045)	g = G/12.92;
    else				g = pow((G+0.055)/1.055,2.4);
    if(B <= 0.04045)	b = B/12.92;
    else				b = pow((B+0.055)/1.055,2.4);
    
    X = r*0.412453 + g*0.357580 + b*0.180423;
    Y = r*0.212671 + g*0.715160 + b*0.072169;
    Z = r*0.019334 + g*0.119193 + b*0.950227;
}

void RGB2LAB(const unsigned char& sR, const unsigned char& sG, const unsigned char& sB, unsigned char& lval, unsigned char& aval, unsigned char& bval)
{
    float X, Y, Z;
    RGB2XYZ(sR, sG, sB, X, Y, Z);
    
    float epsilon = 0.008856;	//actual CIE standard
    float kappa   = 903.3;		//actual CIE standard
    
    float Xr = 0.950456;	//reference white
    float Yr = 1.0;		//reference white
    float Zr = 1.088754;	//reference white
    
    float xr = X/Xr;
    float yr = Y/Yr;
    float zr = Z/Zr;
    
    float fx, fy, fz;
    if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
    else				fx = (kappa*xr + 16.0)/116.0;
    if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
    else				fy = (kappa*yr + 16.0)/116.0;
    if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
    else				fz = (kappa*zr + 16.0)/116.0;
    
    lval = (unsigned char)((116.0*fy-16.0)/100*255+0.5);
    aval = (unsigned char)(500.0*(fx-fy)+128+0.5);
    bval = (unsigned char)(200.0*(fy-fz)+128+0.5);
}



void myrgb2lab1(unsigned char* r,unsigned char* g,unsigned char* b,
        unsigned char* L,unsigned char* A,unsigned char* B,
        int nRows,int nCols
        )
{
    for(int i=0;i<nRows;i++)
        for(int j=0;j<nCols;j++)
            RGB2LAB(r[i*nCols+j],g[i*nCols+j],b[i*nCols+j],L[i*nCols+j],A[i*nCols+j],B[i*nCols+j]);
}



void myrgb2lab2(unsigned char* r,unsigned char* g,unsigned char* b,
        unsigned char* L,unsigned char* A,unsigned char* B,
        int nRows,int nCols, int pw
        )
{
    
    //BILATERAL
    float sigma_s = FLT_MAX;
    float sigma_c = 40;
    float *g_s = (float*) calloc(max(nCols,nRows),sizeof(float));
    float *g_r = (float*) calloc(256,sizeof(float));
    lookup_table(g_s,g_r,max(nCols,nRows),sigma_s,sigma_c);
    
    float * Rf = (float *) calloc(nRows*nCols,sizeof(float));
    float * Gf = (float *) calloc(nRows*nCols,sizeof(float));
    float * Bf = (float *) calloc(nRows*nCols,sizeof(float));
    
    
    for(int i=0;i<nRows;i++) {
        for(int j=0;j<nCols;j++) {
            float count = 0;
            
            int pos_i = i*nCols+j;
            unsigned char valr = r[pos_i];
            unsigned char valg = g[pos_i];
            unsigned char valb = b[pos_i];
            unsigned char val = (valr + valg + valb)/3;
            
            for (int dx=-pw; dx<=pw; dx++) {
                for (int dy=-pw; dy<=pw; dy++) {
                    if ((i+dx<nRows)&&(i+dx>=0)&(j+dy>=0)&(j+dy<nCols)){
                        
                        int pos_d = (i+dx)*nCols+(j+dy);
                        //float kk = kernel[dy+pw+(dx+pw)*(2*pw+1)];
                        unsigned char val2r = r[pos_d];
                        unsigned char val2g = g[pos_d];
                        unsigned char val2b = b[pos_d];
                        unsigned char val2 = (val2r + val2g + val2b)/3;
                        float d_s = g_s[abs(dx)+abs(dy)];
                        float d_rr = g_r[abs(val-val2)];
                        float kk = d_rr*d_s;
                        
                        Rf[pos_i] += (float) (r[pos_d]*kk);
                        Gf[pos_i] += (float) (g[pos_d]*kk);
                        Bf[pos_i] += (float) (b[pos_d]*kk);
                        count += (float) kk;
                    }
                }
            }
            Rf[pos_i]  /= count;
            Gf[pos_i]  /= count;
            Bf[pos_i]  /= count;
            
        }
    }
    
    for(int i=0;i<nRows;i++) {
        for(int j=0;j<nCols;j++) {
            RGB2LAB((unsigned char) Rf[i*nCols+j],(unsigned char) Gf[i*nCols+j],(unsigned char) Bf[i*nCols+j],
                    L[i*nCols+j],A[i*nCols+j],B[i*nCols+j]);
        }
    }
    
    free(Rf);
    free(Gf);
    free(Bf);
    
    
}



#endif