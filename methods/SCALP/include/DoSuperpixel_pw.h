#ifndef DOSUPERPIXEL
#define DOSUPERPIXEL

#include<vector>
#include"preEnforceConnectivity.h"
#include<algorithm>
#include"EnforceConnectivity.h"
#include"point.h"
#include<string.h>

using namespace std;


//Perform weighted kmeans iteratively in the ten dimensional feature space.

void lookup_table2(float * g_s, float * g_r, int max_hw, float sigma_s, float sigma_c) {
    
    //int hs;      = max_hw;
    float  s2  = sigma_s*sigma_s;
    for(int i=0; i< max_hw; i++){
        float v = exp(-0.5*(i*i)/(s2));
        if(v<0.1){
            //hs = i-1;
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


float dist_bresenham_1(int xc,int yc, int xx, int yy, int height,
        float * L1_1, float *L2_1, float * a1_1, float *a2_1, float *b1_1, float *b2_1,
        //float * L1_2, float *L2_2, float * a1_2, float *a2_2, float *b1_2, float *b2_2,
        float centerL1, float centerL2, float centera1, float centera2, float centerb1, float centerb2,
        //float centerx1, float centerx2, float centery1, float centery2,
        float Lab_cc, float xy_cc, float * Lab_2, float *xy_2, float* arr)
{
    int dx, dy, i, e;
    int incx, incy, inc1, inc2;
    int x,y;
    float count = 0.0;
    float max_count = 20;
    
    dx = xc - xx;
    dy = yc - yy;
    
    if(dx < 0)
        dx = -dx;
    if(dy < 0)
        dy = -dy;
    
    incx = 1;
    if(xc < xx)
        incx = -1;
    
    incy = 1;
    if(yc < yy)
        incy = -1;
    
    x=xx;
    y=yy;
    
    float dist = 0;
    int ind;
    float tmp = 0;
    
    if(dx > dy)
    {
        e = 2*dy - dx;
        inc1 = 2*( dy -dx);
        inc2 = 2*dy;
        for(i = 0; i < dx-1; i++)
        {
            if(e >= 0)
            {
                y += incy;
                e += inc1;
            }
            else
                e += inc2;
            
            x += incx;
            
            ind = x*height + y;
            
            
            if (count > max_count)
                return dist/count;
            else {
                
                if (arr[ind] < 0) {
                    tmp = Lab_cc + Lab_2[ind] - 2*(centerL1*L1_1[ind] + centerL2*L2_1[ind] + centera1*a1_1[ind] +
                            centera2*a2_1[ind] + centerb1*b1_1[ind] + centerb2*b2_1[ind]);
                    arr[ind] = tmp;
                }
                dist += arr[ind];
                
                count += 1;
            }
        }
    }
    else
    {
        e = 2*dx - dy;
        inc1 = 2*( dx - dy);
        inc2 = 2*dx;
        
        for(i = 0; i < dy-1; i++)
        {
            if(e >= 0)
            {
                x += incx;
                e += inc1;
            }
            else
                e += inc2;
            
            y += incy;
            
            ind = x*height + y;
            if (count > max_count)
                return dist/count;
            else {
                
                if (arr[ind] < 0) {
                    tmp = Lab_cc + Lab_2[ind] - 2*(centerL1*L1_1[ind] + centerL2*L2_1[ind] + centera1*a1_1[ind] +
                            centera2*a2_1[ind] + centerb1*b1_1[ind] + centerb2*b2_1[ind]);
                    arr[ind] = tmp;
                }
                dist += arr[ind];
                
                count += 1;
                
            }
        }
    }
    
    if (count > 0)
        return dist/count;
    else
        return dist;
    
}



float dist_bresenham_contour(int x1,int y1,int xc,int yc, float * contour, int height)
{
    int dx, dy, i, e;
    int incx, incy, inc1, inc2;
    int x,y;
    float count = 0.0;
    int max_count = 20;
    float max_contour = 0;
    
    dx = xc - x1;
    dy = yc - y1;
    
    if(dx < 0)
        dx = -dx;
    if(dy < 0)
        dy = -dy;
    
    incx = 1;
    if(xc < x1)
        incx = -1;
    
    incy = 1;
    if(yc < y1)
        incy = -1;
    
    x=x1;
    y=y1;
    
    int ind;
    
    
    if(dx > dy)
    {
        e = 2*dy - dx;
        inc1 = 2*( dy -dx);
        inc2 = 2*dy;
        for(i = 0; i < dx-1; i++)
        {
            if(e >= 0)
            {
                y += incy;
                e += inc1;
            }
            else
                e += inc2;
            
            x += incx;
            ind = x*height + y;
            
            if (count > max_count)
                return max_contour;
            else {
                if (contour[ind] > max_contour)
                    max_contour = contour[ind];
                count += 1;
            }
            
            
            
        }
    }
    else
    {
        e = 2*dx - dy;
        inc1 = 2*( dx - dy);
        inc2 = 2*dx;
        
        for(i = 0; i < dy-1; i++)
        {
            if(e >= 0)
            {
                x += incx;
                e += inc1;
            }
            else
                e += inc2;
            
            y += incy;
            
            ind = x*height + y;
            
            if (count > max_count)
                return max_contour;
            else {
                
                if (contour[ind] > max_contour)
                    max_contour = contour[ind];
                count += 1;
            }
        }
    }
    
    return max_contour;
}



void DoSuperpixel(
        float** L1,
        float** L2,
        float** a1,
        float** a2,
        float** b1,
        float** b2,
        float** x1,
        float** x2,
        float** y1,
        float** y2,
        float** W,
        unsigned short int* label,
        point* seedArray,
        int seedNum,
        int nCols,
        int nRows,
        int StepX,
        int StepY,
        int iterationNum,
        int thresholdCoef,
        float * kernel,
        int bres_color,
        int bres_contour,
        float * contour,
        point* seedArray_p,
        int pw,
        unsigned char* R,
        unsigned char* G,
        unsigned char* B
        )
{
    
    float* dist=new float[nCols*nRows]();
    
    float* centerL1=new float[seedNum];
    float* centerL2=new float[seedNum];
    float* centera1=new float[seedNum];
    float* centera2=new float[seedNum];
    float* centerb1=new float[seedNum];
    float* centerb2=new float[seedNum];
    float* centerx1=new float[seedNum];
    float* centerx2=new float[seedNum];
    float* centery1=new float[seedNum];
    float* centery2=new float[seedNum];
    float* WSum=new float[seedNum];
    int* clusterSize=new int[seedNum];
    
    
    float* L1_1=new float[nCols*nRows]();
    float* L1_2=new float[nCols*nRows]();
    float* L2_1=new float[nCols*nRows]();
    float* L2_2=new float[nCols*nRows]();
    float* a1_1=new float[nCols*nRows]();
    float* a1_2=new float[nCols*nRows]();
    float* b1_1=new float[nCols*nRows]();
    float* b1_2=new float[nCols*nRows]();
    float* a2_1=new float[nCols*nRows]();
    float* a2_2=new float[nCols*nRows]();
    float* b2_1=new float[nCols*nRows]();
    float* b2_2=new float[nCols*nRows]();
    float* x1_1=new float[nCols*nRows]();
    float* x1_2=new float[nCols*nRows]();
    float* y1_1=new float[nCols*nRows]();
    float* y1_2=new float[nCols*nRows]();
    float* x2_1=new float[nCols*nRows]();
    float* x2_2=new float[nCols*nRows]();
    float* y2_1=new float[nCols*nRows]();
    float* y2_2=new float[nCols*nRows]();
    
    float* Lab_2=new float[nCols*nRows]();
    float* xy_2=new float[nCols*nRows]();
    
    
    //BILATERAL
    float sigma_s = FLT_MAX;
    float sigma_c = 40;
    float *g_s = (float*) calloc(max(nRows,nCols),sizeof(float));
    float *g_r = (float*) calloc(256,sizeof(float));
    lookup_table2(g_s,g_r,max(nRows,nCols),sigma_s,sigma_c);
    
    
    
    for(int i=0;i<nCols;i++) {
        for(int j=0;j<nRows;j++) {
            float count = 0;
            
            int pos_i = i*nRows+j;
            unsigned char valr = R[pos_i];
            unsigned char valg = G[pos_i];
            unsigned char valb = B[pos_i];
            unsigned char val = (valr + valg + valb)/3;
            
            x1_1[pos_i] = (float) (x1[i][j]);
            y1_1[pos_i] = (float) y1[i][j];
            x2_1[pos_i] = (float) (x2[i][j]);
            y2_1[pos_i] = (float) (y2[i][j]);
            
            xy_2[pos_i] = (float) x1[i][j]*x1[i][j] +  x2[i][j]*x2[i][j] +  y1[i][j]*y1[i][j] +  y2[i][j]*y2[i][j];
            
            
            for (int dx=-pw; dx<=pw; dx++) {
                for (int dy=-pw; dy<=pw; dy++) {
                    if ((i+dx<nCols)&&(i+dx>=0)&(j+dy>=0)&(j+dy<nRows)){
                        
                        int pos_d = (i+dx)*nRows+(j+dy);
                        unsigned char val2r = R[pos_d];
                        unsigned char val2g = G[pos_d];
                        unsigned char val2b = B[pos_d];
                        unsigned char val2 = (val2r + val2g + val2b)/3;
                        float d_s = g_s[abs(dx)+abs(dy)];
                        float d_rr = g_r[abs(val-val2)];
                        float kk = d_rr*d_s;
                        L1_1[pos_i] += (float) (L1[i+dx][j+dy]*kk);
                        L2_1[pos_i] += (float) (L2[i+dx][j+dy]*kk);
                        a1_1[pos_i] += (float) (a1[i+dx][j+dy]*kk);
                        b1_1[pos_i] += (float) (b1[i+dx][j+dy]*kk);
                        a2_1[pos_i] += (float) (a2[i+dx][j+dy]*kk);
                        b2_1[pos_i] += (float) (b2[i+dx][j+dy]*kk);
                        Lab_2[pos_i] += (float) (L1[i+dx][j+dy]*L1[i+dx][j+dy] +
                                L2[i+dx][j+dy]*L2[i+dx][j+dy] +
                                a1[i+dx][j+dy]*a1[i+dx][j+dy] +
                                a2[i+dx][j+dy]*a2[i+dx][j+dy] +
                                b1[i+dx][j+dy]*b1[i+dx][j+dy] +
                                b2[i+dx][j+dy]*b2[i+dx][j+dy])*kk;
                        
                        
                        count += (float) kk;
                    }
                }
            }
            L1_1[pos_i]  /= count;
            L1_2[pos_i]  /= count;
            L2_1[pos_i]  /= count;
            a1_1[pos_i]  /= count;
            b1_1[pos_i]  /= count;
            a2_1[pos_i]  /= count;
            b2_1[pos_i]  /= count;
            
            Lab_2[pos_i] /= count;
            
        }
    }

    
//Initialization
    for(int i=0;i<seedNum;i++)
    {
        centerL1[i]=0;
        centerL2[i]=0;
        centera1[i]=0;
        centera2[i]=0;
        centerb1[i]=0;
        centerb2[i]=0;
        centerx1[i]=0;
        centerx2[i]=0;
        centery1[i]=0;
        centery2[i]=0;
        int x=seedArray[i].x;
        int y=seedArray[i].y;
        int minX=(x-StepX/4<=0)?0:x-StepX/4;
        int minY=(y-StepY/4<=0)?0:y-StepY/4;
        int maxX=(x+StepX/4>=nCols-1)?nCols-1:x+StepX/4;
        int maxY=(y+StepY/4>=nRows-1)?nRows-1:y+StepY/4;
        int Count=0;
        for(int j=minX;j<=maxX;j++)
            for(int k=minY;k<=maxY;k++)
            {
                Count++;
                
                centerL1[i]+=L1[j][k];
                centerL2[i]+=L2[j][k];
                centera1[i]+=a1[j][k];
                centera2[i]+=a2[j][k];
                centerb1[i]+=b1[j][k];
                centerb2[i]+=b2[j][k];
                centerx1[i]+=x1[j][k];
                centerx2[i]+=x2[j][k];
                centery1[i]+=y1[j][k];
                centery2[i]+=y2[j][k];
                
            }
        centerL1[i]/=Count;
        centerL2[i]/=Count;
        centera1[i]/=Count;
        centera2[i]/=Count;
        centerb1[i]/=Count;
        centerb2[i]/=Count;
        centerx1[i]/=Count;
        centerx2[i]/=Count;
        centery1[i]/=Count;
        centery2[i]/=Count;
    }
    
    
    
    float * min_dist_center = (float *)calloc(seedNum,sizeof(float));
  
    float * arr = (float *)malloc(nCols*nRows*sizeof(float));
    
    
    float lambda = 1;
    if (bres_color == 1)
        lambda = 0.5;
    
    //K-means
    for(int iteration=0;iteration<=iterationNum;iteration++) {
        
        for(int i=0;i<nCols;i++)
            for(int j=0;j<nRows;j++)
                dist[i*nRows+j]=FLT_MAX;
        
        
        int minX,minY,maxX,maxY;
        float D;
        for(int i=0;i<seedNum;i++)        {
            
            
            float Lab_cc = centerL1[i]*centerL1[i] + centerL2[i]*centerL2[i] + centera1[i]*centera1[i] +
                    centera2[i]*centera2[i] + centerb1[i]*centerb1[i] + centerb2[i]*centerb2[i];
            
            float xy_cc = centerx1[i]*centerx1[i] + centerx2[i]*centerx2[i] + centery1[i]*centery1[i] + centery2[i]*centery2[i];
            
            int x=seedArray[i].x;
            int y=seedArray[i].y;
            minX=(x-(StepX)<=0)?0:x-StepX;
            minY=(y-(StepY)<=0)?0:y-StepY;
            maxX=(x+(StepX)>=nCols-1)?nCols-1:x+StepX;
            maxY=(y+(StepY)>=nRows-1)?nRows-1:y+StepY;
            
            
            for(int m=minX;m<=maxX;m++)
                for(int n=minY;n<=maxY;n++)
                    arr[m*nRows+n] = -1;
            
            
            for(int m=minX;m<=maxX;m++)
                for(int n=minY;n<=maxY;n++)
                {
                    
                    //Speed up patch-based distance
                    int pos_i = m*nRows+n;
                    
                    D =  Lab_cc + Lab_2[pos_i] - 2*(centerL1[i]*L1_1[pos_i] + centerL2[i]*L2_1[pos_i] +
                            centera1[i]*a1_1[pos_i] + centera2[i]*a2_1[pos_i] + centerb1[i]*b1_1[pos_i] + centerb2[i]*b2_1[pos_i]);
                    D *= lambda;
                    
                    float D_s = xy_cc + xy_2[pos_i] - 2*(centerx1[i]*x1_1[pos_i] + centerx2[i]*x2_1[pos_i] +
                            centery1[i]*y1_1[pos_i] + centery2[i]*y2_1[pos_i]);
                    
                    D += D_s;
                    
                    
                    float D_tmp;
                    if (bres_color && bres_contour) {
                        D_tmp =  (1-lambda)*dist_bresenham_1(m,n,(int) round(seedArray[i].x),(int) round(seedArray[i].y), nRows,
                                L1_1, L2_1, a1_1, a2_1, b1_1, b2_1,
                                centerL1[i], centerL2[i], centera1[i], centera2[i], centerb1[i], centerb2[i],
                                Lab_cc, xy_cc, Lab_2, xy_2,arr);
                        D += D_tmp;
                        float tmp = dist_bresenham_contour(m,n,(int) round(seedArray[i].x), (int) round(seedArray[i].y),contour,nRows);
                        D *= tmp;
                    }
                    
                    else {
                        
                        if (bres_color && bres_contour == 0) {
                            D_tmp =  (1-lambda)*dist_bresenham_1(m,n,(int) round(seedArray[i].x),(int) round(seedArray[i].y), nRows,
                                    L1_1, L2_1, a1_1, a2_1, b1_1, b2_1,
                                    centerL1[i], centerL2[i], centera1[i], centera2[i], centerb1[i], centerb2[i],
                                    Lab_cc, xy_cc, Lab_2, xy_2,arr);
                            D += D_tmp;
                        }
                        else {
                            if (bres_color == 0 && bres_contour) {
                                float tmp = dist_bresenham_contour(m,n,(int) round(seedArray[i].x), (int) round(seedArray[i].y),contour,nRows);
                                D *= tmp;
                            }
                        }
                        
                        
                    }
                    
                    
                    
                    if(D<dist[m*nRows+n])
                    {
                        label[m*nRows+n]=i;
                        dist[m*nRows+n]=D;
                    }
                }
            
            
        }
        
        
        //Update clusters
        if (iteration<iterationNum) {
            for(int i=0;i<seedNum;i++)
            {
                centerL1[i]=0;
                centerL2[i]=0;
                centera1[i]=0;
                centera2[i]=0;
                centerb1[i]=0;
                centerb2[i]=0;
                centerx1[i]=0;
                centerx2[i]=0;
                centery1[i]=0;
                centery2[i]=0;
                WSum[i]=0;
                clusterSize[i]=0;
                seedArray[i].x=0;
                seedArray[i].y=0;
            }
            
            
            
            for(int i=0;i<nCols;i++)
            {
                for(int j=0;j<nRows;j++)
                {
                    int L=label[i*nRows+j];
                    float Weight=W[i][j];
                    
                    centerL1[L]+=Weight*L1[i][j];
                    centerL2[L]+=Weight*L2[i][j];
                    centera1[L]+=Weight*a1[i][j];
                    centera2[L]+=Weight*a2[i][j];
                    centerb1[L]+=Weight*b1[i][j];
                    centerb2[L]+=Weight*b2[i][j];
                    centerx1[L]+=Weight*x1[i][j];
                    centerx2[L]+=Weight*x2[i][j];
                    centery1[L]+=Weight*y1[i][j];
                    centery2[L]+=Weight*y2[i][j];
                    
                    clusterSize[L]++;
                    WSum[L]+=Weight;
                    seedArray[L].x+=i;
                    seedArray[L].y+=j;
                }
            }
            for(int i=0;i<seedNum;i++)
            {
                WSum[i]=(WSum[i]==0)?1:WSum[i];
                clusterSize[i]=(clusterSize[i]==0)?1:clusterSize[i];
            }
            for(int i=0;i<seedNum;i++)
            {
                centerL1[i]/=WSum[i];
                centerL2[i]/=WSum[i];
                centera1[i]/=WSum[i];
                centera2[i]/=WSum[i];
                centerb1[i]/=WSum[i];
                centerb2[i]/=WSum[i];
                centerx1[i]/=WSum[i];
                centerx2[i]/=WSum[i];
                centery1[i]/=WSum[i];
                centery2[i]/=WSum[i];
                seedArray[i].x/=clusterSize[i];
                seedArray[i].y/=clusterSize[i];
            }
            
            
        }
        
        
        
    }
    
    
    
    //EnforceConnection
    int threshold=(nCols*nRows)/(seedNum*thresholdCoef);
    preEnforceConnectivity(label,nCols,nRows);
    EnforceConnectivity(L1_1,L2_1,a1_1,a2_1,b1_1,b2_1,x1_1,x2_1,y1_1,y2_1,W,label,threshold,nCols,nRows);
    
    
    
    
    free(arr);
    free(min_dist_center);
    
    //Clear Memory
    delete []centerL1;
    delete []centerL2;
    delete []centera1;
    delete []centera2;
    delete []centerb1;
    delete []centerb2;
    delete []centerx1;
    delete []centerx2;
    delete []centery1;
    delete []centery2;
    
    delete []L1_1;
    delete []L1_2;
    delete []L2_1;
    delete []L2_2;
    delete []a1_1;
    delete []a1_2;
    delete []a2_1;
    delete []a2_2;
    delete []b1_1;
    delete []b1_2;
    delete []b2_1;
    delete []b2_2;
    delete []y1_1;
    delete []y1_2;
    delete []x1_1;
    delete []x1_2;
    delete []y2_1;
    delete []y2_2;
    delete []x2_1;
    delete []x2_2;
    delete []Lab_2;
    delete []xy_2;
    delete []WSum;
    delete []clusterSize;
    delete []dist;
    return;
}

#endif
