//
//  mexFunction.cpp
//  gradDRW
//
//  Created by ?±ç? on 2018/3/25.
//  Copyright ?© 2018å¹??±ç?. All rights reserved.
//

#include "include.h"
#include "Node.h"
#include "DRW.hpp"
#include "mex.h"
#include <cstdio>
using namespace std;
using namespace cv;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    //----------------------Input Params-------------------------------
    uchar *inData = (uchar*)mxGetPr(prhs[0]); //img
    double *in_seed_x = (double*)mxGetPr(prhs[1]);
    double *in_seed_y = (double*)mxGetPr(prhs[2]);
    double k = max(mxGetM(prhs[ 1 ]), mxGetN(prhs[ 1 ]));
    double *inPB;
    double beta, sigma;
    int iter, isadapt;
    string type_seed, type_grad;
    switch(nrhs)
    {
        case 3:
        {
            beta = 1e-4; //beta
            sigma = 50; //sigma
            iter = 0;
            type_grad = "area";
            isadapt = 1;
            break;
        }
        case 5:
        {
            beta = (double)mxGetScalar( prhs[ 3 ] ); //beta
            sigma = (double)mxGetScalar( prhs[ 4 ] ); //sigma
            iter = 0;
            type_grad = "area";
            isadapt = 1;
            break;
        }
        case 6:
        {
            beta = (double)mxGetScalar( prhs[ 3 ] ); //beta
            sigma = (double)mxGetScalar( prhs[ 4 ] ); //sigma
            iter = (int)mxGetScalar( prhs[ 5 ] ); //iter
            type_grad = "area";
            isadapt = 1;
            break;
        }
        case 9:
        {
            beta = (double)mxGetScalar( prhs[ 3 ] ); //beta
            sigma = (double)mxGetScalar( prhs[ 4 ] ); //sigma
            iter = (int)mxGetScalar( prhs[ 5 ] ); //iter
            inPB = (double*)mxGetPr(prhs[ 6 ]); //pb
            type_grad = mxArrayToString( prhs[ 7 ] ); //type_grad
            isadapt = (int)mxGetScalar( prhs[ 8 ] ); //isadapt   
            break;
        }
        default:
        {
            mexErrMsgTxt("Input Error, Please Cheak Your Input Params!");
        }
    }

    //----------------------Read Image&pb-------------------------------
    int rows = mxGetM(prhs[0]); 
    int cols = mxGetN(prhs[0]); 
    int channel = mxGetNumberOfDimensions(prhs[0]);
    
    cols = cols/channel;
    Mat image(rows,cols,CV_8UC3);//CV_8UC3
    Mat pb(rows,cols,CV_64F);//CV_8UC3
    uchar *InCurRow;
    double *InpbRow;
    for (int i = 0; i < rows; i++)
    {
        InCurRow = (uchar*)image.ptr<uchar>(i);
        InpbRow = (double*)pb.ptr<double>(i);
        for (int j = 0; j < cols; j++)
        {
            for (int l = 0; l < channel; l++)
            {
                *(InCurRow + j * channel + (2 - l)) = *(inData + i + j * rows + l * rows * cols);
            }
            if(nrhs > 6)
                *(InpbRow + j) = *(inPB + i + j * rows);
        }
    }

    //----------------------Read Seed-------------------------------
	vector<int> seed_x, seed_y;
    for(int i=0; i<k; i++)
    {
    	if(in_seed_y[i] - 1 >= rows || in_seed_x[i] - 1 >= cols)
    	{
    		mexErrMsgTxt("Seed Chosen Error, Please Check The Input Seeds");
    	}
    	seed_x.push_back(in_seed_y[i] - 1);
    	seed_y.push_back(in_seed_x[i] - 1);
    }

    //----------------------DRW Algorithm-------------------------------
    Mat img = image;
    DRW DRW_alg(img, k, beta, isadapt);
    DRW_alg.calW(sigma);
    DRW_alg.calP();
    DRW_alg.givenSeeds(seed_x, seed_y);

    /*
    if(nrhs == 8)
    {
        double *wt_x = (double*)mxGetPr(prhs[5]);
        double *wt_y = (double*)mxGetPr(prhs[6]);
        vector<int> water_x, water_y;
        float numsup = mxGetScalar( prhs[7] );
        for(int i=0; i<numsup; i++)
        {
            water_x.push_back(wt_x[i] - 1);
            water_y.push_back(wt_y[i] - 1);
        }
        L_alg.initWater(water_x, water_y);
    }
    */
    
    
    
    
    if(type_grad == "pb")
        DRW_alg.pb = pb;
    else
        DRW_alg.calGradient();

    DRW_alg.initque();
    DRW_alg.doDRW(1e-7);
    DRW_alg.delHole();


    for(int i=0 ;i<iter; i++)
    {
        DRW_alg.updatseeds();
        DRW_alg.initque();
        DRW_alg.doDRW(1e-7);
        DRW_alg.delHole();
    }

    
    Mat labels = DRW_alg.label2Mat(DRW_alg.labels);

           
    

    
    //------------------------Output-------------------------------    
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, DRW_alg.boundary_pixels.size(), mxREAL);
    double *outData = mxGetPr(plhs[0]);
    Mat out;//CV_8UC3
    out = image;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            *(outData + i + j * rows) = (double)labels.at<int>(i, j);
        }
    }
}
