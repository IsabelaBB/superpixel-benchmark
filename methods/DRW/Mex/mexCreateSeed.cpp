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
    double k = mxGetScalar( prhs[ 1 ] ); //k
    string type_seed = mxArrayToString( prhs[ 2 ] ); //type_seed
    double *inPB;
    if(type_seed == "pb")
    {
        inPB = (double*)mxGetPr(prhs[ 3 ]); //pb;
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
            for (int k = 0; k < channel; k++)
            {
                *(InCurRow + j * channel + (2 - k)) = *(inData + i + j * rows + k * rows * cols);
            }
            if(type_seed == "pb")
                *(InpbRow + j) = *(inPB + i + j * rows);
        }

    }
    
    //----------------------DRW Algorithm-------------------------------
    Mat img = image;
    DRW DRW_alg(img, k);
    if(type_seed == "grid")
        DRW_alg.gridSeeds();
    else if(type_seed == "random")
        DRW_alg.randomSeeds();
    else if(type_seed == "pb")
    	DRW_alg.pbAreaSeeds(pb);
    else
        DRW_alg.areaMergeSeeds();

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
    
    //------------------------Output-------------------------------    

    plhs[0] = mxCreateDoubleMatrix(1, DRW_alg.seeds.size(), mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, DRW_alg.seeds.size(), mxREAL);
    double *out_seed_x = mxGetPr(plhs[0]);
    double *out_seed_y = mxGetPr(plhs[1]);
    for(int i=0; i<DRW_alg.seeds.size(); i++)
    {
        *(out_seed_x + i) = DRW_alg.seeds[i].y + 1;
        *(out_seed_y + i) = DRW_alg.seeds[i].x + 1;
    }
    
}
