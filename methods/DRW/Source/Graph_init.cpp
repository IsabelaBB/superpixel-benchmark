//
//  Graph-init.cpp
//  UCM
//
//  Created by 朱磊 on 2018/7/3.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#include "DRW.hpp"
using namespace std;
using namespace cv;


DRW::DRW(Mat img, int k)//构造函数
{
    int rows = img.rows;
    int cols = img.cols;
    int N = rows * cols;
    this->pb = Mat::zeros(rows, cols, CV_64F);
    this->W = vector<vector<double> >(N, vector<double>(8, 0));
    this->P = vector<vector<double> >(N, vector<double>(8, 0));
    this->D = vector<double>(N, 1);
    this->labels = vector<int>(N, 0);
    this->k = k;
    this->Area = vector<double>(N, 0);
    this->max_area = 1;
    this->max_num = 1;
    
    this->label_time = vector<int>(N, 0);
    this->fi = 0;
    this->isadapt = 1;
    
    loops = vector<double>(N, 1);
    real_loops = vector<double>(N, 1);
    cvtColor(img, this->img, COLOR_BGR2Lab); // lab
}

DRW::DRW(Mat img)//构造函数
{
    int rows = img.rows;
    int cols = img.cols;
    int N = rows * cols;
    this->pb = Mat::zeros(rows, cols, CV_64F);
    this->W = vector<vector<double> >(N, vector<double>(8, 0));
    this->P = vector<vector<double> >(N, vector<double>(8, 0));
    this->D = vector<double>(N, 1);
    this->labels = vector<int>(N, 0);
    this->k = 1000;
    this->Area = vector<double>(N, 0);
    this->max_area = 1;
    this->max_num = 1;
    
    this->fi = 0;
    this->isadapt = 1;
    
    this->label_time = vector<int>(N, 0);
    
    loops = vector<double>(N, 1);
    real_loops = vector<double>(N, 1);
    cvtColor(img, this->img, COLOR_BGR2Lab); // lab
}

DRW::DRW(Mat img, int k, double fi)//构造函数
{
    int rows = img.rows;
    int cols = img.cols;
    int N = rows * cols;
    this->pb = Mat::zeros(rows, cols, CV_64F);
    this->W = vector<vector<double> >(N, vector<double>(8, 0));
    this->P = vector<vector<double> >(N, vector<double>(8, 0));
    this->D = vector<double>(N, 1);
    this->labels = vector<int>(N, 0);
    this->k = k;
    this->Area = vector<double>(N, 0);
    this->max_area = 1;
    this->max_num = 1;
    
    this->label_time = vector<int>(N, 0);
    this->fi = fi;
    this->isadapt = 1;
    
    loops = vector<double>(N, 1);
    real_loops = vector<double>(N, 1);
    cvtColor(img, this->img, COLOR_BGR2Lab); // lab
}

DRW::DRW(Mat img, int k, double fi, double isadapt)//构造函数
{
    int rows = img.rows;
    int cols = img.cols;
    int N = rows * cols;
    this->pb = Mat::zeros(rows, cols, CV_64F);
    this->W = vector<vector<double> >(N, vector<double>(8, 0));
    this->P = vector<vector<double> >(N, vector<double>(8, 0));
    this->D = vector<double>(N, 1);
    this->labels = vector<int>(N, 0);
    this->k = k;
    this->Area = vector<double>(N, 0);
    this->max_area = 1;
    this->max_num = 1;
    
    this->label_time = vector<int>(N, 0);
    this->fi = fi;
    this->isadapt = isadapt;
    
    loops = vector<double>(N, 1);
    real_loops = vector<double>(N, 1);
    cvtColor(img, this->img, COLOR_BGR2Lab); // lab
}



//RW、LRW算法相关函数
vector<vector<double> > DRW::calW(int beta)//生成N*N的邻接矩阵
{
    
    //算W
    int rows = img.rows;
    int cols = img.cols;
    int N = rows * cols;
    //int go[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};//定义邻居位置
    int go[8][2] = {{-1, -1}, {-1, 1}, {1, -1}, {1, 1}, {1,0}, {0,1}, {-1,0}, {0,-1}};
    int x, y, idxp, xp, yp, idx;
    double val_w, d, min, max;
    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<cols; j++)
        {
            x = i;
            y = j;
            idx = x*cols + y;
            Vec3b lab = img.at<Vec3b>(x, y);
            d = 0;
            min = 1;
            max = 0;
            for(int k=0; k<8; k++)
            {
                xp = x + go[k][0];
                yp = y + go[k][1];
                idxp = xp*cols + yp;
                if(xp<0 || yp<0 || xp>=rows || yp>=cols)
                    continue;
                Vec3b plab = img.at<Vec3b>(xp, yp);
                val_w = sqrt((pow((double(lab[0]) - double(plab[0]))/100.0, 2) + pow((double(lab[1]) - double(plab[1]))/255.0, 2) + pow((double(lab[2]) - double(plab[2]))/255.0, 2) ) / 3.0);//lab
                W[idx][k] = exp(-beta * val_w);
                d = d + W[idx][k];
            }
            D[idx] = d;//对于rw需要多+sigma
        }
    }
    this->D = D;
    this->W = W;
    return W;
}

vector<vector<double> > DRW::calP()//生成概率邻接矩阵
{
    int rows = img.rows;
    int cols = img.cols;
    int N = rows * cols;
    int go[8][2] = {{-1, -1}, {-1, 1}, {1, -1}, {1, 1}, {1,0}, {0,1}, {-1,0}, {0,-1}};
    int x, y, idxp, xp, yp, idx;
    double val_w, d, min, max;
    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<cols; j++)
        {
            x = i;
            y = j;
            idx = x*cols + y;
            
            for(int k=0; k<8; k++)
            {
                P[idx][k] = W[idx][k] / D[idx];
            }
            loops[idx] = 1;
            real_loops[idx] = 1;
        }
    }
    //P = W;
    return this->P;
}
