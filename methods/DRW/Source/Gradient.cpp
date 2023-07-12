//
//  Gradient.cpp
//  UCM
//
//  Created by 朱磊 on 2018/7/3.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#include "DRW.hpp"
using namespace std;
using namespace cv;

//随机生成中心函数
Mat DRW::calMeanM(Mat src)//求平均值矩阵
{
    //多一行一列
    Mat res;
    double k[4] = {0.25, 0.25, 0.25, 0.25};
    Mat Km = Mat(2,2,CV_64F,k);
    filter2D(src, res, src.depth(), Km);
    return res;
}

double DRW::calDis(Vec3b a, Vec3b b)//求欧氏距离
{
    double res = 0;
    for(int i=0; i<3; i++)
    {
        res = res + pow((double(a[i]) - double(b[i])), 2);
    }
    res = log(1 + sqrt(res));
    return res;
}
void DRW::calArea(Mat src, Mat m)//求Area矩阵
{
    int rows = src.rows;
    int cols = src.cols;
    Vec3b p1, p2, p3, p4;//小正方形四个点
    double a, b, c, s, sp; //a ,b ,c三角形三边, s三角形面积， sp三角形半周长
    int idx;
    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<cols; j++)
        {
            idx = i * cols + j;
            s = 0;
            if(i == 0 || j == 0 || i == rows-1 || j == cols -1)
                continue;
            p1 = m.at<Vec3b>(i, j);
            p2 = m.at<Vec3b>(i+1, j);
            p3 = m.at<Vec3b>(i, j+1);
            p4 = m.at<Vec3b>(i+1, j+1);
            
            //第一个小三角形面积
            a = calDis(p1, p3);
            b = calDis(p3, p4);
            c = calDis(p1, p4);
            sp = 0.5 * (a + b + c);
            s = s + sqrt(sp * (sp - a) * (sp - b) * (sp - c));//海伦公式
            //第二个小三角形面积
            a = calDis(p2, p4);
            b = calDis(p1, p2);
            c = calDis(p1, p4);
            sp = 0.5 * (a + b + c);
            s = s + sqrt(sp * (sp - a) * (sp - b) * (sp - c));//海伦公式
            Area[idx] = s;
            pb.at<double>(i, j) = s;
            if(isnan(s))//如果不行计算另一边的三角
            {
                s = 0;
                //第一个小三角形面积
                a = calDis(p1, p2);
                b = calDis(p1, p3);
                c = calDis(p2, p3);
                sp = 0.5 * (a + b + c);
                s = s + sqrt(sp * (sp - a) * (sp - b) * (sp - c));//海伦公式
                //第二个小三角形面积
                a = calDis(p2, p3);
                b = calDis(p2, p4);
                c = calDis(p3, p4);
                sp = 0.5 * (a + b + c);
                s = s + sqrt(sp * (sp - a) * (sp - b) * (sp - c));//海伦公式
                pb.at<double>(i, j) = s;
            }
            if(isnan(s))
                pb.at<double>(i, j - 1) = s;
        }
    }
}

vector<double> DRW::createA(vector<double> area)//获取A矩阵，保存像素i之前Area的总面积
{
    vector<double> A;
    double sum = 0;
    for(int i=0; i<img.rows; i++)
    {
        for(int j=0; j<img.cols; j++)
        {
            sum += Area[i*img.cols+j];
            //cout<<i<<" "<<j<<" "<<area.at<double>(i, j)<<" "<<sum<<endl;
            A.push_back(sum);
        }
    }
    this->sum_area = sum;
    return A;
}

vector<double> DRW::calGradient()
{
    vector<double> A;
    Mat mean, area;
    mean = calMeanM(this->img);
    calArea(this->img, mean);
    A = createA(Area);
    return A;
}

void DRW::createPB()
{
    Mat mean, area;
    mean = calMeanM(this->img);
    calArea(this->img, mean);
    createA(Area);
    normalize(pb, pb, 0.0, 1.0, NORM_MINMAX);
}
