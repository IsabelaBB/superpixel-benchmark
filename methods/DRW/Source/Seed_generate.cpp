//
//  Seed-generate.cpp
//  UCM
//
//  Created by 朱磊 on 2018/7/3.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#include "DRW.hpp"
using namespace std;
using namespace cv;

vector<Point> DRW::areaMergeSeeds()//二维网格加三维网格
{
    vector<double> A = calGradient();
    
    
    int rows = img.rows;
    int cols = img.cols;
    int N = rows * cols;
    int S = N / (4 * k);//网格面积
    int dx = sqrt(S * rows / cols);
    int dy = dx * cols / rows;
    if(dx <= 1)
        dx = 2;
    if(dy <= 1)
        dy = 2;
    int r_num = rows / dx;
    int c_num = cols / dy;
    
    double add = A[N-1] / (k);
    int x, y, xp, yp;
    double grad = 0;
    double a = 0;
    int count = 0;
    int x_min = 0, y_min = 0;
    int idxp, idx;
    for(int i=0; i<r_num; i++)
    {
        for(int j=0; j<c_num; j++)
        {
            x = dx*i + dx/2;
            y = dy*j + dy/2;
            if(x<0 || y<0 || x>=rows || y>=cols)
                continue;
            if(a == 0)
            {
                x_min = x;
                y_min = y;
                count = 0;
            }
            for(int xp= x-dx/2; xp<x+dx/2; xp++)
            {
                for(int yp=y-dy/2; yp<y+dy/2; yp++)
                {
                    idxp = xp * cols + yp;
                    if(xp<0 || yp<0 || xp>=rows || yp>=cols)
                    {
                        continue;
                    }
                    a += Area[idxp];
                }
            }
            if(a >= add)
            {
                Point p(x, y);
                seeds.push_back(p);
                a = 0;
            }
            else
                count ++;
        }
    }
    
    
    
    for(int i=0; i<seeds.size(); i++)
    {
        double mingrad = 1000000;
        x = seeds[i].x;
        y = seeds[i].y;
        if(x<0 || y<0 || x>=rows || y>=cols)
            continue;
        for(int go_x=-9; go_x<10; go_x++)
        {
            for(int go_y = -9; go_y<10; go_y++)
            {
                xp = x + go_x;
                yp = y + go_y;
                if(xp<1 || yp<1 || xp>=rows-1 || yp>=cols-1)
                    continue;
                grad = 0;
                for(int j=0; j<3; j++)
                {
                    grad += abs(double(img.at<Vec3b>(xp+1, yp)[j]) + double(img.at<Vec3b>(xp-1, yp)[j]) + double(img.at<Vec3b>(xp, yp+1)[j]) + double(img.at<Vec3b>(xp, yp-1)[j]) - 4 * double(img.at<Vec3b>(xp, yp)[j]));
                    grad += 0.5 * abs(double(img.at<Vec3b>(xp, yp)[j]) - double(img.at<Vec3b>(x, y)[j]));
                }
                if(grad < mingrad)
                {
                    mingrad = grad;
                    seeds[i].x = xp;
                    seeds[i].y = yp;
                }
            }
        }
    }
    
    for(int i=0; i<seeds.size(); i++)
    {
        seeds_label.push_back(i + 1);
        idx = seeds[i].x * cols + seeds[i].y;
        labels[idx] = i + 1;
        sp_num.push_back(0);
        sp_area.push_back(0);
    }
    
    real_k = seeds.size();
    return this->seeds;
}

vector<Point> DRW::createSeeds(vector<double> A, int k)//随机生成种子节点
{
    vector<Point> seeds;
    int N = A.size();
    vector<int> temp(N);
    double r;
    int s = 0;
    int x, y;
    int rows = this->img.rows;
    int cols = this->img.cols;
    while(seeds.size() < k)
    {
        r = rand() % int(A[N-1]);
        //cout<<r<<" "<<A[N-1]<<endl;
        s = 0;
        while(A[s] < r)
            s++;
        if(temp[s] == 0)
        {
            x = (s / cols);
            y = s - (x * cols);
            seeds.push_back(Point(x, y));
            temp[s] = 1;
        }
    }
    return seeds;
}

vector<Point> DRW::randomSeeds()
{
    vector<double> A = calGradient();
    seeds = createSeeds(A, k);
    this->seeds = seeds;
    int idx;
    int rows = img.rows;
    int cols = img.cols;
    for(int i=0; i<seeds.size(); i++)
    {
        seeds_label.push_back(i + 1);
        idx = seeds[i].x * cols + seeds[i].y;
        labels[idx] = i + 1;
        sp_num.push_back(0);
        sp_area.push_back(0);
    }
    real_k = seeds.size();
    return seeds;
}

vector<Point> DRW::gridSeeds()//二维网格加三维网格
{
    int rows = img.rows;
    int cols = img.cols;
    int N = rows * cols;
    int S = N / k;//网格面积
    int dx = sqrt(S * rows / cols);
    int dy = dx * cols / rows;

    int r_num = rows / dx;
    int c_num = cols / dy;
    
    int x, y;
    int idxp, idx;
    int count = 0;
    for(int i=0; i<r_num; i++)
    {
        for(int j=0; j<c_num; j++)
        {
            x = dx*i + dx/2;
            y = dy*j + dy/2;
            if(x<0 || y<0 || x>=rows || y>=cols)
                continue;
            Point p(x, y);
            seeds.push_back(p);
        }
    }
    
    for(int i=0; i<seeds.size(); i++)
    {
        seeds_label.push_back(i + 1);
        idx = seeds[i].x * cols + seeds[i].y;
        labels[idx] = i + 1;
        sp_num.push_back(0);
        sp_area.push_back(0);
    }
    
    real_k = seeds.size();
    return seeds;
}

vector<Point> DRW::pbAreaSeeds(Mat pb)
{
    for(int i=0; i<pb.rows; i++)
    {
        for(int j=0; j<pb.cols; j++)
        {
            Area[i*pb.cols +j] = pb.at<double>(i, j);
        }
    }

    vector<double> A = createA(Area);

    int rows = img.rows;
    int cols = img.cols;
    int N = rows * cols;
    int S = N / (4 * k);//网格面积
    int dx = sqrt(S * rows / cols);
    int dy = dx * cols / rows;
    if(dx <= 1)
        dx = 2;
    if(dy <= 1)
        dy = 2;
    int r_num = rows / dx;
    int c_num = cols / dy;
    
    double add = A[N-1] / (k);
    int x, y, xp, yp;
    double grad = 0;
    double a = 0;
    int count = 0;
    int x_min = 0, y_min = 0;
    int idxp, idx;
    for(int i=0; i<r_num; i++)
    {
        for(int j=0; j<c_num; j++)
        {
            x = dx*i + dx/2;
            y = dy*j + dy/2;
            if(x<0 || y<0 || x>=rows || y>=cols)
                continue;
            if(a == 0)
            {
                x_min = x;
                y_min = y;
                count = 0;
            }
            for(int xp= x-dx/2; xp<x+dx/2; xp++)
            {
                for(int yp=y-dy/2; yp<y+dy/2; yp++)
                {
                    idxp = xp * cols + yp;
                    if(xp<0 || yp<0 || xp>=rows || yp>=cols)
                    {
                        continue;
                    }
                    a += Area[idxp];
                }
            }
            if(a >= add)
            {
                Point p(x, y);
                seeds.push_back(p);
                a = 0;
            }
            else
                count ++;
        }
    }
    
    
    
    for(int i=0; i<seeds.size(); i++)
    {
        double mingrad = 1000000;
        x = seeds[i].x;
        y = seeds[i].y;
        if(x<0 || y<0 || x>=rows || y>=cols)
            continue;
        for(int go_x=-9; go_x<10; go_x++)
        {
            for(int go_y = -9; go_y<10; go_y++)
            {
                xp = x + go_x;
                yp = y + go_y;
                if(xp<1 || yp<1 || xp>=rows-1 || yp>=cols-1)
                    continue;
                grad = 0;
                for(int j=0; j<3; j++)
                {
                    grad += abs(double(img.at<Vec3b>(xp+1, yp)[j]) + double(img.at<Vec3b>(xp-1, yp)[j]) + double(img.at<Vec3b>(xp, yp+1)[j]) + double(img.at<Vec3b>(xp, yp-1)[j]) - 4 * double(img.at<Vec3b>(xp, yp)[j]));
                    grad += 0.5 * abs(double(img.at<Vec3b>(xp, yp)[j]) - double(img.at<Vec3b>(x, y)[j]));
                }
                if(grad < mingrad)
                {
                    mingrad = grad;
                    seeds[i].x = xp;
                    seeds[i].y = yp;
                }
            }
        }
    }
    
    for(int i=0; i<seeds.size(); i++)
    {
        seeds_label.push_back(i + 1);
        idx = seeds[i].x * cols + seeds[i].y;
        labels[idx] = i + 1;
        sp_num.push_back(0);
        sp_area.push_back(0);
    }
    
    real_k = seeds.size();
    return this->seeds;
}

void DRW::givenSeeds(vector<int> water_x, vector<int> water_y)
{
    for(int i=0; i<water_x.size(); i++)
    {
        seeds.push_back(Point(water_x[i],water_y[i]));
        seeds_label.push_back(i + 1);
        sp_num.push_back(1);
        sp_area.push_back(1);
    }
    real_k = seeds.size();
    k = seeds.size();
    //updatseeds();
}
