//
//  Auxiliary.cpp
//  UCM
//
//  Created by 朱磊 on 2018/7/3.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#include "DRW.hpp"
using namespace std;
using namespace cv;

//辅助函数
Mat DRW::getoutput(Mat img,Mat labels)
{
    Mat result = img.clone();
    int rows = img.rows;
    int cols = img.cols;
    int go[8][2] = {{1, -1}, {1, 0}, {1, 1}, {0, -1}, {0, 1}, {-1, -1}, {-1, 0}, {-1, 1}};
    for(int i=0; i<img.rows; i++)
    {
        for(int j=0; j<img.cols; j++)
        {
            for(int k=0; k<8; k++)
            {
                int x = i + go[k][0];
                int y = j + go[k][1];
                if(x<0 || y<0 || x>=rows || y>=cols)
                    continue;
                if(labels.at<int>(x, y) != labels.at<int>(i, j))
                {
                    result.at<Vec3b>(i, j)[0] = 0;
                    result.at<Vec3b>(i, j)[1] = 0;
                    result.at<Vec3b>(i, j)[2] = 255;
                    break;
                }
                
            }
        }
    }
    return result;
}

void DRW::delHole()
{
    int rows = img.rows;
    int cols = img.cols;
    vector<vector<Point> > node_k(real_k);
    vector<set<int> > nbor(real_k);
    int x, y, xp, yp, idx, idxp;
    int go[8][2] = {{-1, -1}, {-1, 1}, {1, -1}, {1, 1}, {1,0}, {0,1}, {-1,0}, {0,-1}};
    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<cols; j++)
        {
            x = i;
            y = j;
            idx = i * cols + j;
            int node_label = labels[idx] - 1;
            node_k[node_label].push_back(Point(i, j));
            for(int k=0; k<8; k++)
            {
                xp = x + go[k][0];
                yp = y + go[k][1];
                idxp = xp * cols + yp;
                if(xp<0 || yp<0 || xp>=rows || yp>=cols)
                    continue;
                if(labels[idx] != labels[idxp])
                {
                    nbor[node_label].insert(labels[idxp]);
                }
            }
        }
    }
    for(int i=0; i<real_k; i++)
    {
        if(sp_num[i] <= 20)
        {
            int max = 0;
            int lab = 0;
            set<int>::iterator j = nbor[i].begin();
            while(j != nbor[i].end())
            {
                if(sp_num[*j - 1] > max)
                {
                    lab = *j;
                    max = sp_num[*j - 1];
                }
                j++;
            }
            for(int k=0; k<node_k[i].size(); k++)
            {
                labels[node_k[i][k].x * cols + node_k[i][k].y] = lab;
                node_k[lab - 1].push_back(Point(node_k[i][k].x, node_k[i][k].y));
            }
        }
    }
    
    
    
    //重新编号
    int index = 0;
    vector<int>::iterator itb, ite, it;
    int p;
    for(int i=0; i<real_k; i++)
    {
        if(sp_num[i] <= 5)
            continue;
        itb = labels.begin();
        ite = labels.end();
        while(itb != ite)
        {
            it = find(itb, ite, i+1);
            if(it == ite)
                break;
            p = it - labels.begin();
            labels[p] = index + 1;
            // *it = index + 1;
            itb = it + 1;
        }
        index ++;
    }
    
    
    real_k = index;
    
    
    
}

void DRW::calBoundary()
{
    int rows = img.rows;
    int cols = img.cols;
    int go[8][2] = {{1, -1}, {1, 0}, {1, 1}, {0, -1}, {0, 1}, {-1, -1}, {-1, 0}, {-1, 1}};
    int x, y, xp, yp, idx, idxp;
    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<cols; j++)
        {
            x = i;
            y = j;
            idx = x * cols + y;
            for(int k=0; k<8; k++)
            {
                xp = x + go[k][0];
                yp = y + go[k][1];
                idxp = xp * cols + yp;
                if(xp < 0 || xp >= rows || yp < 0 || yp >= cols)
                    continue;
                if(labels[idx] != labels[idxp])
                {
                    boundary_pixels.push_back(Point(x, y));
                    break;
                }
            }
        }
    }
}

Mat DRW::label2Mat(vector<int> labels)
{
    Mat res = Mat::zeros(img.rows, img.cols, CV_64F);
    for(int i=0; i<img.rows; i++)
    {
        for(int j=0; j<img.cols; j++)
        {
            res.at<int>(i, j) = labels[i*img.cols + j];
        }
    }
    return res;
}

Mat DRW::showseeds(Mat source_img)
{
    Mat result = source_img.clone();
    for(int i=0; i<seeds.size();i++)
    {
        result.at<Vec3b>(seeds[i].x, seeds[i].y)[0] = 0;
        result.at<Vec3b>(seeds[i].x, seeds[i].y)[1] = 0;
        result.at<Vec3b>(seeds[i].x, seeds[i].y)[2] = 255;
    }
    return result;
}
