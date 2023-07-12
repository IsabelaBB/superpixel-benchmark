//
//  myAlg.cpp
//  fastLRW
//
//  Created by 朱磊 on 2018/3/12.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#include "DRW.hpp"
using namespace std;
using namespace cv;

//优先队列比较函数
bool cmp::operator()(const qNode& one, const qNode& two)
{
    return one.posb < two.posb;//for increasing order of distances
}
bool cmp::operator()(const initNode& one, const initNode& two)
{
    return one.posb < two.posb;//for increasing order of distances
}


double DRW::caldE(Point pos, int lb)//求熵率
{
    double eH = 0, eB = 0;
    double eh = 0;
    int rows = this->img.rows;
    int cols = this->img.cols;
    double N = rows * cols;
    //int go[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    int go[8][2] = {{-1, -1}, {-1, 1}, {1, -1}, {1, 1}, {1,0}, {0,1}, {-1,0}, {0,-1}};
    
    int x = pos.x;
    int y = pos.y;
    int idx = x*cols + y;
    int xp, yp, idxp;
    int count = 0;
    for(int i=0; i<8; i++)
    {
        eh = 0;
        xp = go[i][0] + x;
        yp = go[i][1] + y;
        idxp = xp*cols + yp;
        if(xp < 0 || yp <0 || xp >=rows || yp >= cols)  
            continue;
        if(labels[idxp] != lb)
            continue;
        if(loops[idxp] <= 0 || loops[idx] <= 0)
            continue;
        
        //带权重
        
        //cout<< Area[idx]<<" "<<Area[idxp]<<endl;
        /*
        if(loops[idx] - W[idx][i] / D[idx] > 0)
            eh -= (loops[idx] - W[idx][i] / D[idx]) * log(loops[idx] - W[idx][i] / D[idx]) * exp(-fi * (label_time[idxp] + 1)) / (Area[idx] + 1);
        if(loops[idxp] - W[idx][i] / D[idxp] > 0)
            eh -= (loops[idxp] - W[idx][i] / D[idxp]) * log(loops[idxp] - W[idx][i] / D[idxp]) * exp(-fi * label_time[idxp]) / (Area[idxp] + 1);
        eh += (loops[idx])*log(loops[idx]) * exp(-fi * (label_time[idxp] + 1)) / (Area[idx] + 1);
        eh += (loops[idxp])*log(loops[idxp]) * exp(-fi * label_time[idxp]) / (Area[idxp] + 1);
        eh -= (W[idx][i] / D[idx]) * log(W[idx][i] / D[idx]) * exp(-fi * (label_time[idxp] + 1)) / (Area[idx] + 1);
        eh -= (W[idx][i] / D[idxp]) * log(W[idx][i] / D[idxp]) * exp(-fi * (label_time[idxp])) / (Area[idxp] + 1);
        */
        
        //use pb
        //cout<<pb.at<double>(x, y)<<" "<<pb.at<double>(xp, yp)<<endl;
        
        
        /*
        if(loops[idx] - W[idx][i] / D[idx] > 0)
            eh -= (loops[idx] - W[idx][i] / D[idx]) * log(loops[idx] - W[idx][i] / D[idx]) * exp(-fi * (label_time[idxp] + 1)) / (1 + pb.at<double>(x, y));
        if(loops[idxp] - W[idx][i] / D[idxp] > 0)
            eh -= (loops[idxp] - W[idx][i] / D[idxp]) * log(loops[idxp] - W[idx][i] / D[idxp]) * exp(-fi * label_time[idxp]) / (1 + pb.at<double>(xp, yp));
        eh += (loops[idx])*log(loops[idx]) * exp(-fi * (label_time[idxp] + 1)) / (pb.at<double>(x, y) + 1);
        eh += (loops[idxp])*log(loops[idxp]) * exp(-fi * label_time[idxp]) / (pb.at<double>(xp, yp) + 1);
        eh -= (W[idx][i] / D[idx]) * log(W[idx][i] / D[idx]) * exp(-fi * (label_time[idxp] + 1)) / (1 + pb.at<double>(x, y));
        eh -= (W[idx][i] / D[idxp]) * log(W[idx][i] / D[idxp]) * exp(-fi * (label_time[idxp])) / (1 + pb.at<double>(xp, yp));
        */

        if(loops[idx] - W[idx][i] / D[idx] > 0)
            eh -= (loops[idx] - W[idx][i] / D[idx]) * log(loops[idx] - W[idx][i] / D[idx]) * exp(-fi * (label_time[idxp] + 1));
        if(loops[idxp] - W[idx][i] / D[idxp] > 0)
            eh -= (loops[idxp] - W[idx][i] / D[idxp]) * log(loops[idxp] - W[idx][i] / D[idxp]) * exp(-fi * label_time[idxp]);
        eh += (loops[idx])*log(loops[idx]) * exp(-fi * (label_time[idxp] + 1));
        eh += (loops[idxp])*log(loops[idxp]) * exp(-fi * label_time[idxp]);
        eh -= (W[idx][i] / D[idx]) * log(W[idx][i] / D[idx]) * exp(-fi * (label_time[idxp] + 1));
        eh -= (W[idx][i] / D[idxp]) * log(W[idx][i] / D[idxp]) * exp(-fi * (label_time[idxp]));

        /*
        //不带权重
        if(loops[idx] - W[idx][i] / D[idx] > 0)
            eh -= (loops[idx] - W[idx][i] / D[idx]) * log(loops[idx] - W[idx][i] / D[idx]));
        if(loops[idxp] - W[idx][i] / D[idxp] > 0)
            eh -= (loops[idxp] - W[idx][i] / D[idxp]) * log(loops[idxp] - W[idx][i] / D[idxp]);
        eh += (loops[idx])*log(loops[idx]);
        eh += (loops[idxp])*log(loops[idxp]);
        eh -= (W[idx][i] / D[idx]) * log(W[idx][i] / D[idx]);
        eh -= (W[idx][i] / D[idxp]) * log(W[idx][i] / D[idxp]);
        */
        
        //cout<<exp(-0.05 * label_time[idxp])<<" "<<label_time[idxp]<<endl;
        //cout<<exp(-0.0001 * label_time[idxp])<<" "<<label_time[idxp]<<endl;
        //eH += (eh * exp(-0.0001 * label_time[idxp]));
        eH +=eh;
        count ++;
    }
    eH /= (1 + pb.at<double>(x, y));
    if(isadapt == 0)
        eH = eH / count;
    //double S1 = (sp_num[lb - 1] + 1e-10) / (N / k + 1e-10);
    double S1 = (sp_num[lb - 1] + 1e-10) / double(N + 1e-10);
    double S2 = sp_area[lb - 1] / sum_area;
    double S4 = sp_area[lb - 1] / double(sp_num[lb - 1]);
    if(S4 == 0)
        S4 = 1;
    eB -= S1;
    eB -= S2;
    //eB -= abs(S4 - Area[idx]) / max_area;
    
    //eB -= 1e-1 * ((Area.at<double>(pos.x, pos.y) + 1e-10) / max_area) * log((Area.at<double>(pos.x, pos.y) + 1e-10) / max_area);
    //return eH - 1e-2 * (sp_num[lb - 1]/max_area) * log(sp_num[lb - 1]/max_area);
    //return eH + 1e-1 * eB - 1e-2 * S1;
    //return eH + 10 * eB;
    return eH;
    //return eH + 1e-1 * eB;
    //return eH - 1e-2 * S2;
}

priority_queue<qNode, vector<qNode>, cmp> DRW::initque()//初始化优先队列
{
    int rows = this->img.rows;
    int cols = this->img.cols;
    int x, y, xp, yp, idx, idxp, xpp, ypp, idxpp;
    //int go[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};;//定义邻居位置
    int go[8][2] = {{-1, -1}, {-1, 1}, {1, -1}, {1, 1}, {1,0}, {0,1}, {-1,0}, {0,-1}};
    double e = 0;
    for(int i=0; i<seeds.size(); i++)
    {
        x = seeds[i].x;
        y = seeds[i].y;
        idx = x*cols + y;
        this->labels[idx] = seeds_label[i];//初始化聚类中心类别
        sp_num[seeds_label[i] - 1] += 1;
        sp_area[seeds_label[i] - 1] += Area[idx];
        label_time[idx] = 0;
        //cur_posb.ref<double>(idx, seeds_label[i] - 1) = 1;
        for(int j=0; j<8; j++)//对其邻居节点建立node并入队
        {
            xp = x + go[j][0];
            yp = y + go[j][1];
            idxp = xp * cols + yp;
            if(xp<0 || yp<0 || xp>=rows || yp>=cols)
                continue;
            if(labels[idxp] != 0)
                continue;
            //[idxp] -= W[idx][j]/D[idxp];
            //loops[idx] -= W[idx][j]/D[idx];
            e = caldE(Point(xp, yp), seeds_label[i]);
            //loops[idxp] += W[idx][j]/D[idxp];
            //loops[idx] += W[idx][j]/D[idx];
            qNode qn(Point(xp, yp), Point(x, y), seeds_label[i], e);
            pque.push(qn);
        }
    }
    return this->pque;
}

vector<int> DRW::doDRW(double sigma)//执行RW算法
{
    int rows = this->img.rows;
    int cols = this->img.cols;
    //int go[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    int go[8][2] = {{-1, -1}, {-1, 1}, {1, -1}, {1, 1}, {1,0}, {0,1}, {-1,0}, {0,-1}};
    int c_label, x, y, xp, yp ,idx, idxp, lazyp, rd, mul, cx, cy, idxc;
    int xpp, ypp, idxpp;
    double c_posb, dis, d, nb_pos;
    qNode c_node;
    Point pos;
    Point cpos;
    while(pque.size() != 0)
    {
        c_node = pque.top();
        pque.pop();
        pos = c_node.pos;
        c_label = c_node.label;
        cpos = c_node.cpos;
        x = pos.x;
        y = pos.y;
        idx = x * cols + y;
        cx = cpos.x;
        cy = cpos.y;
        idxc = cx * cols + cy;
        //cout<<c_posb<<" "<<cur_posb.ref<double>(x*cols + y, c_label-1)<<endl;
        if(labels[idx] == 0)//该节点未被标注
        {
            idx = x * cols + y;
            //添加lazy项
            d = D[idx];
            lazyp = (sigma / (d + sigma)) * 100;
            rd = rand()%100 + 1;
            if(rd < lazyp)
                mul = sigma;
            else
                mul = 1;
            labels[idx] = c_label;
            sp_num[c_label - 1] += 1;
            sp_area[c_label - 1] += Area[idx];
            label_time[idx] = sp_num[c_label - 1];
            if(sp_area[c_label - 1] > max_area)
                max_area = sp_area[c_label - 1];
            if(sp_num[c_label - 1] > max_num)
                max_num = sp_num[c_label - 1];
            //loops[idx] -= P[idx][nb_pos];
            //loops[idxc] -= P[idx][nb_pos];
            for(int j=0; j<8; j++)
            {
                xp = x + go[j][0];
                yp = y + go[j][1];
                idxp = xp * cols + yp;
                if(xp<0 || yp<0 || xp>=rows || yp>=cols)
                    continue;
                if(labels[idxp] != c_label)
                    continue;
                loops[idxp] -= W[idx][j] / D[idxp];
                loops[idx] -= W[idx][j] / D[idx];
            }
            for(int j=0; j<8; j++)
            {
                xp = x + go[j][0];
                yp = y + go[j][1];
                idxp = xp * cols + yp;
                if(xp<0 || yp<0 || xp>=rows || yp>=cols)
                    continue;
                if(labels[idxp] != 0)//邻居节点已被标记
                    continue;
                //loops[idxp] -= W[idx][j]/D[idxp];
                //loops[idx] -= W[idx][j]/D[idx];
                double e = caldE(Point(xp, yp), c_label);
                //loops[idxp] += W[idx][j]/D[idxp];
                //loops[idx] += W[idx][j]/D[idx];
                qNode qn(Point(xp, yp), Point(x,y), c_label, mul * e);
                pque.push(qn);//入队
            }
        }
    }
    return labels;
}

vector<Point> DRW::updatseeds()//更新聚类中心
{
    int rows = img.rows;
    int cols = img.cols;
    int N = rows * cols;
    int c_count = this->k;
    int S = N / c_count;//网格面积
    seeds.clear();
    sp_num.clear();
    sp_area.clear();
    seeds_label.clear();
    
    vector<vector<Point>> node_k(real_k);
    vector<double> avg_x(real_k);
    vector<double> avg_y(real_k);
    vector<double> avg_l(real_k);
    vector<double> avg_a(real_k);
    vector<double> avg_b(real_k);
    int idx;
    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<cols; j++)
        {
            idx = i * cols + j;
            int node_label = labels[idx] - 1;
            avg_x[node_label] += i;
            avg_y[node_label] += j;
            avg_l[node_label] += img.at<Vec3b>(i, j)[0];
            avg_a[node_label] += img.at<Vec3b>(i, j)[1];
            avg_b[node_label] += img.at<Vec3b>(i, j)[2];
            node_k[node_label].push_back(Point(i, j));
        }
    }
    for(int i=1; i<=real_k; i++)
    {
        double min = 100000000;
        Point pos;
        //if(node_k[i-1].size() < 10)
        //    continue;
        avg_x[i-1] /= node_k[i-1].size();
        avg_y[i-1] /= node_k[i-1].size();
        avg_l[i-1] /= node_k[i-1].size();
        avg_a[i-1] /= node_k[i-1].size();
        avg_b[i-1] /= node_k[i-1].size();
        for(int j=0; j<node_k[i-1].size(); j++)
        {
            double dis = 0.1 * sqrt(pow(avg_x[i-1] - node_k[i-1][j].x, 2) + pow(avg_y[i-1] - node_k[i-1][j].y, 2));
            dis += sqrt(pow(avg_l[i-1] - img.at<Vec3b>(node_k[i-1][j].x, node_k[i-1][j].y)[0], 2) + pow(avg_a[i-1] - img.at<Vec3b>(node_k[i-1][j].x, node_k[i-1][j].y)[1], 2) + pow(avg_b[i-1] - img.at<Vec3b>(node_k[i-1][j].x, node_k[i-1][j].y)[2], 2));
            if(dis < min)
            {
                min = dis;
                pos = node_k[i-1][j];
            }
        }
        seeds.push_back(pos);
    }
    
    for(int i=0; i<seeds.size(); i++)
    {
        seeds_label.push_back(i + 1);
        idx = seeds[i].x * cols + seeds[i].y;
        labels[idx] = i + 1;
        sp_num.push_back(1);
        sp_area.push_back(1);
    }
    
    labels = vector<int>(N, 0);
    
    this->max_area = 1;
    this->max_num = 1;
    this->label_time = vector<int>(N, 0);
    loops = vector<double>(N, 1);
    real_loops = vector<double>(N, 1);
    
    real_k = seeds.size();
    
    return this->seeds;
}



