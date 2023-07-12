//
//  myAlg.hpp
//  fastLRW
//
//  Created by 朱磊 on 2018/3/12.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#ifndef myAlg_hpp
#include "include.h"
#include "Node.h"
#define myAlg_hpp

using namespace std;
using namespace cv;
struct cmp
{
    bool operator()(const qNode& one, const qNode& two);
    bool operator()(const initNode& one, const initNode& two);
};

class DRW
{
public:
    Mat img;
    vector<vector<double>> W;//邻接权重
    vector<vector<double>> P;//转移概率
    vector<double> D;//度
    priority_queue<qNode, vector<qNode>, cmp> pque;//优先队列
    vector<Point> seeds;//聚类中心(seeds节点)
    vector<int> seeds_label;//聚类中心类别(seeds节点)
    vector<double> loops;//节点的度（算法中修改）
    vector<double> real_loops;//节点的度（算法中不变）
    vector<int> sp_num;//各超像素中像素个数
    vector<int> sp_area;//各超像素的曲面面积
    vector<int> labels;//各元素标签
    
    vector<int> label_time;//标注时间
    
    vector<double> Area;//三维曲面面积
    double sum_area;//三维曲面总面积
    int k; //指定超像素个数
    int real_k;//真实超像素个数
    double max_area;//最大超像素曲面面积
    double max_num;//最大超像素个数
    
    double fi;//动态节点所占比重
    double isadapt;
    
    vector<Point> boundary_pixels;
    
    Mat pb;
    
    //构造函数
    DRW(){};
    DRW(Mat img);
    DRW(Mat img, int k);
    DRW(Mat img, int k, double fi);
    DRW(Mat img, int k, double fi, double isadapt);
    
    //DRW函数
    vector<vector<double>> calW(int beta);//生成邻接矩阵
    vector<vector<double>> calP();//生成概率邻接矩阵
    vector<Point> initgrid(Mat img);//网格初始化聚类中心
    double caldE(Point pos, int lb);//求游走熵增
    priority_queue<qNode, vector<qNode>, cmp> initque();//初始化优先队列
    vector<int> doDRW(double sigma);//执行DRW算法,sigma表示自循环概率
    vector<Point> updatseeds();//更新聚类中心
    
    
    
    //新型初始化
    Mat calMeanM(Mat src);//求平局像素值
    double calDis(Vec3b a, Vec3b b);//求距离
    void calArea(Mat src, Mat m);//求各像素曲面面积
    vector<double> createA(vector<double> area);//求总面积
    vector<double> calGradient(); //求梯度信息
    void createPB(); //求PB矩阵
    
    
    //节点初始化
    vector<Point> areaMergeSeeds();//初始化3D2Dseed节点
    vector<Point> gridSeeds();
    void givenSeeds(vector<int> water_x, vector<int> water_y);
    vector<Point> pbAreaSeeds(Mat pb);
    vector<Point> createSeeds(vector<double> A, int k);
    vector<Point> randomSeeds();
    //辅助函数
    void calBoundary();
    Mat label2Mat(vector<int> labels);
    Mat getoutput(Mat img, Mat labels);//获取输出结果
    Mat showseeds(Mat source_img);//获取输出种子节点分布
    void delHole();//去除过小超像素
    void outputMat(Mat m);//输出矩阵信息
};
#endif /* myAlg_hpp */

