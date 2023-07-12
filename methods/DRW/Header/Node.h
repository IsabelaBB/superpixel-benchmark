//
//  Node.hpp
//  fastLRW
//
//  Created by 朱磊 on 2018/3/9.
//  Copyright © 2018年 朱磊. All rights reserved.
//

#ifndef Node_hpp
#define Node_hpp
using namespace cv;
using namespace std;
class qNode
{
public:
    Point pos;
    Point cpos;
    int label;
    double posb;
    qNode(){};
    qNode(Point p, Point c, int lab, double pb):pos(p),cpos(c),label(lab),posb(pb){};
};

class initNode
{
public:
    Point pos;
    double posb;
    initNode(){};
    initNode(Point p, double pb):pos(p), posb(pb){};
};

class Superpixels
{
public:
    vector<Point> b_pos;
    vector<Point> p_pos;
    set<int> nb;
    int lab;
};

class Boundarys
{
public:
    vector<Point> b_pos;
    double pb;
    int s_lab;
    int e_lab;
};

class TreeNode
{
public:
    int lson;
    int rson;
    double prob;
    vector<Point> p_pos;
};

#endif /* Node_hpp */
