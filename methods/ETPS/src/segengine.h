/*
    This file is part of spixel.

    Spixel is free software : you can redistribute it and / or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include "stdafx.h"
#include "structures.h"
#include "tsdeque.h"
#include <map>

#define ADD_LEVEL_PARAM_DOUBLE(field, node, name) \
    updateDouble[name] = [](SPSegmentationParameters& params, double val) { params.field = val; };\
    AddLevelParamFromNodeDouble(node, name);

#define ADD_LEVEL_PARAM_INT(field, node, name) \
    updateInt[name] = [](SPSegmentationParameters& params, int val) { params.field = val; };\
    AddLevelParamFromNodeInt(node, name);



using namespace cv;
using namespace std;

void UpdateFromNode(double& val, const FileNode& node);
void UpdateFromNode(int& val, const FileNode& node);
void UpdateFromNode(bool& val, const FileNode& node);

template <typename T> 
void AddLevelParamFromNode(const FileNode& parentNode, const string& nodeName,
    vector<pair<string, vector<T>>>& paramsList)
{
    FileNode n = parentNode[nodeName];

    if (n.empty()) return;

    vector<T> levParams;

    if (n.type() != FileNode::SEQ) {
        levParams.push_back((T)n);
    } else {
        for (FileNode n1 : n) {
            levParams.push_back((T)n1);
        }
    }
    paramsList.push_back(pair<string, vector<T>>(nodeName, levParams));
}


struct SPSegmentationParameters {
    int superpixelNum = 300;        // Number of superpixels (the actual number can be different)
    
    double appWeight = 1.0;
    double regWeight = 1.0;
    double lenWeight = 1.0;
    double sizeWeight = 1.0;        
    double dispWeight = 2000.0;     // \lambda_{disp}
    double smoWeight = 0.2;         // \lambda_{smo}
    double smoWeightCo = 0.2;       // \lambda_{smo}
    double smoWeightHi = 0.2;       // \lambda_{smo}
    double priorWeight = 0.2;       // \lambda_{prior}
    double occPriorWeight = 15.0;   // \lambda_{occ}
    double hiPriorWeight = 5.0;     // \lambda_{hinge}
    double noDisp = 9.0;            // \lambda_{d}
    double inlierThreshold = 3.0;
    int peblThreshold = 0;          // planeEstimationBundaryLengthThreshold 
    double updateThreshold = 0.01;

    int iterations = 1;
    int maxUpdates = 400000;
    int minPixelSize = 1;
    int maxPixelSize = 16;
    int reSteps = 10;

    bool instantBoundary = false;   // Boundary re-estimation on each step of iteration
    bool stereo = false;
    bool computeSGM = false;
    bool batchProcessing = true;
    bool inpaint = false;           // use opencv's inpaint method to fill gaps in
                                    // disparity image
    bool debugOutput = false;
    bool timingOutput = true;
    int randomSeed = 0;

    vector<pair<string, vector<double>>> levelParamsDouble;
    vector<pair<string, vector<int>>> levelParamsInt;
    map<string, function<void(SPSegmentationParameters&, double)>> updateDouble;
    map<string, function<void(SPSegmentationParameters&, int)>> updateInt;

    void SetLevelParams(int level);

    void Read(const FileNode& node)
    {
        UpdateFromNode(superpixelNum, node["superpixelNum"]);
        ADD_LEVEL_PARAM_DOUBLE(appWeight, node, "appWeight");
        ADD_LEVEL_PARAM_DOUBLE(regWeight, node, "regWeight");
        ADD_LEVEL_PARAM_DOUBLE(lenWeight, node, "lenWeight");
        ADD_LEVEL_PARAM_DOUBLE(sizeWeight, node, "sizeWeight");
        ADD_LEVEL_PARAM_DOUBLE(dispWeight, node, "dispWeight");
        ADD_LEVEL_PARAM_DOUBLE(smoWeight, node, "smoWeight");
        ADD_LEVEL_PARAM_DOUBLE(smoWeightCo, node, "smoWeightCo");
        ADD_LEVEL_PARAM_DOUBLE(smoWeightHi, node, "smoWeightHi");
        ADD_LEVEL_PARAM_DOUBLE(priorWeight, node, "priorWeight");
        ADD_LEVEL_PARAM_DOUBLE(occPriorWeight, node, "occPriorWeight");
        ADD_LEVEL_PARAM_DOUBLE(hiPriorWeight, node, "hiPriorWeight");
        UpdateFromNode(noDisp, node["noDisp"]);
        UpdateFromNode(stereo, node["stereo"]);
        UpdateFromNode(computeSGM, node["computeSGM"]);
        UpdateFromNode(batchProcessing, node["batchProcessing"]);
        UpdateFromNode(inpaint, node["inpaint"]);
        UpdateFromNode(instantBoundary, node["instantBoundary"]);
        UpdateFromNode(iterations, node["iterations"]);
        ADD_LEVEL_PARAM_INT(reSteps, node, "reSteps");
        UpdateFromNode(inlierThreshold, node["inlierThreshold"]);
        UpdateFromNode(maxUpdates, node["maxUpdates"]);
        UpdateFromNode(minPixelSize, node["minPixelSize"]);
        UpdateFromNode(maxPixelSize, node["maxPixelSize"]);
        ADD_LEVEL_PARAM_INT(peblThreshold, node, "peblThreshold");
        UpdateFromNode(updateThreshold, node["updateThreshold"]);
        UpdateFromNode(debugOutput, node["debugOutput"]);
        UpdateFromNode(timingOutput, node["timingOutput"]);
        UpdateFromNode(randomSeed, node["randomSeed"]);
        SetLevelParams(0);
    }


private:
    void AddLevelParamFromNodeDouble(const FileNode& parentNode, const string& nodeName)
    {
        AddLevelParamFromNode<double>(parentNode, nodeName, levelParamsDouble);
    }

    void AddLevelParamFromNodeInt(const FileNode& parentNode, const string& nodeName)
    {
        AddLevelParamFromNode<int>(parentNode, nodeName, levelParamsInt);
    }

};

static void read(const FileNode& node, SPSegmentationParameters& x, const SPSegmentationParameters& defaultValue = SPSegmentationParameters());


class BInfoMatrix  {
private:
    BInfo** data;
    int rows;
    int cols;
public:
    BInfoMatrix() : data(nullptr), rows(0), cols(0) { }
    ~BInfoMatrix() { Release(); }

    void Resize(int _rows, int _cols)
    {
        Release();
        rows = _rows;
        cols = _cols;
        data = new BInfo*[rows*cols];
        memset(data, 0, sizeof(BInfo*)*rows*cols);
    }

    BInfo& operator ()(int r, int c) { return *(*data + r*cols + c); }

private:
    void Release()
    {
        if (data != nullptr) {
            BInfo* end = *data + rows*cols, *p = *data;
            for (; p != end; p++) { if (p != nullptr) delete p; }
            delete[] data;
        }
    }
};


class SPSegmentationEngine {
private:
    struct PerformanceInfo {
        double init = 0.0;
        double imgproc = 0.0;
        double ransac = 0.0;
        vector<double> levelTimes;
        vector<int> levelIterations;
        double total = 0.0;
        vector<double> levelMaxEDelta;
        vector<int> levelMaxPixelSize;
    };

    PerformanceInfo performanceInfo;

    // Parameters
    SPSegmentationParameters params;
    double planeSmoothWeight = 1;   // Calculated from params in initialization
    double planeSmoothWeightCo = 0.1;
    double planeSmoothWeightHi = 0.1;
    int initialMaxPixelSize;        // Calculated in initialization

    // Original image to process
    Mat origImg;

    // Depth image (required for stereo)
    Mat1d depthImg;
    Mat1d depthImgAdj;

    // Image to process (in lab color space)
    Mat img;

    // Support structures
    Matrix<Pixel> pixelsImg;    // pixels matrix, dimension varies, depends on level
    Matrix<Pixel*> ppImg;       // matrix of dimension of img, pointers to pixelsImg pixels (for stereo)
    vector<Superpixel*> superpixels;
public:
    SPSegmentationEngine(SPSegmentationParameters params, Mat img, Mat depthImg = Mat());
    virtual ~SPSegmentationEngine();

    void ProcessImage();
    void ProcessImageStereo();
    Mat GetSegmentedImage();
    Mat GetSegmentedImagePlain();
    Mat GetSegmentedImageStereo();
    Mat GetSegmentation() const;
    Mat GetDisparity() const;
    string GetSegmentedImageInfo();
    void PrintDebugInfo();
    void PrintDebugInfo2();
    void PrintDebugInfoStereo();
    void PrintPerformanceInfo();
    int GetNoOfSuperpixels() const;
    double ProcessingTime() { return performanceInfo.total; }
private:
    void Initialize(Superpixel* spGenerator(int));
    void InitializeStereo();
    void InitializeStereoEnergies();
    void InitializePPImage();
    void UpdatePPImage();
    int IterateMoves(int level);
    void ReEstimatePlaneParameters();
    void EstimatePlaneParameters();
    bool SplitPixels(int& newMaxPixelSize);
    void Reset();
    void UpdateBoundaryData();
    void UpdateBoundaryData2();
    void UpdateInlierSums();
    void UpdatePlaneParameters();
    void UpdateStereoSums();
    void UpdateHingeStereoSums();
    void UpdateDispSums();

    void DebugNeighborhoods();
    void DebugBoundary();
    void DebugDispSums();

    bool TryMovePixel(Pixel* p, Pixel* q, PixelMoveData& psd);
    bool TryMovePixelStereo(Pixel* p, Pixel* q, PixelMoveData& psd);

    int Iterate(Deque<Pixel*>& list, Matrix<bool>& inList);
};


SPSegmentationParameters ReadParameters(const string& fileName, const SPSegmentationParameters& defaultValue = SPSegmentationParameters());



