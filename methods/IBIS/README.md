# IBIS : Iterative Boundaries implicit Identification for Segmentation

This repo. is dedicated to the IBIS super-pixel method submitted at the IEEE CVPRW : CVPM 2018 CVPR workshop.

### IBIS: Real Time Superpixels Segmentation using a Fraction of an Image

#### Abstract : 

> Segmentation is a critical step for many computer vision applications. Among them, the remote photoplethysmography technique is significantly impacted by the quality of region of interest segmentation. With the heart-rate estimation accuracy, the processing time is obviously a key issue for real-time monitoring. Recent face detection algorithms can perform real-time processing, however for unsupervised algorithms, i.e. without any subject detection based on supervised learning, existing methods are not able to achieve real-time on regular platform. In this paper, we propose a new method to perform real-time unsupervised remote photoplethysmograhy based on efficient temporally propagated superpixels segmentation. The proposed method performs the segmentation step by implicitly identifying the superpixel boundaries. Hence, only a fraction of the image is used to perform the segmentation which reduces greatly the computational burden of the process. The segmentation quality remains comparable to state of the art methods while computational time is divided by a factor up to 8 times. The efficiency of the superpixel segmentation allow us to propose a real-time unsupervised rPPG algorithm considering frames of 640x480, RGB, at 25 frames per second on a single core platform. We obtained real-time processing for 93% of precision at 2.5 beat per minute using our inhouse video database.

![alt text](https://github.com/xapha/IBIS/blob/master/intro.png "intro figure")

## IBIS implementation

This implementation is done in C++ and was tested on unix environment.

For benchmark comparison, please refer to the *ibis.h* file for options description.

## Dependences :

You will need *cmake* and *openCV* ( tested with version 3.4.1: https://github.com/opencv/opencv/releases/tag/3.4.1 ) to run this code.
For paralell execution, you will need *openMP* as well.

## Compilation


```Shell Session
git clone "https://github.com/xapha/IBIS.git"
cd IBIS
mkdir build && cd build
mkdir results
cmake -D CMAKE_BUILD_TYPE=Release ..
make
./IBIS
```
Temporal version of IBIS for video stream:
[IBIS temporal](https://github.com/xapha/IBIS_Temporal)
