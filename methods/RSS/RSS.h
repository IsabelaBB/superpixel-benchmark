/**
 * Copyright (c) 2016, Dengfeng Chai
 * Contact: chaidf@zju.edu.cn
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef RSS_H
#define RSS_H

#pragma once

#include <cstring>
#include <cmath>
#include <deque>
#include <vector>

using namespace std;

#define maxcost 1024

typedef unsigned char uchar;

class RSS
{
/*
public:
	RSS();
	~RSS();*/
private:

	int nchannels;//channels
	int nrows;//rows
	int ncols;//columns
	int nframes;//for video, supervoxel
	int npixels;
	int nseeds;//number of seeds calculated
	int w;//width of superpixel
	int wz;//for video, supervoxel
	float lambdan;//normalized scale balancing color and distance
	float lambdanz;//normalized scale balancing color and distance

	uchar *dataND; // dataND[npixels*nchannels]
	int *label;

	vector<int> seeds;
	deque<int> cluster[maxcost];

	vector<int> seedsP[4];//4 cpus kernels
	deque<int> clusterP[4][maxcost];

public:
//superpixel
	void Superpixels(uchar *pdata, int *plabel, int rows, int cols, int channels, int K, int func, float lambda, bool perturbseeds, bool flagP);

	void SelectSeeds();
	void PerturbSeeds_Diff();
	void PerturbSeeds_Range();
	void PerformRSS_Diff();
	void PerformRSS_Range();

//parallel version
	void SelectPSeeds();
	void PerturbPSeeds_Diff();
	void PerturbPSeeds_Range();
	void PerformPRSS_Diff();
	void PerformPRSS_Range();

//supervoxel
	void Supervoxels(uchar *pdata, int *plabel, int frames, int rows, int cols, int channels, int K, int func, float lambda, float lambdaz, bool perturbseeds, bool flagP);

	void SelectSeedsV();
	void PerturbSeedsV_Diff();
	void PerturbSeedsV_Range();
	void PerformRSSV_Diff();
	void PerformRSSV_Range();

//parallel version
	void SelectPSeedsV();
	void PerturbPSeedsV_Diff();
	void PerturbPSeedsV_Range();
	void PerformPRSSV_Diff();
	void PerformPRSSV_Range();

};

#endif // RSS_h