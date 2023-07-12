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

#include "RSS.h"

/*
RSS::RSS()
{
}

RSS::~RSS()
{
}*/

void RSS::Superpixels(uchar *pdata, int *plabel, int rows, int cols, int channels, int K, int func, float lambda, bool perturbseeds, bool flagP)
{
	nrows = rows;
	ncols = cols;
	nchannels = channels;
	npixels = nrows*ncols;

    w = (int) (0.5f + sqrt(rows*cols / (float) K));

	lambdan = lambda;
	
	dataND = pdata;
	label = plabel;

    if (func==0)
    {
		if (flagP)
		{
			SelectPSeeds();
			if (perturbseeds)
				PerturbPSeeds_Diff();
			PerformPRSS_Diff();
		}
		else
		{
			SelectSeeds();
			if (perturbseeds)
				PerturbSeeds_Diff();
			PerformRSS_Diff();
		}
    }
    else if (func==1)
    {
		if (flagP)
		{
			SelectPSeeds();
			if (perturbseeds)
				PerturbPSeeds_Range();
			PerformPRSS_Range();
		}
		else
		{
			SelectSeeds();
			if (perturbseeds)
				PerturbSeeds_Range();
			PerformRSS_Range();
		}
    }
}

void RSS::SelectSeeds()
{
	int xstrips, ystrips;
	int xerr, yerr;
	double xerrperstrip, yerrperstrip;
	int xoff, yoff;
	int xe, ye;
	int seedx, seedy;

	xstrips = (0.5 + double(ncols) / double(w));
	ystrips = (0.5 + double(nrows) / double(w));

	xerr = ncols - w*xstrips; if (xerr < 0){ xstrips--; xerr = ncols - w*xstrips; }
	yerr = nrows - w*ystrips; if (yerr < 0){ ystrips--; yerr = nrows - w*ystrips; }

	xerrperstrip = double(xerr) / double(xstrips);
	yerrperstrip = double(yerr) / double(ystrips);

	xoff = w / 2;
	yoff = w / 2;

	seeds.clear();
	
	ye=0;
	seedy=yoff;
	for (int y = 0; y < ystrips; y++)
	{
		xe=0;
		seedx=xoff;
		for (int x = 0; x < xstrips; x++)
		{
			int index = seedy*ncols + seedx;
			seeds.push_back(index);
			xe+=xerrperstrip;
			seedx+=(w+xerrperstrip);
		}
		ye+=yerrperstrip;
		seedy+=(w+yerrperstrip);
	}
	nseeds=seeds.size();
}

//parallel version
void RSS::SelectPSeeds()
{
	int xstrips, ystrips;
	int xerr, yerr;
	double xerrperstrip, yerrperstrip;
	int xoff, yoff;

	xstrips = (0.5 + double(ncols) / double(w));
	ystrips = (0.5 + double(nrows) / double(w));

	xerr = ncols - w*xstrips; if (xerr < 0){ xstrips--; xerr = ncols - w*xstrips; }
	yerr = nrows - w*ystrips; if (yerr < 0){ ystrips--; yerr = nrows - w*ystrips; }

	xerrperstrip = double(xerr) / double(xstrips);
	yerrperstrip = double(yerr) / double(ystrips);

	xoff = w / 2;
	yoff = w / 2;

	for (int c = 0; c < 4; c++)
		seedsP[c].clear();

	int ye = 0;
	int seedy = yoff;
	for (int y = 0; y < ystrips; y++)
	{
		int xe = 0;
		int seedx = xoff;
		for (int x = 0; x < xstrips; x++)
		{
			int index = seedy*ncols + seedx;
			seedsP[(y % 2) * 2 + x % 2].push_back(index);
			xe += xerrperstrip;
			seedx += (w + xerrperstrip);
		}
		ye += yerrperstrip;
		seedy += (w + yerrperstrip);
	}
}

void RSS::PerturbSeeds_Diff()
{
	for (int k = 0; k < nseeds; k++)
	{
		int index = seeds[k];
		int seedy = index / ncols;
		int seedx = index % ncols;
		int newx = seedx;
		int newy = seedy;

		uchar dif, difmin = 255;

		for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
		{
			int centerx = seedx + j;
			int centery = seedy + i;
			if (centerx < 0 || centerx >= ncols) continue;
			if (centery < 0 || centery >= nrows) continue;
			int pixelindex = centery*ncols + centerx;
            uchar *ppixel = dataND+pixelindex*nchannels;

            dif=0;
			for (int ii = -1; ii <= 1; ii++)
			for (int jj = -1; jj <= 1; jj++)
			{
				if (ii == 0 && jj == 0) continue;

				int x = centerx + jj;
				int y = centery + ii;
				if (x < 0 || x >= ncols) continue;
				if (y < 0 || y >= nrows) continue;
				int index_ = y*ncols + x;
				uchar *p=dataND+index_*nchannels;
				
				for (int b = 1; b < nchannels; b++)
				{
                    int d = (p[b]>ppixel[b]) ? p[b]-ppixel[b] : ppixel[b]-p[b];
                    if (dif < d)
                        dif = d;
				}
			}

			if (difmin > dif)
			{
				newx = centerx;
				newy = centery;
				difmin = dif;
			}
		}
		if (newy != seedy || newx != seedx)
			seeds[k] = newy*ncols + newx;
	}
}

void RSS::PerturbPSeeds_Diff()
{
#pragma omp parallel for
	for (int c = 0; c < 4; c++)
	{
		vector<int> &seeds = seedsP[c];
		int nseeds = seeds.size();
		for (int k = 0; k < nseeds; k++)
		{
			int index = seeds[k];
			int seedy = index / ncols;
			int seedx = index % ncols;
			int newx = seedx;
			int newy = seedy;

			uchar dif, difmin = 255;

			for (int i = -1; i <= 1; i++)
			for (int j = -1; j <= 1; j++)
			{
				int centerx = seedx + j;
				int centery = seedy + i;
				if (centerx < 0 || centerx >= ncols) continue;
				if (centery < 0 || centery >= nrows) continue;
				int pixelindex = centery*ncols + centerx;
				uchar *ppixel = dataND+pixelindex*nchannels;

				dif=0;
				for (int ii = -1; ii <= 1; ii++)
				for (int jj = -1; jj <= 1; jj++)
				{
					if (ii == 0 && jj == 0) continue;

					int x = centerx + jj;
					int y = centery + ii;
					if (x < 0 || x >= ncols) continue;
					if (y < 0 || y >= nrows) continue;
					int index_ = y*ncols + x;
					uchar *p=dataND+index_*nchannels;
					
					for (int b = 1; b < nchannels; b++)
					{
						int d = (p[b]>ppixel[b]) ? p[b]-ppixel[b] : ppixel[b]-p[b];
						if (dif < d)
							dif = d;
					}
				}

				if (difmin > dif)
				{
					newx = centerx;
					newy = centery;
					difmin = dif;
				}
			}
			if (newy != seedy || newx != seedx)
				seeds[k] = newy*ncols + newx;
		}
	}
}

void RSS::PerturbSeeds_Range()
{
	for (int k = 0; k < nseeds; k++)
	{
		int index = seeds[k];
		int seedy = index / ncols;
		int seedx = index % ncols;
		int newx = seedx;
		int newy = seedy;

		uchar rang, rangmin = 255;

		for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
		{
			int centerx = seedx + j;
			int centery = seedy + i;
			if (centerx < 0 || centerx >= ncols) continue;
			if (centery < 0 || centery >= nrows) continue;

			uchar *pMax = new uchar[nchannels];
			uchar *pMin = new uchar[nchannels];

			for (int b = 1; b < nchannels; b++)
			{
				pMax[b] = 0;
				pMin[b] = 255;
			}

			for (int ii = -1; ii <= 1; ii++)
			for (int jj = -1; jj <= 1; jj++)
			{
				int x = centerx + jj;
				int y = centery + ii;
				if (x < 0 || x >= ncols) continue;
				if (y < 0 || y >= nrows) continue;
				int index_ = y*ncols + x;

				uchar *p=dataND+index_*nchannels;
				for (int b = 1; b < nchannels; b++)
				{
					if (p[b]>pMax[b]) pMax[b]=p[b];
					if (p[b]<pMin[b]) pMin[b]=p[b];
				}
			}

			rang=pMax[1]-pMin[1];
			for (int b = 2; b < nchannels; b++)
			{
				uchar d = pMax[b] - pMin[b];
				if (rang < d) rang = d;
			}

			if (rangmin > rang)
			{
				newx = centerx;
				newy = centery;
				rangmin = rang;
			}
		}
		if (newy != seedy || newx != seedx)
			seeds[k] = newy*ncols + newx;
	}
}

void RSS::PerturbPSeeds_Range()
{
#pragma omp parallel for
	for (int c = 0; c < 4; c++)
	{
		vector<int> &seeds = seedsP[c];
		int nseeds = seeds.size();
		for (int k = 0; k < nseeds; k++)
		{
			int index = seeds[k];
			int seedy = index / ncols;
			int seedx = index % ncols;
			int newx = seedx;
			int newy = seedy;

			uchar rang, rangmin = 255;

			for (int i = -1; i <= 1; i++)
			for (int j = -1; j <= 1; j++)
			{
				int centerx = seedx + j;
				int centery = seedy + i;
				if (centerx < 0 || centerx >= ncols) continue;
				if (centery < 0 || centery >= nrows) continue;

				uchar *pMax = new uchar[nchannels];
				uchar *pMin = new uchar[nchannels];

				for (int b = 1; b < nchannels; b++)
				{
					pMax[b] = 0;
					pMin[b] = 255;
				}

				for (int ii = -1; ii <= 1; ii++)
				for (int jj = -1; jj <= 1; jj++)
				{
					int x = centerx + jj;
					int y = centery + ii;
					if (x < 0 || x >= ncols) continue;
					if (y < 0 || y >= nrows) continue;
					int index_ = y*ncols + x;

					uchar *p=dataND+index_*nchannels;
					for (int b = 1; b < nchannels; b++)
					{
						if (p[b]>pMax[b]) pMax[b]=p[b];
						if (p[b]<pMin[b]) pMin[b]=p[b];
					}
				}

				rang=pMax[1]-pMin[1];
				for (int b = 2; b < nchannels; b++)
				{
					uchar d = pMax[b] - pMin[b];
					if (rang < d) rang = d;
				}

				if (rangmin > rang)
				{
					newx = centerx;
					newy = centery;
					rangmin = rang;
				}
			}
			if (newy != seedy || newx != seedx)
				seeds[k] = newy*ncols + newx;
		}
	}
}

void RSS::PerformRSS_Diff()
{
	int dx[8] = { 1, 1, 0, -1, -1, -1, 0, 1 };
	int dy[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
	int dindex[8];//shift of 8 neighbors
    uchar *seedcolor = new uchar[nseeds*nchannels];
    int *seedx = new int[nseeds];
    int *seedy = new int[nseeds];
    int *pathDiff = new int[npixels];
	float *cost = new float[npixels];

	for (int i = 0; i < 8; i++)
	{
		dindex[i] = dy[i] * ncols + dx[i];
	}
	
	for (int i = 0; i < npixels; i++)
	{
		cost[i] = maxcost;
	}

	for (int i = 0; i < maxcost; i++)
	{
		cluster[i].clear();
	}

	for (int i = 0; i < nseeds; i++)
	{
		int index = seeds[i];
		int y = index / ncols;
		int x = index % ncols;

        memcpy(seedcolor+i*nchannels, dataND+index*nchannels, nchannels);
		seedx[i] = x;
		seedy[i] = y;
        
        pathDiff[index] = 0;
		cost[index] = 0;
		label[index] = i + 1;
		cluster[0].push_back(index);
	}
    
	for (int i = 0; i < maxcost; i++)
	{
        deque<int> &Clusteri = cluster[i];
        while (Clusteri.size()>0)
        {
            //update the shortest path cost
            int index = Clusteri.front();
            Clusteri.pop_front();
            int y = index / ncols;
            int x = index % ncols;

            uchar *ppixel = dataND+index*nchannels;
            if (ppixel[0] == 255)
            {
                continue;
            }

            ppixel[0] = 255;
            
            int pixell = label[index];

            uchar *pcolor = seedcolor+(pixell-1)*nchannels;

            int sx = seedx[pixell-1];
            int sy = seedy[pixell-1];

            for (int nk = 0; nk < 8; nk++)
            {
                int nx = x + dx[nk];
                int ny = y + dy[nk];

                if (nx < 0 || nx >= ncols) continue;
                if (ny < 0 || ny >= nrows) continue;

                int indexnk = index + dindex[nk];

                uchar *ppixelnk = dataND+indexnk*nchannels;
                if (ppixelnk[0] == 255)
                {
                    continue;
                }

                int difmax=0;
                for (int b = 1; b < nchannels; b++)
                {
                    int dif = (pcolor[b]>ppixelnk[b]) ? pcolor[b]-ppixelnk[b] : ppixelnk[b]-pcolor[b];
                    if (difmax < dif)
                        difmax = dif;
                }
                if (lambdan>0)
                {
                    int dx = (sx>nx) ? sx-nx : nx-sx;
                    int dy = (sy>ny) ? sy-ny : ny-sy;
                    int dif = (dx>dy) ? dx : dy;
                    dif = lambdan * dif / w;
                    if (difmax < dif)
                        difmax = dif;
                }
                if (difmax < pathDiff[index])
                    difmax = pathDiff[index];

                int diffT = difmax;
                if (diffT < cost[indexnk])
                {
                    pathDiff[indexnk] = diffT;
                    cost[indexnk] = diffT;
                    label[indexnk] = pixell;
                    cluster[diffT].push_back(indexnk);
                }
            }
        }
	}
    delete seedcolor;
    delete seedx;
    delete seedy;
    delete pathDiff;
	delete cost;
}

void RSS::PerformPRSS_Diff()
{
	int dx[8] = { 1, 1, 0, -1, -1, -1, 0, 1 };
	int dy[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
	int dindex[8];//shift of 8 neighbors
    uchar **seedcolor = new uchar*[4];
    int **seedx = new int*[4];
    int **seedy = new int*[4];

    int *pathDiff = new int[npixels];
	float *cost = new float[npixels];

#pragma omp parallel for
	for (int i = 0; i < 8; i++)
	{
		dindex[i] = dy[i] * ncols + dx[i];
	}
	
#pragma omp parallel for
	for (int i = 0; i < npixels; i++)
	{
		cost[i] = maxcost;
	}

#pragma omp parallel for
	for (int c = 0; c < 4; c++)
	{
		for (int i = 0; i < maxcost; i++)
		{
			clusterP[c][i].clear();
		}
	}

#pragma omp parallel for
	for (int c = 0; c < 4; c++)
	{
		vector<int> &seeds = seedsP[c];
		int nseeds = seeds.size();
		seedcolor[c] = new uchar[nseeds*nchannels];
		seedx[c] = new int[nseeds];
		seedy[c] = new int[nseeds];
		for (int i = 0; i < nseeds; i++)
		{
			int index = seeds[i];
			int y = index / ncols;
			int x = index % ncols;

			memcpy(seedcolor[c]+i*nchannels, dataND+index*nchannels, nchannels);
			seedx[c][i] = x;
			seedy[c][i] = y;
			
			pathDiff[index] = 0;
			cost[index] = 0;
			label[index] = i * 4 + c + 1;
			clusterP[c][0].push_back(index);
		}
	}

	for (int i = 0; i < maxcost; i++)
	{
#pragma omp parallel for
		for (int c = 0; c < 4; c++)
		{
			//std::cout << c;
			deque<int> &Clusteri = clusterP[c][i];

			while (Clusteri.size()>0)
			{
				//update the shortest path cost
				int index = Clusteri.front();
				Clusteri.pop_front();
				int y = index / ncols;
				int x = index % ncols;

				uchar *ppixel = dataND+index*nchannels;
				if (ppixel[0] == 255)
				{
					continue;
				}

				ppixel[0] = 255;
				
				int pixell = label[index];
				int pixell4 = (pixell-1-c)/4;

				uchar *pcolor = seedcolor[c]+pixell4*nchannels;

				int sx = seedx[c][pixell4];
				int sy = seedy[c][pixell4];

				for (int nk = 0; nk < 8; nk++)
				{
					int nx = x + dx[nk];
					int ny = y + dy[nk];

					if (nx < 0 || nx >= ncols) continue;
					if (ny < 0 || ny >= nrows) continue;

					int indexnk = index + dindex[nk];

					uchar *ppixelnk = dataND+indexnk*nchannels;
					if (ppixelnk[0] == 255)
					{
						continue;
					}

					int difmax=0;
					for (int b = 1; b < nchannels; b++)
					{
						int dif = (pcolor[b]>ppixelnk[b]) ? pcolor[b]-ppixelnk[b] : ppixelnk[b]-pcolor[b];
						if (difmax < dif)
							difmax = dif;
					}
					if (lambdan>0)
					{
						int dx = (sx>nx) ? sx-nx : nx-sx;
						int dy = (sy>ny) ? sy-ny : ny-sy;
						int dif = (dx>dy) ? dx : dy;
						dif = lambdan * dif / w;
						if (difmax < dif)
							difmax = dif;
					}
					if (difmax < pathDiff[index])
						difmax = pathDiff[index];

					int diffT = difmax;
					if (diffT < cost[indexnk])
					{
						pathDiff[indexnk] = diffT;
						cost[indexnk] = diffT;
						label[indexnk] = pixell;
						clusterP[c][diffT].push_back(indexnk);
					}
				}
			}
		}
	}
	for(int i = 0; i < 4; i++)
	{
		delete seedcolor[i];
		delete seedx[i];
		delete seedy[i];
	}	
    delete seedcolor;
    delete seedx;
    delete seedy;
    delete pathDiff;
	delete cost;
}

void RSS::PerformRSS_Range()
{
	int dx[8] = { 1, 1, 0, -1, -1, -1, 0, 1 };
	int dy[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
	int dindex[8];//shift of 8 neighbors
	
    uchar *pathColorMax = new uchar[npixels * nchannels];
    uchar *pathColorMin = new uchar[npixels * nchannels];

    int *pathXYMax = new int[npixels * 2];
    int *pathXYMin = new int[npixels * 2];

	float *cost = new float[npixels];

	for (int i = 0; i < 8; i++)
	{
		dindex[i] = dy[i] * ncols + dx[i];
	}
	
	for (int i = 0; i < npixels; i++)
	{
		cost[i] = maxcost;
	}

	for (int i = 0; i < maxcost; i++)
	{
		cluster[i].clear();
	}

	for (int i = 0; i < nseeds; i++)
	{
		int index = seeds[i];
		int y = index / ncols;
		int x = index % ncols;

        memcpy(pathColorMax+index*nchannels, dataND+index*nchannels, nchannels);
        memcpy(pathColorMin+index*nchannels, dataND+index*nchannels, nchannels);
		
		pathXYMax[index * 2] = x;
		pathXYMax[index * 2 + 1] = y;
		pathXYMin[index * 2] = x;
		pathXYMin[index * 2 + 1] = y;

		cost[index] = 0;
		label[index] = i + 1;
		cluster[0].push_back(index);
	}

	for (int i = 0; i < maxcost; i++)
	{
        deque<int> &Clusteri = cluster[i];
        while (Clusteri.size()>0)
        {
            //update the shortest path cost
            int index = Clusteri.front();
            Clusteri.pop_front();
            int y = index / ncols;
            int x = index % ncols;

            uchar *ppixel = dataND+index*nchannels;
            if (ppixel[0] == 255)
            {
                continue;
            }

            ppixel[0] = 255;
            
            int pixell = label[index];

            uchar *pMax = pathColorMax+index*nchannels;
            uchar *pMin = pathColorMin+index*nchannels;

            int *pXYMax = pathXYMax + index * 2;
            int *pXYMin = pathXYMin + index * 2;

			uchar *pColorMaxnk = new uchar[nchannels];
			uchar *pColorMinnk = new uchar[nchannels];
			int pXYMaxnk[2];
			int pXYMinnk[2];

            for (int nk = 0; nk < 8; nk++)
            {
                int nx = x + dx[nk];
                int ny = y + dy[nk];

                if (nx < 0 || nx >= ncols) continue;
                if (ny < 0 || ny >= nrows) continue;

                int indexnk = index + dindex[nk];
                
                uchar *ppixelnk = dataND+indexnk*nchannels;
                if (ppixelnk[0] == 255)
                {
                    continue;
                }

                int difmax = 0;
                for (int b = 1; b < nchannels; b++)
                {
                    pColorMaxnk[b] = (pMax[b]>ppixelnk[b]) ? pMax[b] : ppixelnk[b];
                    pColorMinnk[b] = (pMin[b]<ppixelnk[b]) ? pMin[b] : ppixelnk[b];
                    int dif = pColorMaxnk[b] - pColorMinnk[b];
                    if (difmax < dif)
                        difmax = dif;
                }

                pXYMaxnk[0] = (pXYMax[0]>nx) ? pXYMax[0] : nx;
                pXYMaxnk[1] = (pXYMax[1]>ny) ? pXYMax[1] : ny;
                pXYMinnk[0] = (pXYMin[0]<nx) ? pXYMin[0] : nx;
                pXYMinnk[1] = (pXYMin[1]<ny) ? pXYMin[1] : ny;

                if (lambdan > 0)
                {
                    int dif = pXYMaxnk[0] - pXYMinnk[0];
                    if (dif < pXYMaxnk[1] - pXYMinnk[1])
                        dif = pXYMaxnk[1] - pXYMinnk[1];;

                    dif = lambdan * dif / w;
                    if (difmax < dif)
                        difmax = dif;
                }

                int diffT = difmax;
                if (diffT < cost[indexnk])
                {
                    cost[indexnk] = diffT;
                    label[indexnk] = pixell;
                    memcpy(pathColorMax + indexnk*nchannels, pColorMaxnk, nchannels);
                    memcpy(pathColorMin + indexnk*nchannels, pColorMinnk, nchannels);

                    pathXYMax[indexnk * 2    ] = pXYMaxnk[0];
                    pathXYMax[indexnk * 2 + 1] = pXYMaxnk[1];
                    pathXYMin[indexnk * 2    ] = pXYMinnk[0];
                    pathXYMin[indexnk * 2 + 1] = pXYMinnk[1];
                    
                    cluster[diffT].push_back(indexnk);
                }
            }
			delete pColorMaxnk;
			delete pColorMinnk;
        }
	}
    delete pathColorMax;
    delete pathColorMin;
    delete pathXYMax;
    delete pathXYMin;
	delete cost;
}

void RSS::PerformPRSS_Range()
{
	int dx[8] = { 1, 1, 0, -1, -1, -1, 0, 1 };
	int dy[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
	int dindex[8];//shift of 8 neighbors
	
    uchar *pathColorMax = new uchar[npixels * nchannels];
    uchar *pathColorMin = new uchar[npixels * nchannels];

    int *pathXYMax = new int[npixels * 2];
    int *pathXYMin = new int[npixels * 2];

	float *cost = new float[npixels];

#pragma omp parallel for
	for (int i = 0; i < 8; i++)
	{
		dindex[i] = dy[i] * ncols + dx[i];
	}
	
#pragma omp parallel for
	for (int i = 0; i < npixels; i++)
	{
		cost[i] = maxcost;
	}

#pragma omp parallel for
	for (int c = 0; c < 4; c++)
	{
		for (int i = 0; i < maxcost; i++)
		{
			clusterP[c][i].clear();
		}
	}
#pragma omp parallel for
	for (int c = 0; c < 4; c++)
	{
		vector<int> &seeds = seedsP[c];
		int nseeds = seeds.size();
		for (int i = 0; i < nseeds; i++)
		{
			int index = seeds[i];
			int y = index / ncols;
			int x = index % ncols;

			memcpy(pathColorMax+index*nchannels, dataND+index*nchannels, nchannels);
			memcpy(pathColorMin+index*nchannels, dataND+index*nchannels, nchannels);
			
			pathXYMax[index * 2] = x;
			pathXYMax[index * 2 + 1] = y;
			pathXYMin[index * 2] = x;
			pathXYMin[index * 2 + 1] = y;

			cost[index] = 0;
			label[index] = i * 4 + c + 1;
			clusterP[c][0].push_back(index);
		}
	}

	for (int i = 0; i < maxcost; i++)
	{
#pragma omp parallel for
		for (int c = 0; c < 4; c++)
		{
			deque<int> &Clusteri = clusterP[c][i];

			while (Clusteri.size()>0)
			{
				//update the shortest path cost
				int index = Clusteri.front();
				Clusteri.pop_front();
				int y = index / ncols;
				int x = index % ncols;

				uchar *ppixel = dataND+index*nchannels;
				if (ppixel[0] == 255)
				{
					continue;
				}

				ppixel[0] = 255;
				
				int pixell = label[index];

				uchar *pMax = pathColorMax+index*nchannels;
				uchar *pMin = pathColorMin+index*nchannels;

				int *pXYMax = pathXYMax + index * 2;
				int *pXYMin = pathXYMin + index * 2;

				uchar *pColorMaxnk = new uchar[nchannels];
				uchar *pColorMinnk = new uchar[nchannels];
				int pXYMaxnk[2];
				int pXYMinnk[2];

				for (int nk = 0; nk < 8; nk++)
				{
					int nx = x + dx[nk];
					int ny = y + dy[nk];

					if (nx < 0 || nx >= ncols) continue;
					if (ny < 0 || ny >= nrows) continue;

					int indexnk = index + dindex[nk];
					
					uchar *ppixelnk = dataND+indexnk*nchannels;
					if (ppixelnk[0] == 255)
					{
						continue;
					}

					int difmax = 0;
					for (int b = 1; b < nchannels; b++)
					{
						pColorMaxnk[b] = (pMax[b]>ppixelnk[b]) ? pMax[b] : ppixelnk[b];
						pColorMinnk[b] = (pMin[b]<ppixelnk[b]) ? pMin[b] : ppixelnk[b];
						int dif = pColorMaxnk[b] - pColorMinnk[b];
						if (difmax < dif)
							difmax = dif;
					}

					pXYMaxnk[0] = (pXYMax[0]>nx) ? pXYMax[0] : nx;
					pXYMaxnk[1] = (pXYMax[1]>ny) ? pXYMax[1] : ny;
					pXYMinnk[0] = (pXYMin[0]<nx) ? pXYMin[0] : nx;
					pXYMinnk[1] = (pXYMin[1]<ny) ? pXYMin[1] : ny;

					if (lambdan > 0)
					{
						int dif = pXYMaxnk[0] - pXYMinnk[0];
						if (dif < pXYMaxnk[1] - pXYMinnk[1])
							dif = pXYMaxnk[1] - pXYMinnk[1];;

						dif = lambdan * dif / w;
						if (difmax < dif)
							difmax = dif;
					}

					int diffT = difmax;
					if (diffT < cost[indexnk])
					{
						cost[indexnk] = diffT;
						label[indexnk] = pixell;
						memcpy(pathColorMax + indexnk*nchannels, pColorMaxnk, nchannels);
						memcpy(pathColorMin + indexnk*nchannels, pColorMinnk, nchannels);

						pathXYMax[indexnk * 2    ] = pXYMaxnk[0];
						pathXYMax[indexnk * 2 + 1] = pXYMaxnk[1];
						pathXYMin[indexnk * 2    ] = pXYMinnk[0];
						pathXYMin[indexnk * 2 + 1] = pXYMinnk[1];
						
						clusterP[c][diffT].push_back(indexnk);
					}
				}
				delete pColorMaxnk;
				delete pColorMinnk;
			}
		}
	}
    delete pathColorMax;
    delete pathColorMin;
    delete pathXYMax;
    delete pathXYMin;
	delete cost;
}
