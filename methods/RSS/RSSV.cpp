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

void RSS::Supervoxels(uchar *pdata, int *plabel, int frames, int rows, int cols, int channels, int K, int func, float lambda, float lambdaz, bool perturbseeds, bool flagP)
{
	nframes = frames;
	nrows = rows;
	ncols = cols;
	nchannels = channels;
	npixels = nframes*nrows*ncols;

    w = (int) (0.5f + pow(npixels / float(K), 1.0 / 3.0));
	wz = w;

	lambdan = lambda;
	lambdanz = lambdaz;
	
	dataND = pdata;
	label = plabel;

    if (func==0)
    {
		if (flagP)
		{
			SelectPSeedsV();
			if (perturbseeds)
				PerturbPSeedsV_Diff();
			PerformPRSSV_Diff();
		}
		else
		{
			SelectSeedsV();
			if (perturbseeds)
				PerturbSeedsV_Diff();
			PerformRSSV_Diff();
		}
    }
    else if (func==1)
    {
		if (flagP)
		{
			SelectPSeedsV();
			if (perturbseeds)
				PerturbPSeedsV_Range();
			PerformPRSSV_Range();
		}
		else
		{
			SelectSeedsV();
			if (perturbseeds)
				PerturbSeedsV_Range();
			PerformRSSV_Range();
		}
    }
}

void RSS::SelectSeedsV()
{
	int xstrips, ystrips, zstrips;
	int xerr, yerr, zerr;
	double xerrperstrip, yerrperstrip, zerrperstrip;
	int xoff, yoff, zoff;

	xstrips = (0.5 + double(ncols) / double(w));
	ystrips = (0.5 + double(nrows) / double(w));
	zstrips = (0.5 + double(nframes) / double(wz));

	xerr = ncols - w*xstrips; if (xerr < 0){ xstrips--; xerr = ncols - w*xstrips; }
	yerr = nrows - w*ystrips; if (yerr < 0){ ystrips--; yerr = nrows - w*ystrips; }
	zerr = nframes - wz*zstrips; if (zerr < 0){ zstrips--; zerr = nframes - wz*zstrips; }

	xerrperstrip = double(xerr) / double(xstrips);
	yerrperstrip = double(yerr) / double(ystrips);
	zerrperstrip = double(zerr) / double(zstrips);

	xoff = w / 2;
	yoff = w / 2;
	zoff = wz / 2;

	seeds.clear();

	int ze=0;
	int seedz=zoff;
	for (int z = 0; z < zstrips; z++)
	{
		int ye=0;
		int seedy=yoff;
		for (int y = 0; y < ystrips; y++)
		{
			int xe=0;
			int seedx=xoff;
			for (int x = 0; x < xstrips; x++)
			{
				int index = (seedz*nrows + seedy)*ncols + seedx;
				seeds.push_back(index);
				xe+=xerrperstrip;
				seedx+=(w+xerrperstrip);
			}
			ye+=yerrperstrip;
			seedy+=(w+yerrperstrip);
		}
		ze+=zerrperstrip;
		seedz+=(w+zerrperstrip);
	}
	nseeds=seeds.size();
}

//parallel version
void RSS::SelectPSeedsV()
{
	int xstrips, ystrips, zstrips;
	int xerr, yerr, zerr;
	double xerrperstrip, yerrperstrip, zerrperstrip;
	int xoff, yoff, zoff;

	xstrips = (0.5 + double(ncols) / double(w));
	ystrips = (0.5 + double(nrows) / double(w));
	zstrips = (0.5 + double(nframes) / double(wz));

	xerr = ncols - w*xstrips; if (xerr < 0){ xstrips--; xerr = ncols - w*xstrips; }
	yerr = nrows - w*ystrips; if (yerr < 0){ ystrips--; yerr = nrows - w*ystrips; }
	zerr = nframes - wz*zstrips; if (zerr < 0){ zstrips--; zerr = nframes - wz*zstrips; }

	xerrperstrip = double(xerr) / double(xstrips);
	yerrperstrip = double(yerr) / double(ystrips);
	zerrperstrip = double(zerr) / double(zstrips);

	xoff = w / 2;
	yoff = w / 2;
	zoff = wz / 2;

	for (int c = 0; c < 4; c++)
		seedsP[c].clear();

	int ze=0;
	int seedz=zoff;
	for (int z = 0; z < zstrips; z++)
	{
		int ye=0;
		int seedy=yoff;
		for (int y = 0; y < ystrips; y++)
		{
			int xe=0;
			int seedx=xoff;
			for (int x = 0; x < xstrips; x++)
			{
				int index = (seedz*nrows + seedy)*ncols + seedx;
				seedsP[(y % 2) * 2 + x % 2].push_back(index);
				xe+=xerrperstrip;
				seedx+=(w+xerrperstrip);
			}
			ye+=yerrperstrip;
			seedy+=(w+yerrperstrip);
		}
		ze+=zerrperstrip;
		seedz+=(w+zerrperstrip);
	}
}

void RSS::PerturbSeedsV_Diff()
{
	for (int k = 0; k < nseeds; k++)
	{
		int index = seeds[k];
		int seedx = index % ncols;
		int seedy = index / ncols;
		int seedz = seedy / nrows;
		seedy = seedy % nrows;
		int nx = seedx;
		int ny = seedy;
		int nz = seedz;

		uchar dif, difmin = 255;

		for (int t = -1; t <= 1; t++)
		for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
		{
			int centerx = seedx + j;
			int centery = seedy + i;
			int centerz = seedz + t;
			if (centerx < 0 || centerx >= ncols) continue;
			if (centery < 0 || centery >= nrows) continue;
			if (centerz < 0 || centerz >= nframes) continue;
			int pixelindex = (centerz*nrows+centery)*ncols + centerx;

            uchar *ppixel = dataND+pixelindex*nchannels;

            dif=0;
			for (int tt = -1; tt <= 1; tt++)
			for (int ii = -1; ii <= 1; ii++)
			for (int jj = -1; jj <= 1; jj++)
			{
				if (tt == 0 && ii == 0 && jj == 0) continue;
				int x = centerx + jj;
				int y = centery + ii;
				int z = centerz + tt;
				if (x < 0 || x >= ncols) continue;
				if (y < 0 || y >= nrows) continue;
				if (z < 0 || z >= nframes) continue;
				int index_ = (z*nrows+y)*ncols + x;

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
				nx = centerx;
				ny = centery;
				nz = centerz;
				difmin = dif;
			}
		}
		if (nz != seedz || ny != seedy || nx != seedx)
			seeds[k] = (nz*nrows+ ny)*ncols + nx;
	}
}

//parallel version
void RSS::PerturbPSeedsV_Diff()
{
#pragma omp parallel for
	for (int c = 0; c < 4; c++)
	{
		vector<int> &seeds = seedsP[c];
		int nseeds = seeds.size();
		for (int k = 0; k < nseeds; k++)
		{
			int index = seeds[k];
			int seedx = index % ncols;
			int seedy = index / ncols;
			int seedz = seedy / nrows;
			seedy = seedy % nrows;
			int nx = seedx;
			int ny = seedy;
			int nz = seedz;

			uchar dif, difmin = 255;

			for (int t = -1; t <= 1; t++)
			for (int i = -1; i <= 1; i++)
			for (int j = -1; j <= 1; j++)
			{
				int centerx = seedx + j;
				int centery = seedy + i;
				int centerz = seedz + t;
				if (centerx < 0 || centerx >= ncols) continue;
				if (centery < 0 || centery >= nrows) continue;
				if (centerz < 0 || centerz >= nframes) continue;
				int pixelindex = (centerz*nrows+centery)*ncols + centerx;

				uchar *ppixel = dataND+pixelindex*nchannels;

				dif=0;
				for (int tt = -1; tt <= 1; tt++)
				for (int ii = -1; ii <= 1; ii++)
				for (int jj = -1; jj <= 1; jj++)
				{
					if (tt == 0 && ii == 0 && jj == 0) continue;
					int x = centerx + jj;
					int y = centery + ii;
					int z = centerz + tt;
					if (x < 0 || x >= ncols) continue;
					if (y < 0 || y >= nrows) continue;
					if (z < 0 || z >= nframes) continue;
					int index_ = (z*nrows+y)*ncols + x;

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
					nx = centerx;
					ny = centery;
					nz = centerz;
					difmin = dif;
				}
			}
			if (nz != seedz || ny != seedy || nx != seedx)
				seeds[k] = (nz*nrows+ ny)*ncols + nx;
		}
	}
}

void RSS::PerturbSeedsV_Range()
{
	for (int k = 0; k < nseeds; k++)
	{
		int index = seeds[k];
		int seedx = index % ncols;
		int seedy = index / ncols;
		int seedz = seedy / nrows;
		seedy = seedy % nrows;
		int nx = seedx;
		int ny = seedy;
		int nz = seedz;

		uchar rang, rangmin = 255;

		for (int t = -1; t <= 1; t++)
		for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
		{
			int centerx = seedx + j;
			int centery = seedy + i;
			int centerz = seedz + t;
			if (centerx < 0 || centerx >= ncols) continue;
			if (centery < 0 || centery >= nrows) continue;
			if (centerz < 0 || centerz >= nframes) continue;

			uchar *pMax = new uchar[nchannels];
			uchar *pMin = new uchar[nchannels];

			for (int b = 1; b < nchannels; b++)
			{
				pMax[b] = 0;
				pMin[b] = 255;
			}

			for (int tt = -1; tt <= 1; tt++)
			for (int ii = -1; ii <= 1; ii++)
			for (int jj = -1; jj <= 1; jj++)
			{
				int x = centerx + jj;
				int y = centery + ii;
				int z = centerz + tt;
				if (x < 0 || x >= ncols) continue;
				if (y < 0 || y >= nrows) continue;
				if (z < 0 || z >= nframes) continue;
				int index_ = (z*nrows+y)*ncols + x;

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
				nx = centerx;
				ny = centery;
				nz = centerz;
				rangmin = rang;
			}
		}
		if (nz != seedz || ny != seedy || nx != seedx)
			seeds[k] = (nz*nrows+ ny)*ncols + nx;
	}
}

void RSS::PerturbPSeedsV_Range()
{
#pragma omp parallel for
	for (int c = 0; c < 4; c++)
	{
		vector<int> &seeds = seedsP[c];
		int nseeds = seeds.size();
		for (int k = 0; k < nseeds; k++)
		{
			int index = seeds[k];
			int seedx = index % ncols;
			int seedy = index / ncols;
			int seedz = seedy / nrows;
			seedy = seedy % nrows;
			int nx = seedx;
			int ny = seedy;
			int nz = seedz;

			uchar rang, rangmin = 255;

			for (int t = -1; t <= 1; t++)
			for (int i = -1; i <= 1; i++)
			for (int j = -1; j <= 1; j++)
			{
				int centerx = seedx + j;
				int centery = seedy + i;
				int centerz = seedz + t;
				if (centerx < 0 || centerx >= ncols) continue;
				if (centery < 0 || centery >= nrows) continue;
				if (centerz < 0 || centerz >= nframes) continue;

				uchar *pMax = new uchar[nchannels];
				uchar *pMin = new uchar[nchannels];

				for (int b = 1; b < nchannels; b++)
				{
					pMax[b] = 0;
					pMin[b] = 255;
				}

				for (int tt = -1; tt <= 1; tt++)
				for (int ii = -1; ii <= 1; ii++)
				for (int jj = -1; jj <= 1; jj++)
				{
					int x = centerx + jj;
					int y = centery + ii;
					int z = centerz + tt;
					if (x < 0 || x >= ncols) continue;
					if (y < 0 || y >= nrows) continue;
					if (z < 0 || z >= nframes) continue;
					int index_ = (z*nrows+y)*ncols + x;

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
					nx = centerx;
					ny = centery;
					nz = centerz;
					rangmin = rang;
				}
			}
			if (nz != seedz || ny != seedy || nx != seedx)
				seeds[k] = (nz*nrows+ ny)*ncols + nx;
		}
	}
}

void RSS::PerformRSSV_Diff()
{
    uchar *seedcolor = new uchar[nseeds*nchannels];
    int *seedx = new int[nseeds];
    int *seedy = new int[nseeds];
    int *seedz = new int[nseeds];
    int *pathcost = new int[npixels];
	float *cost = new float[npixels];

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
		int x = index % ncols;
		int y = index / ncols;
		int z = y / nrows;
		y = y % nrows;

        memcpy(seedcolor+i*nchannels, dataND+index*nchannels, nchannels);
		seedx[i] = x;
		seedy[i] = y;
		seedz[i] = z;
        
        pathcost[index] = 0;
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
			int x = index % ncols;
			int y = index / ncols;
			int z = y / nrows;
			y = y % nrows;

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
            int sz = seedz[pixell-1];

			for (int dz = -1; dz < 2; dz++)
			for (int dy = -1; dy < 2; dy++)
			for (int dx = -1; dx < 2; dx++)
			{
				if (dx == 0 && dy == 0 && dz == 0) continue;

				int nx = x + dx;
				int ny = y + dy;
				int nz = z + dz;
				if (nx < 0 || nx >= ncols) continue;
				if (ny < 0 || ny >= nrows) continue;
				if (nz < 0 || nz >= nframes) continue;

				int indexnk = (nz*nrows + ny)*ncols + nx;

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
                    int dz = (sz>nz) ? sz-nz : nz-sz;
                    int dif = (dx>dy) ? dx : dy;
                    dif = lambdan * dif / w;
                    int difz = lambdanz * dz / wz;

                    if (dif < difz)
                        dif = difz;
                    if (difmax < dif)
                        difmax = dif;
                }
                if (difmax < pathcost[index])
                    difmax = pathcost[index];

                int diffT = difmax;
                if (diffT < cost[indexnk])
                {
                    pathcost[indexnk] = diffT;
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
    delete seedz;
    delete pathcost;
	delete cost;
}

void RSS::PerformPRSSV_Diff()
{
    uchar **seedcolor = new uchar*[4];
    int **seedx = new int*[4];
    int **seedy = new int*[4];
    int **seedz = new int*[4];

    int *pathDiff = new int[npixels];
	float *cost = new float[npixels];

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
		seedz[c] = new int[nseeds];
		for (int i = 0; i < nseeds; i++)
		{
			int index = seeds[i];
			int x = index % ncols;
			int y = index / ncols;
			int z = y / nrows;
			y = y % nrows;

			memcpy(seedcolor[c]+i*nchannels, dataND+index*nchannels, nchannels);
			seedx[c][i] = x;
			seedy[c][i] = y;
			seedy[c][i] = z;
			
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
			deque<int> &Clusteri = clusterP[c][i];

			while (Clusteri.size()>0)
			{
				//update the shortest path cost
				int index = Clusteri.front();
				Clusteri.pop_front();
				int x = index % ncols;
				int y = index / ncols;
				int z = y / nrows;
				y = y % nrows;

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
				int sz = seedz[c][pixell4];

				for (int dz = -1; dz < 2; dz++)
				for (int dy = -1; dy < 2; dy++)
				for (int dx = -1; dx < 2; dx++)
				{
					if (dx == 0 && dy == 0 && dz == 0) continue;

					int nx = x + dx;
					int ny = y + dy;
					int nz = z + dz;
					if (nx < 0 || nx >= ncols) continue;
					if (ny < 0 || ny >= nrows) continue;
					if (nz < 0 || nz >= nframes) continue;

					int indexnk = (nz*nrows + ny)*ncols + nx;

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
						int dz = (sz>nz) ? sz-nz : nz-sz;
						int dif = (dx>dy) ? dx : dy;
						dif = lambdan * dif / w;
						int difz = lambdanz * dz / wz;

						if (dif < difz)
							dif = difz;
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
		delete seedz[i];
	}	
    delete seedcolor;
    delete seedx;
    delete seedy;
    delete seedz;
    delete pathDiff;
	delete cost;
}

void RSS::PerformRSSV_Range()
{
    uchar *pathColorMax = new uchar[npixels * nchannels];
    uchar *pathColorMin = new uchar[npixels * nchannels];
    uchar *pColorMaxnk = new uchar[nchannels];
    uchar *pColorMinnk = new uchar[nchannels];

    int *pathXYZMax = new int[npixels * 3];
    int *pathXYZMin = new int[npixels * 3];
	int pXYZMaxnk[3];
	int pXYZMinnk[3];

	float *cost = new float[npixels];

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
		int x = index % ncols;
		int y = index / ncols;
		int z = y / nrows;
		y = y % nrows;

		memcpy(pathColorMax+index*nchannels, dataND+index*nchannels, nchannels);
		memcpy(pathColorMin+index*nchannels, dataND+index*nchannels, nchannels);
		
		pathXYZMax[index * 3] = x;
		pathXYZMax[index * 3 + 1] = y;
		pathXYZMax[index * 3 + 2] = z;
		pathXYZMin[index * 3] = x;
		pathXYZMin[index * 3 + 1] = y;
		pathXYZMin[index * 3 + 2] = z;

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
			int x = index%ncols;
			int y = index/ncols;
			int z = y/nrows;
			y = y%nrows;

			uchar *ppixel = dataND+index*nchannels;
			if (ppixel[0] == 255)
			{
				continue;
			}

			ppixel[0] = 255;
			
			int pixell = label[index];

			uchar *pMax = pathColorMax+index*nchannels;
			uchar *pMin = pathColorMin+index*nchannels;

			int *pXYZMax = pathXYZMax + index * 3;
			int *pXYZMin = pathXYZMin + index * 3;

			for (int dz = -1; dz < 2; dz++)
			for (int dy = -1; dy < 2; dy++)
			for (int dx = -1; dx < 2; dx++)
			{

				if (dx == 0 && dy == 0 && dz == 0) continue;

				int nx = x + dx;
				int ny = y + dy;
				int nz = z + dz;
				if (nx < 0 || nx >= ncols) continue;
				if (ny < 0 || ny >= nrows) continue;
				if (nz < 0 || nz >= nframes) continue;

				int indexnk = (nz*nrows + ny)*ncols + nx;

				uchar *ppixelnk = dataND+indexnk*nchannels;
				if (ppixelnk[0] == 255)
				{
					continue;
				}

				pXYZMaxnk[0] = (pXYZMax[0]>nx) ? pXYZMax[0] : nx;
				pXYZMaxnk[1] = (pXYZMax[1]>ny) ? pXYZMax[1] : ny;
				pXYZMaxnk[2] = (pXYZMax[2]>nz) ? pXYZMax[2] : nz;
				pXYZMinnk[0] = (pXYZMin[0]<nx) ? pXYZMin[0] : nx;
				pXYZMinnk[1] = (pXYZMin[1]<ny) ? pXYZMin[1] : ny;
				pXYZMinnk[2] = (pXYZMin[2]<nz) ? pXYZMin[2] : nz;

				float difmax = pXYZMaxnk[0] - pXYZMinnk[0];
				if (difmax < pXYZMaxnk[1] - pXYZMinnk[1])
					difmax = pXYZMaxnk[1] - pXYZMinnk[1];;

				difmax = lambdan * difmax / w;

				float difmaxz = pXYZMaxnk[2] - pXYZMinnk[2];
				difmaxz = lambdanz * difmaxz / wz;

				if (difmax < difmaxz)
					difmax = difmaxz;

				for (int b = 1; b < nchannels; b++)
				{
					pColorMaxnk[b] = (pMax[b]>ppixelnk[b]) ? pMax[b] : ppixelnk[b];
					pColorMinnk[b] = (pMin[b]<ppixelnk[b]) ? pMin[b] : ppixelnk[b];
					int dif = pColorMaxnk[b] - pColorMinnk[b];
					if (difmax < dif)
						difmax = dif;
				}

				int diffT = difmax;
				if (diffT >= maxcost) diffT = maxcost - 1;
				if (diffT < cost[indexnk])
				{
					cost[indexnk] = diffT;
					label[indexnk] = pixell;
					memcpy(pathColorMax + indexnk*nchannels, pColorMaxnk, nchannels);
					memcpy(pathColorMin + indexnk*nchannels, pColorMinnk, nchannels);

					pathXYZMax[indexnk * 3    ] = pXYZMaxnk[0];
					pathXYZMax[indexnk * 3 + 1] = pXYZMaxnk[1];
					pathXYZMax[indexnk * 3 + 2] = pXYZMaxnk[2];
					pathXYZMin[indexnk * 3    ] = pXYZMinnk[0];
					pathXYZMin[indexnk * 3 + 1] = pXYZMinnk[1];
					pathXYZMin[indexnk * 3 + 2] = pXYZMinnk[2];
					
					cluster[diffT].push_back(indexnk);
				}
			}
		}
	}

	delete pathColorMax;
	delete pathColorMin;
	delete pathXYZMax;
	delete pathXYZMin;
	delete pColorMaxnk;
	delete pColorMinnk;
	delete cost;
}

void RSS::PerformPRSSV_Range()
{
    uchar *pathColorMax = new uchar[npixels * nchannels];
    uchar *pathColorMin = new uchar[npixels * nchannels];

    int *pathXYZMax = new int[npixels * 3];
    int *pathXYZMin = new int[npixels * 3];

	float *cost = new float[npixels];

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
			int x = index % ncols;
			int y = index / ncols;
			int z = y / nrows;
			y = y % nrows;

			memcpy(pathColorMax+index*nchannels, dataND+index*nchannels, nchannels);
			memcpy(pathColorMin+index*nchannels, dataND+index*nchannels, nchannels);
			
			pathXYZMax[index * 3] = x;
			pathXYZMax[index * 3 + 1] = y;
			pathXYZMax[index * 3 + 2] = z;
			pathXYZMin[index * 3] = x;
			pathXYZMin[index * 3 + 1] = y;
			pathXYZMin[index * 3 + 2] = z;

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
				int x = index%ncols;
				int y = index/ncols;
				int z = y/nrows;
				y = y%nrows;

				uchar *ppixel = dataND+index*nchannels;
				if (ppixel[0] == 255)
				{
					continue;
				}

				ppixel[0] = 255;
				
				int pixell = label[index];

				uchar *pMax = pathColorMax+index*nchannels;
				uchar *pMin = pathColorMin+index*nchannels;

				int *pXYZMax = pathXYZMax + index * 3;
				int *pXYZMin = pathXYZMin + index * 3;

				uchar *pColorMaxnk = new uchar[nchannels];
				uchar *pColorMinnk = new uchar[nchannels];
				int pXYZMaxnk[3];
				int pXYZMinnk[3];

				for (int dz = -1; dz < 2; dz++)
				for (int dy = -1; dy < 2; dy++)
				for (int dx = -1; dx < 2; dx++)
				{

					if (dx == 0 && dy == 0 && dz == 0) continue;

					int nx = x + dx;
					int ny = y + dy;
					int nz = z + dz;
					if (nx < 0 || nx >= ncols) continue;
					if (ny < 0 || ny >= nrows) continue;
					if (nz < 0 || nz >= nframes) continue;

					int indexnk = (nz*nrows + ny)*ncols + nx;

					uchar *ppixelnk = dataND+indexnk*nchannels;
					if (ppixelnk[0] == 255)
					{
						continue;
					}

					pXYZMaxnk[0] = (pXYZMax[0]>nx) ? pXYZMax[0] : nx;
					pXYZMaxnk[1] = (pXYZMax[1]>ny) ? pXYZMax[1] : ny;
					pXYZMaxnk[2] = (pXYZMax[2]>nz) ? pXYZMax[2] : nz;
					pXYZMinnk[0] = (pXYZMin[0]<nx) ? pXYZMin[0] : nx;
					pXYZMinnk[1] = (pXYZMin[1]<ny) ? pXYZMin[1] : ny;
					pXYZMinnk[2] = (pXYZMin[2]<nz) ? pXYZMin[2] : nz;

					float difmax = pXYZMaxnk[0] - pXYZMinnk[0];
					if (difmax < pXYZMaxnk[1] - pXYZMinnk[1])
						difmax = pXYZMaxnk[1] - pXYZMinnk[1];;

					difmax = lambdan * difmax / w;

					float difmaxz = pXYZMaxnk[2] - pXYZMinnk[2];
					difmaxz = lambdanz * difmaxz / wz;

					if (difmax < difmaxz)
						difmax = difmaxz;

					for (int b = 1; b < nchannels; b++)
					{
						pColorMaxnk[b] = (pMax[b]>ppixelnk[b]) ? pMax[b] : ppixelnk[b];
						pColorMinnk[b] = (pMin[b]<ppixelnk[b]) ? pMin[b] : ppixelnk[b];
						int dif = pColorMaxnk[b] - pColorMinnk[b];
						if (difmax < dif)
							difmax = dif;
					}

					int diffT = difmax;
					if (diffT >= maxcost) diffT = maxcost - 1;
					if (diffT < cost[indexnk])
					{
						cost[indexnk] = diffT;
						label[indexnk] = pixell;
						memcpy(pathColorMax + indexnk*nchannels, pColorMaxnk, nchannels);
						memcpy(pathColorMin + indexnk*nchannels, pColorMinnk, nchannels);

						pathXYZMax[indexnk * 3    ] = pXYZMaxnk[0];
						pathXYZMax[indexnk * 3 + 1] = pXYZMaxnk[1];
						pathXYZMax[indexnk * 3 + 2] = pXYZMaxnk[2];
						pathXYZMin[indexnk * 3    ] = pXYZMinnk[0];
						pathXYZMin[indexnk * 3 + 1] = pXYZMinnk[1];
						pathXYZMin[indexnk * 3 + 2] = pXYZMinnk[2];
						
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
	delete pathXYZMax;
	delete pathXYZMin;
	delete cost;
}