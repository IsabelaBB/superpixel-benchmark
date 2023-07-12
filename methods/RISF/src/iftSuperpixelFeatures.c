#include "iftSuperpixelFeatures.h"

// ---------- Private Methods Declaration ----------

/**
 * @brief Computes histogram assignment of each pixel.
 * @author Felipe Lemes Galvao
 * @date January 5, 2017
 */
iftImage * iftColorHistogramAssignPixelBin(const iftImage *img, int binsPerColor);

// ---------- Private Methods Implementation ----------

iftImage * iftColorHistogramAssignPixelBin(const iftImage *img, int binsPerColor)
{
  iftMImage *mimg = NULL;
  if (iftIsColorImage(img))
    mimg = iftImageToMImage(img, LABNorm_CSPACE);
  else
    mimg = iftImageToMImage(img, GRAYNorm_CSPACE);
  for (int b = 0; b < mimg->m; ++b) {
    float minB = mimg->val[0][b];
    for (int p = 0; p < mimg->n; ++p) {
      float val = mimg->val[p][b];
      if (val < minB)
        minB = val;
    }
    for (int p = 0; p < mimg->n; ++p) {
      mimg->val[p][b] += minB;
      if (mimg->val[p][b] < 0.0f)
        mimg->val[p][b] = 0.0f;
    }
    float max=IFT_INFINITY_FLT_NEG;
    for (int p=0; p < mimg->n; p++)
      if (mimg->val[p][b]>max)
	max = mimg->val[p][b];
    if (max > IFT_EPSILON)
      for (int p=0; p < mimg->n; p++)
	mimg->val[p][b] /= max;
  }

  iftImage *res = iftCreateImage(img->xsize, img->ysize, img->zsize);
  for (int p = 0; p < res->n; ++p) {
    int bin = 0;
    int mult = 1;
    for (int b = 0; b < mimg->m; ++b) {
      int dimBin = floor(mimg->val[p][b] * binsPerColor);
      /*if (dimBin < 0) {
        printf("dimBin = floor (%f * %d)\n", mimg->val[p][b], binsPerColor);
      }*/
      if (dimBin >= binsPerColor) // fix when val = 1.0
        dimBin = binsPerColor - 1;
      assert(dimBin >= 0);
      bin += dimBin * mult;
      mult *= binsPerColor;
    }
    res->val[p] = bin;
  }

  iftDestroyMImage(&mimg);
  return res;
}

// ---------- Public Methods Implementation ----------

iftMatrix * iftComputeSuperpixelFeaturesByColorSpaceMedoid(const iftImage *superpixelLabelMap, const iftImage *img, iftColorSpace colorSpace)
{
  iftMImage *mimg = iftImageToMImage(img, colorSpace);
  iftMatrix *res = iftComputeSuperpixelFeaturesByMImageMedoid(mimg, superpixelLabelMap);
  iftDestroyMImage(&mimg);
  return res;
}

iftMatrix * iftComputeSuperpixelFeaturesByMImageMedoid(const iftMImage *mimg, const iftImage *superpixelLabelMap)
{
  assert(mimg != NULL && superpixelLabelMap != NULL);
  assert(mimg->xsize == superpixelLabelMap->xsize);
  assert(mimg->ysize == superpixelLabelMap->ysize);
  assert(mimg->zsize == superpixelLabelMap->zsize);
  // TODO assert(check label map consistency) 

  // Piggyback from the mean calculation
  iftMatrix *mx = iftComputeSuperpixelFeaturesByMImageMean(mimg, superpixelLabelMap);

  // Keep track of currently known closest pixel for each superpixel
  int *medoidCandidate = iftAllocIntArray(mx->nrows);
  float *candidateDist = iftAllocFloatArray(mx->nrows);
  for (int i = 0; i < mx->nrows; ++i) {
    medoidCandidate[i] = -1;
    candidateDist[i] = FLT_MAX;
  }

  // Try to update closest pixel for each pixel
  for (int p = 0; p < mimg->n; ++p) {
    int row = superpixelLabelMap->val[p] - 1;
    if (row < 0)
      continue;
    // Compute distance from its superpixel mean
    float dist = 0.0;
    for (int f = 0; f < mimg->m; ++f) {
      float a = mimg->val[p][f];
      float b = mx->val[iftGetMatrixIndex(mx,f,row)];
      dist += (a-b)*(a-b);
    }
    // Update on improvement
    if (dist < candidateDist[row]) {
      medoidCandidate[row] = p;
      candidateDist[row] = dist;
    }
  }

  // Overwrite mean matrix with medoids to avoid new allocation
  for (int row = 0; row < mx->nrows; ++row) {
    for (int f = 0; f < mimg->m; ++f) {
      float feat = mimg->val[medoidCandidate[row]][f];
      mx->val[iftGetMatrixIndex(mx,f,row)] = feat;
    }
  }

  free(medoidCandidate);
  free(candidateDist);
  return(mx);
}

iftMatrix * iftComputeSuperpixelFeaturesByMImageMean(const iftMImage *mimg, const iftImage *superpixelLabelMap)
{
  assert(mimg != NULL && superpixelLabelMap != NULL);
  assert(mimg->xsize == superpixelLabelMap->xsize);
  assert(mimg->ysize == superpixelLabelMap->ysize);
  assert(mimg->zsize == superpixelLabelMap->zsize);
  // TODO assert(check label map consistency) 

  int nSuperpixels = iftMaximumValue(superpixelLabelMap);
  iftMatrix *mx = iftCreateMatrix(mimg->m, nSuperpixels);
  int *pixelCount = iftAllocIntArray(nSuperpixels);

  /* Sum pixel features for each superpixel with 
   *   naive float accumulation, O(n) error. */
  for (int p = 0; p < mimg->n; ++p) {
    if (superpixelLabelMap->val[p] == 0)
      continue;
    int row = superpixelLabelMap->val[p] - 1;
    for (int f = 0; f < mimg->m; ++f)
      mx->val[iftGetMatrixIndex(mx,f,row)] += mimg->val[p][f];
    pixelCount[row] += 1;
  }

  for (int row = 0; row < nSuperpixels; ++row)
    for (int f = 0; f < mimg->m; ++f)
      mx->val[iftGetMatrixIndex(mx,f,row)] /= pixelCount[row];

  free(pixelCount);
  return mx;
}

iftMatrix * iftComputeSuperpixelFeaturesByColorHistogram(const iftImage *img, const iftImage *labelMap, int binsPerColor)
{
  int nSuperpixels = iftMaximumValue(labelMap);
  int nDims = iftIsColorImage(img) ? 3 : 1;
  int nBins = lround(pow(binsPerColor, nDims));
  iftMatrix *feats = iftCreateMatrix(nBins, nSuperpixels);
  iftImage *binImage = iftColorHistogramAssignPixelBin(img, binsPerColor);

  for (int p = 0; p < labelMap->n; ++p) {
    int label = labelMap->val[p] - 1;
    iftMatrixRowPointer(feats, label)[binImage->val[p]] += 1.0;
  }

  iftDestroyImage(&binImage);
  return feats;
}

