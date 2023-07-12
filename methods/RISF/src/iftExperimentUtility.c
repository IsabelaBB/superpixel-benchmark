#include "iftExperimentUtility.h"

iftMatrix *iftCountLabelsInSuperpixels(const iftImage *spLabelMap, const iftImage *baseLabelMap, bool isNormalized)
{
  assert(spLabelMap != NULL);
  assert(baseLabelMap != NULL);
  assert(iftIsDomainEqual(spLabelMap, baseLabelMap));

  // Extract implicit info from label maps
  int nSuperpixels = iftMaximumValue(spLabelMap);
  int nLabels = iftMaximumValue(baseLabelMap) + 1; // + background 0
  iftMatrix *labelCount = iftCreateMatrix(nLabels, nSuperpixels);

  // Count the labels in each superpixel
  for (int p = 0; p < spLabelMap->n; ++p) {
    // Skip pixel with no associated superpixel
    if (baseLabelMap->val[p] == 0)
      continue;

    int superpixel = spLabelMap->val[p] - 1;
    int label = baseLabelMap->val[p];
    int index = iftGetMatrixIndex(labelCount, label, superpixel);
    labelCount->val[index] += 1.0;
  }

  if (isNormalized) {
    for (int row = 0; row < labelCount->nrows; ++row) {
      float sum = 0.0;
      for (int col = 0; col < labelCount->ncols; ++col) {
        int index = iftGetMatrixIndex(labelCount, col, row);
        sum += labelCount->val[index];
      }
      for (int col = 0; col < labelCount->ncols; ++col) {
        int index = iftGetMatrixIndex(labelCount, col, row);
        labelCount->val[index] /= sum;
      }
    }
  }

  return labelCount;
}

// ---------- Public Methods Implementation ----------

iftImage *iftOverlaySegmentationBorders(iftImage *img, iftImage *segLabelMap, iftColor YCbCr)
{
  assert(img != NULL);
  assert(segLabelMap != NULL);
  assert(iftIsDomainEqual(img, segLabelMap));

  iftImage *res = iftCopyImage(img);
  iftAdjRel *A = iftCircular(1.0);
  iftAdjRel *B = iftCircular(0.0);
  iftDrawBorders(res, segLabelMap, A, YCbCr, B);

  iftDestroyAdjRel(&A);
  iftDestroyAdjRel(&B);
  return res;
}

iftImage *iftSuperpixelToMajoritySegmentation(const iftImage *segLabelMap, const iftImage *gtLabelMap)
{
  assert(segLabelMap != NULL);
  assert(gtLabelMap != NULL);
  assert(iftIsDomainEqual(segLabelMap, gtLabelMap));

  // Copy only to get domain, we will overwrite all values
  iftImage *res = iftCopyImage(segLabelMap);

  // Extract implicit info from label maps. 
  iftMatrix *labelCount = iftCountLabelsInSuperpixels(segLabelMap, gtLabelMap, false);

  // Find the mode of each superpixel
  int *superpixelMode = iftAllocIntArray(labelCount->nrows); 
  for (int i = 0; i < labelCount->nrows; ++i) {
    float *superpixelRow = iftMatrixRowPointer(labelCount, i);
    superpixelMode[i] = iftFArgmax(superpixelRow, labelCount->ncols);
  }

  // Map each superpixel to its mode
  for (int p = 0; p < res->n; ++p) {
    // Pixels without superpixels are assigned background by default
    if (segLabelMap->val[p] == 0) {
      res->val[p] = 0;
      continue;
    }

    int superpixel = segLabelMap->val[p] - 1;
    res->val[p] = superpixelMode[superpixel];
  }

  iftDestroyMatrix(&labelCount);
  free(superpixelMode);
  return res;
}

float iftAchievableSegmentationAccuracy(const iftImage *segLabelMap, const iftImage *gtLabelMap)
{
  assert(segLabelMap != NULL);
  assert(gtLabelMap != NULL);
  assert(iftIsDomainEqual(segLabelMap, gtLabelMap));

  iftImage *bestSeg = iftSuperpixelToMajoritySegmentation(segLabelMap, gtLabelMap);
  int correctPixels = 0;
  for (int p = 0; p < gtLabelMap->n; ++p)
    if (bestSeg->val[p] == gtLabelMap->val[p])
      correctPixels += 1;

  float ASA = (float) correctPixels / (float) gtLabelMap->n;
  iftDestroyImage(&bestSeg);
  
  return ASA;
}

float iftBoundaryPrecision(const iftImage *segBorders, const iftImage *gtBorders, float toleranceDist) 
{
  assert(segBorders != NULL);
  assert(gtBorders != NULL);
  assert(iftIsDomainEqual(segBorders, gtBorders));
  assert(toleranceDist >= 0.0);

  iftAdjRel *A = NULL;
  if (iftIs3DImage(segBorders))
    A = iftSpheric(toleranceDist);
  else
    A = iftCircular(toleranceDist);

  int truePositives = 0;
  int nPositives = 0;
  for (int p = 0; p < segBorders->n; ++p) {
    if (segBorders->val[p] == 0)
      continue;

    nPositives += 1;

    iftVoxel u = iftGetVoxelCoord(segBorders, p);
    for (int i = 0; i < A->n; ++i) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, i);
      if (!iftValidVoxel(segBorders, v))
        continue;

      int q = iftGetVoxelIndex(segBorders, v);
      if (gtBorders->val[q] != 0) {
        truePositives += 1;
        break;
      }
    }
  }

  assert(nPositives > 0.0);
  return (float)truePositives / (float)nPositives;
}

float iftBoundaryFScore(const iftImage *segBorders, const iftImage *gtBorders, float toleranceDist)
{
  assert(segBorders != NULL);
  assert(gtBorders != NULL);
  assert(iftIsDomainEqual(segBorders, gtBorders));
  assert(toleranceDist >= 0.0);

  float recall = iftBoundaryRecall((iftImage *) gtBorders,
      (iftImage *) segBorders, toleranceDist);
  float precision = iftBoundaryPrecision(segBorders, gtBorders, toleranceDist);

  if (recall * precision < IFT_EPSILON)
    return 0.0f;

  return (2.0f * precision * recall) / (precision + recall);
}

iftImage *iftForceLabelMapConnectivity(iftImage *labelMap, int minSize)
{
  // Init
  iftImage *res = iftCreateImage(labelMap->xsize, labelMap->ysize, labelMap->zsize);
  iftAdjRel *A = iftIs3DImage(labelMap) ? iftSpheric(1.0f) : iftCircular(1.0f);
  iftFIFO *Q = iftCreateFIFO(res->n);

  // Each connected component is given a unique label if their size > minSize
  int label = 1;
  for (int pTop = 0; pTop < res->n; ++pTop) {
    // Connected component already processed
    if (res->val[pTop] > 0)
      continue;

    // Initialization
    iftResetFIFO(Q);
    iftInsertFIFO(Q, pTop);
    int size = 0;
    int neighLabel = 0;
    iftSet *componentVoxels = NULL;

    // Cover entire component with BFS
    res->val[pTop] = label;
    while (!iftEmptyFIFO(Q)) {
      // Process voxel
      int p = iftRemoveFIFO(Q);
      iftInsertSet(&componentVoxels, p);
      size += 1;

      // Visit neighbors
      iftVoxel u = iftGetVoxelCoord(res, p);
      for (int i = 1; i < A->n; ++i) {
        iftVoxel v = iftGetAdjacentVoxel(A, u, i);
        if (!iftValidVoxel(res, v))
          continue;
        int q = iftGetVoxelIndex(res, v);
        if (res->val[q] == 0) {
          if (labelMap->val[p] == labelMap->val[q]) {
            res->val[q] = label;
            iftInsertFIFO(Q, q);
          }
        } else if (res->val[q] != label) {
          neighLabel = res->val[q];
        }
      }
    }

    // Relabel with neighbor if too small (except first one)
    if (size < minSize && neighLabel > 0) {
      while (componentVoxels != NULL) {
        int p = iftRemoveSet(&componentVoxels);
        res->val[p] = neighLabel;
      }
    } else {
      label += 1;
      iftDestroySet(&componentVoxels);
    }
  }

  iftDestroyAdjRel(&A);
  iftDestroyFIFO(&Q);

  return res;
}

extern bool _iftHasCSVHeader(const char *csv_pathname, char separator);
extern void _iftCountNumOfRowsAndColsFromCSVFile(const char *csv_pathname, long *nrows, long *ncols, char separator);
extern iftCSV *_iftCreateCSVWithoutStringAllocation(long nrows, long ncols);
iftCSV *iftReadCSVWithHeader(const char *csv_pathname, const char separator, bool *has_header) 
{
    if (!iftFileExists(csv_pathname))
        iftError("The CSV file pathname \"%s\" does not exist!", "iftReadCSV", csv_pathname);
    
    char strSeparator[2] = {separator, '\0'};
    
    *has_header = _iftHasCSVHeader(csv_pathname, separator);
    
    long nrows, ncols;
    _iftCountNumOfRowsAndColsFromCSVFile(csv_pathname, &nrows, &ncols, separator);
    
    iftCSV *csv = _iftCreateCSVWithoutStringAllocation(nrows, ncols);
    
    FILE *fp = fopen(csv_pathname, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "_iftCountNumOfRowsAndColsFromCSVFile", csv_pathname);
    
    // copies the values from the CSV file
    iftSList *SL = NULL;
    char *line = iftGetLine(fp);
    
    long i = 0;
    while (line != NULL) {
        SL = iftSplitString(line, strSeparator);
        
        for (long j = 0; j < csv->ncols; j++) {
            csv->data[i][j] = iftRemoveSListHead(SL); // just points to string
            // removes the '\n' and '\r' from the paths
            iftRightTrim(csv->data[i][j], '\n');
            iftRightTrim(csv->data[i][j], '\r');
        }
        i++;
        iftDestroySList(&SL);
        
        iftFree(line);
        line = iftGetLine(fp);
    }
    
    fclose(fp);
    
    return csv;
}

