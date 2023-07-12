/**
 * @file
 * @brief Computes metrics of an existing segmentation.
 * @author Felipe Lemes Galvao
 * @date December, 2017
 */

#include <ift.h>
#include <iftExperimentUtility.h>

void fixDoubleSlashes(char *str)
{
  int len = strlen(str);
  for (int i = 1; i < len; ++i) {
    if ((str[i] == '/' || str[i] == '\\') && str[i] == str[i-1]) {
      memmove(&(str[i]), &(str[i+1]), strlen(&(str[i])));
      len -= 1;
      i -= 1;
    }
  }
}

void stripPathToSuffix(char *str)
{
  // Find last /
  int i = 0;
  int lastSlash = 0;
  while (str[i] != '\0') {
    if (str[i] == '/')
      lastSlash = i;
    i++;
  }
  strcpy(str, str + lastSlash + 1);
}

void loadSegmentation(const char *path, iftImage **labelMap, iftImage **borders, bool isBorderSeg)
{
  if (isBorderSeg) {
    *borders =  iftReadImageByExt(path);
    // Temporary workaround because iftBorderImageToLabelImage
    //   changes the input value inside the iftWaterGray function.
    iftImage *tmp = iftCopyImage(*borders);
    *labelMap = iftBorderImageToLabelImage(tmp);
    iftDestroyImage(&tmp);
  } else {
    iftImage *tmp = iftReadImageByExt(path);
    *labelMap = iftForceLabelMapConnectivity(tmp, 0);
    int labelDiff = iftMaximumValue(*labelMap) - iftMaximumValue(tmp);
    if (labelDiff > 0) {
      fprintf(stderr, "Warning: supplied label map is not 4-connected.\n");
      fprintf(stderr, "Disconnected components have been relabeled (%d new labels)\n", labelDiff);
    }
    iftDestroyImage(&tmp);
    *borders = iftBorderImage(*labelMap, 0);
  }
}

float calcBinomial(float val)
{
  return ((val * (val - 1.0f))/2.0f);
}

/*
 * Reference:
 *   "A Measure for Objective Evaluation of Image Segmentation Algorithms", Unnikrishnan et al. (2005)
 *   Section 3.2
 */
float iftSegmentationRandIndex(iftImage *seg, iftImage *gt)
{
  assert(seg != NULL);
  assert(gt != NULL);
  assert(seg->xsize == gt->xsize);
  assert(seg->ysize == gt->ysize);
  assert(seg->zsize == gt->zsize);

  int nLabelsSeg = iftMaximumValue(seg);
  int nLabelsGt = iftMaximumValue(gt);
  iftMatrix *contingencyTable = iftCreateMatrix(nLabelsSeg, nLabelsGt); // ncols, nrows

  // Fill contingency table
  for (int p = 0; p < seg->n; ++p) {
    int col = seg->val[p] - 1;
    int row = gt->val[p] - 1;
    iftMatrixElem(contingencyTable, col, row) += 1.0f;
  }

  float nAgreeOnBoth = 0;
  for (int i = 0; i < contingencyTable->n; ++i)
    nAgreeOnBoth += calcBinomial(contingencyTable->val[i]);

  // Agree on seg and disagree on gt
  float nAgreeSeg = -nAgreeOnBoth;
  for (int segLabel = 0; segLabel < contingencyTable->nrows; ++segLabel) {
    float sumWithSegLabel = 0.0f;
    for (int col = 0; col < contingencyTable->ncols; ++col)
      sumWithSegLabel += iftMatrixElem(contingencyTable, col, segLabel);
    nAgreeSeg += calcBinomial(sumWithSegLabel);
  } 

  // Agree on gt and disagree on seg
  float nAgreeGt = -nAgreeOnBoth;
  for (int gtLabel = 0; gtLabel < contingencyTable->ncols; ++gtLabel) {
    float sumWithGtLabel = 0.0f;
    for (int row = 0; row < contingencyTable->nrows; ++row)
      sumWithGtLabel += iftMatrixElem(contingencyTable, gtLabel, row);
    nAgreeGt += calcBinomial(sumWithGtLabel);
  } 

  // Disagree on both
  float nPairs = calcBinomial((float)seg->n);
  float nDisagreeOnBoth = nPairs - nAgreeOnBoth - nAgreeSeg - nAgreeGt;

  float randIndex = (nAgreeOnBoth + nDisagreeOnBoth) / nPairs;

  iftDestroyMatrix(&contingencyTable);
  return randIndex;
}

/*
 * We convention on using segmentation as base and ground-truth as cover,
 *   but both directions are valid metrics.
 *
 * Reference:
 *  "Contour Detection and Hierarchical Image Segmentation", Arbelaez et al. (2011)
 *  - Section 2.3.3
 */
float iftSegmentationCoveringMetric(iftImage *base, iftImage *cover)
{
  assert(base != NULL);
  assert(cover != NULL);
  assert(base->xsize == cover->xsize);
  assert(base->ysize == cover->ysize);
  assert(base->zsize == cover->zsize);

  // -- Compute all pairs of overlap between regions
  int nLabelsBase = iftMaximumValue(base);
  int nLabelsCover = iftMaximumValue(cover);
  iftMatrix *andMx = iftCreateMatrix(nLabelsBase, nLabelsCover); // ncols, nrows
  iftMatrix *orMx = iftCreateMatrix(nLabelsBase, nLabelsCover);
  iftMatrix *overlapMx = iftCreateMatrix(nLabelsBase, nLabelsCover);

  // Compute intersection and union terms for each pair separately
  for (int p = 0; p < base->n; ++p) {
    int baseLabel = base->val[p] - 1;
    int coverLabel = cover->val[p] - 1;

    // Intersection
    iftMatrixElem(andMx, baseLabel, coverLabel) += 1.0f;

    // Union
    for (int row = 0; row < orMx->nrows; ++row)
      iftMatrixElem(orMx, baseLabel, row) += 1.0f;
    for (int col = 0; col < orMx->ncols; ++col)
      iftMatrixElem(orMx, col, coverLabel) += 1.0f;
    iftMatrixElem(orMx, baseLabel, coverLabel) -= 1.0f;
  }

  // Compute final overlap
  for (int i = 0; i < andMx->n; ++i)
    overlapMx->val[i] = andMx->val[i] / orMx->val[i];

  // -- Compute covering based on overlap
  int *baseRegionSize = iftAllocIntArray(nLabelsBase);
  for (int i = 0; i < base->n; ++i) {
    int label = base->val[i] - 1;
    baseRegionSize[label] += 1;
  }

  // Covering computation (Eq. 8 on ref)
  float covering = 0.0f;
  for (int baseLabel = 0; baseLabel < nLabelsBase; ++baseLabel) {
    float maxOverlap = 0.0f;
    for (int row = 0; row < overlapMx->nrows; ++row) {
      float pairOverlap = iftMatrixElem(overlapMx, baseLabel, row);
      if (pairOverlap > maxOverlap)
	maxOverlap = pairOverlap;
    }
    covering += maxOverlap * baseRegionSize[baseLabel];
  }
  covering /= ((float) base->n);

  iftDestroyMatrix(&andMx);
  iftDestroyMatrix(&orMx);
  iftDestroyMatrix(&overlapMx);
  free(baseRegionSize);

  return covering;
}

/*
 * Reference:
 *  "Comparing Clusterings - An Axiomatic View", Meila (2005)
 *  - Section 3 (notation in Section 2)
 *
 * Mutual information (MI) is defined incorrectly there 
 *   (marginal probabilities should be inverted).
 * Either way we use alternative definition for MI as H(X) + H(Y) - H(X,Y).
 */
float iftSegmentationVariationOfInformation(iftImage *seg, iftImage *gt)
{
  assert(seg != NULL);
  assert(gt != NULL);
  assert(seg->xsize == gt->xsize);
  assert(seg->ysize == gt->ysize);
  assert(seg->zsize == gt->zsize);
  

  // Get cluster and cluster intersection sizes
  int nLabelsSeg = iftMaximumValue(seg);
  int nLabelsGt = iftMaximumValue(gt);
  float *segRegionSize = iftAllocFloatArray(nLabelsSeg);
  float *gtRegionSize = iftAllocFloatArray(nLabelsGt);
  iftMatrix *jointSize = iftCreateMatrix(nLabelsSeg, nLabelsGt);

  for (int p = 0; p < seg->n; ++p) {
    int segLabel = seg->val[p] - 1;
    int gtLabel = gt->val[p] - 1;

    segRegionSize[segLabel] += 1.0f;
    gtRegionSize[gtLabel] += 1.0f;
    iftMatrixElem(jointSize, segLabel, gtLabel) += 1.0f;
  }

  // Entropy
  float segEntropy = 0.0f;
  for (int segLabel = 0; segLabel < nLabelsSeg; ++segLabel) {
    float frac = segRegionSize[segLabel] / ((float) seg->n);
    segEntropy -= (frac * logf(frac));
  }
  float gtEntropy = 0.0f;
  for (int gtLabel = 0; gtLabel < nLabelsGt; ++gtLabel) {
    float frac = gtRegionSize[gtLabel] / ((float) gt->n);
    gtEntropy -= (frac * logf(frac));
  }


  // Mutual information
  float jointEntropy = 0.0f;
  for (int segLabel = 0; segLabel < nLabelsSeg; ++segLabel) {
    for (int gtLabel = 0; gtLabel < nLabelsGt; ++gtLabel) {
      float frac = iftMatrixElem(jointSize, segLabel, gtLabel) / ((float) seg->n);
      if (frac > IFT_EPSILON)
        jointEntropy -= (frac * logf(frac));
    }
  }
  float mutualInformation = segEntropy + gtEntropy - jointEntropy;

  float variationOfInformation = segEntropy + gtEntropy - 2 * mutualInformation;

  iftDestroyMatrix(&jointSize);
  free(segRegionSize);
  free(gtRegionSize);

  return variationOfInformation;
}

#define SEG_ARG "--input-segmentation"
#define GT_ARG "--input-groundtruth"
#define ISBORDERSEG_ARG "--is-border-seg"
#define ISBORDERGT_ARG "--is-border-gt"
#define PRINT_ARG "--print-opt"
#define UNION_ARG "--border-union"
#define TOLERANCE_ARG "--boundary-tolerance"
#define BR_ARG "--boundary-recall"
#define BP_ARG "--boundary-precision"
#define BF_ARG "--boundary-fscore"
#define UE_ARG "--undersegmentation"
#define ASA_ARG "--achievable-seg"
#define COMP_ARG "--compactness"
#define SPMETRICS_ARG "--std-superpixels"
#define DICE_ARG "--dice"
#define TOP_ARG "--topology"
#define COVER_ARG "--covering"
#define VOI_ARG "--variation-of-info"
#define PRI_ARG "--prob-rand-index" 
#define SEGMETRICS_ARG "--large-segment-metrics"

iftDict *iftGetArguments(int argc, const char *argv[]);

int main(int argc, char* argv[])
{
  iftDict* args = iftGetArguments(argc, (const char **) argv);
  char *tmpPath;

  bool isBorderGt = iftDictContainKey(ISBORDERGT_ARG, args, NULL);
  bool isBorderSeg = iftDictContainKey(ISBORDERSEG_ARG, args, NULL);
  int printOpt = (iftDictContainKey(PRINT_ARG, args, NULL)) ?
    iftGetLongValFromDict(PRINT_ARG, args) : 1;
  bool doBF = iftDictContainKey(BF_ARG, args, NULL) || iftDictContainKey(SEGMETRICS_ARG, args, NULL);
  bool doBR = iftDictContainKey(BR_ARG, args, NULL) || iftDictContainKey(SPMETRICS_ARG, args, NULL) || doBF;
  bool doBP = iftDictContainKey(BP_ARG, args, NULL) || doBF;
  bool doUnion = iftDictContainKey(UNION_ARG, args, NULL) && (doBR || doBP || doBF);
  bool doUE = iftDictContainKey(UE_ARG, args, NULL) || iftDictContainKey(SPMETRICS_ARG, args, NULL);
  bool doASA = iftDictContainKey(ASA_ARG, args, NULL) || iftDictContainKey(SPMETRICS_ARG, args, NULL);
  bool doComp = iftDictContainKey(COMP_ARG, args, NULL) || iftDictContainKey(SPMETRICS_ARG, args, NULL);
  bool doDice = iftDictContainKey(DICE_ARG, args, NULL);
  int nObjects = doDice ? iftGetLongValFromDict(DICE_ARG, args) : 0;
  bool doTop = iftDictContainKey(TOP_ARG, args, NULL);
  bool doCover = iftDictContainKey(COVER_ARG, args, NULL) || iftDictContainKey(SEGMETRICS_ARG, args, NULL);
  bool doVOI = iftDictContainKey(VOI_ARG, args, NULL) || iftDictContainKey(SEGMETRICS_ARG, args, NULL);
  bool doPRI = iftDictContainKey(PRI_ARG, args, NULL) || iftDictContainKey(SEGMETRICS_ARG, args, NULL);

  // Load groundtruth paths based on prefix
  tmpPath = iftGetStrValFromDict(GT_ARG, args);
  fixDoubleSlashes(tmpPath); // Regex won't work otherwise
  char gtRegex[256];
  char *folderPath = iftParentDir(tmpPath);
  stripPathToSuffix(tmpPath);
  sprintf(gtRegex, "^%s.*", tmpPath);
  iftDir *gtDir = iftLoadFilesFromDirByRegex(folderPath, gtRegex); 
  if (gtDir->nfiles <= 0 && printOpt < 3)
    iftError("Invalid groundtruth prefix.", "main");
  bool isMultiGt = gtDir->nfiles > 1;
  free(tmpPath);
  free(folderPath);

  // Only print header for CSV
  if (printOpt == 3 || printOpt == 4) {
    char stdStr[5] = "";
    if (printOpt == 4)
      strcpy(stdStr, ",std");
    printf("nSp%s", stdStr);
    if (doBR)
      printf(",BR%s", stdStr);
    if (isMultiGt && doUnion && doBR)
      printf(",UBR%s", stdStr);
    if (doBP)
      printf(",BP%s", stdStr);
    if (isMultiGt && doUnion && doBP)
      printf(",UBP%s", stdStr);
    if (doBF)
      printf(",BF%s", stdStr);
    if (isMultiGt && doUnion && doBF)
      printf(",UBF%s", stdStr);
    if (doASA)
      printf(",ASA%s", stdStr);
    if (doUE)
      printf(",UE%s", stdStr);
    if (doComp)
      printf(",Comp%s", stdStr);
    if (doTop)
      printf(",Top%s", stdStr);
    if (doCover)
      printf(",Covering%s", stdStr);
    if (doVOI)
      printf(",VOI%s", stdStr);
    if (doPRI)
      printf(",PRI%s", stdStr);
    if (doDice) {
      printf(",DiceAvg%s", stdStr);
      for (int i = 0; i < nObjects; ++i)
        printf(",DiceObj%d%s", i + 1, stdStr);
    }

    iftDestroyDir(&gtDir);
    return 0;
  }

  // Load existing segmentation
  iftImage *segLabelMap = NULL;
  iftImage *segBorders = NULL;
  tmpPath = iftGetStrValFromDict(SEG_ARG, args);
  loadSegmentation(tmpPath, &segLabelMap, &segBorders, isBorderSeg);
  bool is3D = iftIs3DImage(segLabelMap);
  free(tmpPath);

  // Compute groundtruth based metrics
  const float boundaryTolerance = (iftDictContainKey(TOLERANCE_ARG, args, NULL)) ?
	  iftGetDblValFromDict(TOLERANCE_ARG, args) : 2.0;
  float boundaryRecall = 0.0;
  float boundaryPrecision = 0.0;
  float boundaryFScore = 0.0;
  float ASA = 0.0;
  float underSegmentationError = 0.0;
  float covering = 0.0;
  float VOI = 0.0;
  float PRI = 0.0;
  float *dice = NULL;
  if (doDice) {
    dice = iftAllocFloatArray(nObjects + 1); 
  }
  iftImage *gtBordersUnion = NULL;
  if (isMultiGt && doUnion)
    gtBordersUnion = iftCreateImage(segLabelMap->xsize,
        segLabelMap->ysize, segLabelMap->zsize);
  for (int g = 0; g < gtDir->nfiles; ++g) {
    iftImage *gtLabelMap = NULL;
    iftImage *gtBorders = NULL;
    loadSegmentation(gtDir->files[g]->path, &gtLabelMap,
        &gtBorders, isBorderGt);
    
    if (isMultiGt && doUnion) {
      for (int p = 0; p < gtBorders->n; ++p)
        if (gtBorders->val[p] != 0)
          gtBordersUnion->val[p] = gtBorders->val[p];
    }

    // Segmentation metrics
    if (doBR)
      boundaryRecall += 
        iftBoundaryRecall(gtBorders, segBorders, boundaryTolerance);
    if (doBP)
      boundaryPrecision +=
        iftBoundaryPrecision(segBorders, gtBorders, boundaryTolerance);
    if (doBF)
      boundaryFScore +=
        iftBoundaryFScore(segBorders, gtBorders, boundaryTolerance);
    if (doASA)
      ASA += iftAchievableSegmentationAccuracy(segLabelMap, gtLabelMap);
    if (doUE)
      underSegmentationError += 
        iftUnderSegmentation(gtLabelMap, segLabelMap);
    if (doDice) {
      iftImage *bestSegmentation = iftSuperpixelToMajoritySegmentation(segLabelMap, gtLabelMap);
      float *diceAux = iftFScoreMultiLabel(bestSegmentation, gtLabelMap, iftMaximumValue(gtLabelMap));
      int nObjectsUpdate = iftMin(nObjects, iftMaximumValue(gtLabelMap)) + 1;
      for (int i = 0; i < nObjectsUpdate; ++i)
        dice[i] += diceAux[i];
      iftDestroyImage(&bestSegmentation);
      free(diceAux);
    }
    if (doCover)
      covering += iftSegmentationCoveringMetric(segLabelMap, gtLabelMap);
    if (doVOI)
      VOI += iftSegmentationVariationOfInformation(segLabelMap, gtLabelMap);
    if (doPRI)
      PRI += iftSegmentationRandIndex(segLabelMap, gtLabelMap);

    iftDestroyImage(&gtLabelMap);
    iftDestroyImage(&gtBorders);
  }
  boundaryRecall /= gtDir->nfiles;
  boundaryPrecision /= gtDir->nfiles;
  boundaryFScore /= gtDir->nfiles;
  ASA /= gtDir->nfiles;
  underSegmentationError /= gtDir->nfiles;
  covering /= gtDir->nfiles;
  VOI /= gtDir->nfiles;
  PRI /= gtDir->nfiles;
  if (doDice)
    for (int i = 0; i < nObjects + 1; ++i)
      dice[i] /= gtDir->nfiles;

  // Compute border metrics with combination of borders for multiple gts
  float brUnion, bpUnion, bfUnion;
  if (isMultiGt && doUnion) {
    if (doBR)
      brUnion = iftBoundaryRecall(gtBordersUnion,
          segBorders, boundaryTolerance);
    if (doBP)
      bpUnion = iftBoundaryPrecision(segBorders,
          gtBordersUnion, boundaryTolerance);
    if (doBF)
      bfUnion = iftBoundaryFScore(segBorders,
          gtBordersUnion, boundaryTolerance);
  } else {
    brUnion = bpUnion = bfUnion = 0.0;
  }

  // Superpixel properties
  int nSp = iftMaximumValue(segLabelMap);
  float comp = (!is3D && !doComp) ? 0.0 : iftCompactness2D(segLabelMap);
  // TopologyMeasure is exceeding memory consumption on large images
  float top = doTop ? iftTopologyMeasure(segLabelMap) : 0.0;

  // Print results
  switch (printOpt) {
    case 1: // Human readable stats
      printf("Number of superpixels: %d\n", nSp);
      if (doBR)
        printf("Boundary recall: %f\n", boundaryRecall);
      if (doBP)
        printf("Boundary precision: %f\n", boundaryPrecision);
      if (doBF)
        printf("Boundary FScore: %f\n", boundaryFScore);
      if (isMultiGt && doUnion) {
	if (doBR)
          printf("Union Boundary recall: %f\n", brUnion);
	if (doBP)
          printf("Union Boundary precision: %f\n", bpUnion);
	if (doBF)
          printf("Union Boundary FScore: %f\n", bfUnion);
      }
      if (doASA)
        printf("ASA: %f\n", ASA);
      if (doUE)
        printf("Undersegmentation error: %f\n", underSegmentationError);
      if (doComp) {
        if (!is3D)
          printf("Compactness: %f\n", comp);
	else
          printf("Compactness: not available for 3D images\n");
      }
      if (doTop)
        printf("Topology: %f\n", top);
      if (doCover)
        printf("Covering: %f\n", covering);
      if (doVOI)
        printf("VOI: %f\n", VOI);
      if (doPRI)
        printf("PRI: %f\n", PRI);
      if (doDice) {
        printf("Average dice: %f\n", dice[0]);
        for (int i = 0; i < nObjects; ++i)
          printf("Object %d dice: %f\n", i+1, dice[i+1]);
      }
      break;
    case 2: // CSV stats
      printf("%d", nSp);
      if (doBR)
        printf(",%f", boundaryRecall);
      if (isMultiGt && doUnion && doBR)
        printf(",%f", brUnion);
      if (doBP)
        printf(",%f", boundaryPrecision);
      if (isMultiGt && doUnion && doBP)
        printf(",%f", bpUnion);
      if (doBF)
        printf(",%f", boundaryFScore);
      if (isMultiGt && doUnion && doBF)
        printf(",%f", bfUnion);
      if (doASA)
        printf(",%f", ASA);
      if (doUE)
	printf(",%f", underSegmentationError);
      if (doComp)
        printf(",%f", comp);
      if (doTop)
        printf(",%f", top);
      if (doCover)
        printf(",%f", covering);
      if (doVOI)
        printf(",%f", VOI);
      if (doPRI)
        printf(",%f", PRI);
      if (doDice)
        for (int i = 0; i < nObjects + 1; ++i)
          printf(",%f", dice[i]);
      break;
    default:
      break;
  }

  // Clean up
  iftDestroyImage(&segLabelMap);
  iftDestroyImage(&segBorders);
  if (doDice)
    free(dice);
  iftDestroyDir(&gtDir);

  return 0;
}

iftDict *iftGetArguments(int argc, const char *argv[]) {
  iftCmdLineOpt cmd_line_opts[] = {
    {.short_name = "-i", .long_name = SEG_ARG, .has_arg = true,
      .arg_type = IFT_STR_TYPE, .required = true, 
      .help = "Input segmentation path."},
    {.short_name = "-g", .long_name = GT_ARG, .has_arg = true,
      .arg_type = IFT_STR_TYPE, .required = true, 
      .help = "Input ground-truth path."},
    {.short_name = "-ib", .long_name = ISBORDERSEG_ARG, .has_arg = false,
      .required = false, .help = "Interpret segmentation as border map."},
    {.short_name = "-gb", .long_name = ISBORDERGT_ARG, .has_arg = false,
      .required = false, .help = "Interpret ground-truth as border map."},
    {.short_name = "-p", .long_name = PRINT_ARG, .has_arg = true,
      .arg_type = IFT_LONG_TYPE, .required = false,
      .help = "1 = Human, 2 = CSV, 3,4 = CSV_header (4 to add std field)."},
    {.short_name = "-u", .long_name = UNION_ARG, .has_arg = false,
      .required = false, .help = "Compute border metrics for union of all gt images."},
    {.short_name = "-t", .long_name = TOLERANCE_ARG, .has_arg = true,
      .arg_type = IFT_DBL_TYPE, .required = false,
      .help = "Maximum distance to count as valid border pixel (default = 2.0)."},
    {.short_name = "-br", .long_name = BR_ARG, .has_arg = false,
      .required = false, .help = "Compute boundary recall."},
    {.short_name = "-bp", .long_name = BP_ARG, .has_arg = false,
      .required = false, .help = "Compute boundary precision."},
    {.short_name = "-bf", .long_name = BF_ARG, .has_arg = false,
      .required = false, .help = "Compute boundary f-score."},
    {.short_name = "-ue", .long_name = UE_ARG, .has_arg = false,
      .required = false, .help = "Compute undersegmentation error."},
    {.short_name = "-as", .long_name = ASA_ARG, .has_arg = false,
      .required = false, .help = "Compute achievable segmentation accuracy."},
    {.short_name = "-c", .long_name = COMP_ARG, .has_arg = false,
      .required = false, .help = "Compute compactness."},
    {.short_name = "-sp", .long_name = SPMETRICS_ARG, .has_arg = false,
      .required = false, .help = "Compute standard superpixel metrics (equivalent to -br -ue -as -c)."},
    {.short_name = "-d", .long_name = DICE_ARG, .has_arg = true,
      .arg_type = IFT_DBL_TYPE, .required = false,
      .help = "Compute dice for [param] objects."},
    {.short_name = "-to", .long_name = TOP_ARG, .has_arg = false,
      .required = false, .help = "Compute topology metric."},
    {.short_name = "-cv", .long_name = COVER_ARG, .has_arg = false,
      .required = false, .help = "Compute covering."},
    {.short_name = "-vi", .long_name = VOI_ARG, .has_arg = false,
      .required = false, .help = "Compute variation of information."},
    {.short_name = "-ri", .long_name = PRI_ARG, .has_arg = false,
      .required = false, .help = "Compute probabilistic rand index."},
    {.short_name = "-ls", .long_name = SEGMETRICS_ARG, .has_arg = false,
      .required = false, .help = "Compute standard large segment metrics (equivalent to -cv -vi -ri)."}
  };
  int n_opts = sizeof(cmd_line_opts) / sizeof (iftCmdLineOpt);

  char program_description[IFT_STR_DEFAULT_SIZE] = 
    "This is a general demo to compute image segmentation metrics.";

  // Parser Setup
  iftCmdLineParser *parser =
    iftCreateCmdLineParser(program_description, n_opts, cmd_line_opts);
  iftDict *args = iftParseCmdLine(argc, argv, parser);
  iftDestroyCmdLineParser(&parser);

  return args;
}

