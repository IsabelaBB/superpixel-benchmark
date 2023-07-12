/**
 * @file
 * @brief General RISF program.
 * @author Felipe Lemes Galvao
 *
 * New version of iftRISF_segmentation with better argument
 *   handling and correct name (RISF).
 */

#include <ift.h>
#include "iftSuperpixelFeatures.h"
#include "iftExperimentUtility.h"
#include "iftRISF.h"

/* TODO fix memory leaks when reading strings from dictionary. */

#define SAMPLING_ARG "--sampling"
#define SEED_RECOMP_ARG "--seed-recomp"
#define PATH_COST_ARG "--path-cost"
#define INPUT_ARG "--input-image"
#define OUTPUT_ARG "--output-image"
#define NSP_ARG "--superpixel-num"
#define ALPHA_ARG "--alpha"
#define BETA_ARG "--beta"
#define ITER_ARG "--iters-num"
#define PRINT_ARG "--print-opt"
#define PREVSEG_ARG "--prev-seg"
#define FORMAT_ARG "--output-format"

iftDict *iftGetArguments(int argc, const char *argv[]);

int main(int argc, char* argv[])
{
  iftDict* args = iftGetArguments(argc, (const char **) argv);
  char *tmpPath;

  int printOpt = (iftDictContainKey(PRINT_ARG, args, NULL)) ?
    iftGetLongValFromDict(PRINT_ARG, args) : 1;
  bool isHierarchic = iftDictContainKey(PREVSEG_ARG, args, NULL);

  // Load Image
  tmpPath = iftGetStrValFromDict(INPUT_ARG, args);
  iftImage *img = iftReadImageByExt(tmpPath);
  free(tmpPath);
  bool is3D = iftIs3DImage(img);

  // Prepare igraph reference image TODO Move to a separate function
  iftImage *refImg = NULL;
  if (isHierarchic) {
    // Load pre-computed superpixel segmentation (previous RISF level)
    tmpPath = iftGetStrValFromDict(PREVSEG_ARG, args);
    refImg = iftReadImageByExt(tmpPath);
    free(tmpPath);
    int minLabel = iftMinimumValue(refImg);
    if (minLabel <= 0)
      for (int p = 0; p < refImg->n; ++p)
        refImg->val[p] += (1 - minLabel);
  } else {
    // Make superpixel graph act as a pixel-based one
    refImg = iftSelectImageDomain(img->xsize, img->ysize, img->zsize);
    for (int i = 0; i < refImg->n; ++i)
      refImg->val[i] = i + 1;
  }

  // Build igraph
  iftSuperpixelIGraph *igraph = iftInitSuperpixelIGraph(refImg);
  iftAdjRel *A = is3D ? iftSpheric(1.0) : iftCircular(1.0);
  if (isHierarchic)
    iftSetSuperpixelIGraphExplicitRAG(igraph, A);
  else
    iftSetSuperpixelIGraphImplicitAdjacency(igraph, A);

  // Set superpixel parametric features
  float alpha = (iftDictContainKey(ALPHA_ARG, args, NULL)) ?
    iftGetDblValFromDict(ALPHA_ARG, args) : 
    iftIsColorImage(img) ? 0.2 : 40.96;
  float beta = (iftDictContainKey(BETA_ARG, args, NULL)) ?
    iftGetDblValFromDict(BETA_ARG, args) : 12.0;
  int featOpt = 1;
  iftMatrix *pFeats = NULL;
  switch (featOpt) {
    case 1: // Color mean
    default:
      {
        if (iftIsColorImage(img))
          pFeats = iftComputeSuperpixelFeaturesByColorSpaceMean(igraph->refImg, img, LABNorm_CSPACE);  
        else {
          iftSetCbCr(img, (iftMaxImageRange(iftImageDepth(img))+1)/2);
          pFeats = iftComputeSuperpixelFeaturesByColorSpaceMean(igraph->refImg, img, LABNorm_CSPACE);  
        }
        iftSetSuperpixelIGraphFeatures(igraph, pFeats, iftDistance1, alpha, beta);
      }
      break;
    case 2: // 20^3-bins LAB histogram XiSquare dist
      pFeats = iftComputeSuperpixelFeaturesByColorHistogram(img, igraph->refImg,  20);
      for (int node = 0; node < pFeats->nrows; ++node)
        iftNormalizeFeatures(iftMatrixRowPointer(pFeats, node), pFeats->ncols);
      iftSetSuperpixelIGraphFeatures(igraph, pFeats, iftDistance10, alpha, beta);
      break;
  }

  // Set RISF parameters
  iftISFSampling sampling = (iftDictContainKey(SAMPLING_ARG, args, NULL)) ?
    iftGetLongValFromDict(SAMPLING_ARG, args) : 1;
  iftISFSeedRecomp seedRecomp = (iftDictContainKey(SEED_RECOMP_ARG, args, NULL)) ?
    iftGetLongValFromDict(SEED_RECOMP_ARG, args) : 2;
  iftISFPathCost pathCost = (iftDictContainKey(PATH_COST_ARG, args, NULL)) ?
    iftGetLongValFromDict(PATH_COST_ARG, args) : 1;
  int nSuperpixels = iftGetLongValFromDict(NSP_ARG, args);
  int nIters = (iftDictContainKey(ITER_ARG, args, NULL)) ?
    iftGetLongValFromDict(ITER_ARG, args) : 10;
  int nSmoothIters = 0; // For now, no smoothing over superpixels
  float smallSpThreshold = 0.00;
  if (sampling == IFT_ISF_GEODESIC_SAMPLING) {
    //nSuperpixels *= 1.2;
    //smallSpThreshold = 0.05; // Default = 0.05
  }
  iftISFOptions *optISF = iftInitISFOptions(nSuperpixels, nIters,
      nSmoothIters, smallSpThreshold, sampling, pathCost, seedRecomp);

  // Actually compute RISF
  timer *t1 = iftTic();
  iftImage *superpixelLabelMap = iftSuperpixelSegmentationByRISF(igraph, optISF, img);
  timer *t2 = iftToc();
  if (printOpt == 1)
    printf("Processing time: %fms\n", iftCompTime(t1, t2));
  else if (printOpt == 2)
    printf("%f", iftCompTime(t1, t2));

  // Write segmentation result to file
  tmpPath = iftGetStrValFromDict(OUTPUT_ARG, args);
  int outputFormat = (iftDictContainKey(FORMAT_ARG, args, NULL)) ?
    iftGetLongValFromDict(FORMAT_ARG, args) : 0;
  switch (outputFormat) {
    case 1:
      iftWriteImageP2(superpixelLabelMap, tmpPath);
      break;
    case 2:
      iftWriteImageP5(superpixelLabelMap, tmpPath);
      break;
    default:
    case 0:
      iftWriteImageByExt(superpixelLabelMap, tmpPath);
      break;
  }
  free(tmpPath);

  iftDestroyImage(&img);
  iftDestroyImage(&refImg);
  iftDestroySuperpixelIGraph(&igraph);
  iftDestroyAdjRel(&A);
  iftDestroyMatrix(&pFeats);
  iftDestroyISFOptions(&optISF);
  iftDestroyImage(&superpixelLabelMap);

  return 0;
}

#define SAMPLING_HELP "Initial sampling method:\n\
    1 = Grid (DEFAULT)\n\
    2 = Mixed\n\
    3 = Random\n\
    4 = Geodesic\n\
    5 = Grid on Mask"

#define SEED_RECOMP_HELP "Seed recomputation method:\n\
    1 = Medoid - node with feature vector closest to mean\n\
    2 = Centroid\n\
    DEFAULT based on path cost, medoid for root and centroid for mean."

#define PATH_COST_HELP "Non-trivial path cost function:\n\
    1 = Root (DEFAULT)\n\
    2 = Mean"

#define ALPHA_HELP "Linear weight for parametric distance.\n\
    Smaller is more compact. Recommended range [0.12,0.5] for natural\n\
    images. Default is 0.2."

#define BETA_HELP "Exponential weight for parametric distance\n\
    Default is 12."

#define PRINT_HELP "Printing options:\n\
    0 = Nothing is printed\n\
    1 = Human-readable basic stats (DEFAULT)\n\
    2 = CSV formatted stats"

#define OUTPUT_FORMAT_HELP "Label map output format:\n\
    0 = Deduce by extension (DEFAULT)\n\
    1 = PGM P2 - ASCII\n\
    2 = PGM P5 - binary"

iftDict *iftGetArguments(int argc, const char *argv[]) {
  iftCmdLineOpt cmd_line_opts[] = {
    {.short_name = "-s", .long_name = SAMPLING_ARG, .has_arg = true,
      .arg_type = IFT_LONG_TYPE, .required = false,
      .help = SAMPLING_HELP},
    {.short_name = "-r", .long_name = SEED_RECOMP_ARG, .has_arg = true,
      .arg_type = IFT_LONG_TYPE, .required = false,
      .help = SEED_RECOMP_HELP},
    {.short_name = "-c", .long_name = PATH_COST_ARG, .has_arg = true,
      .arg_type = IFT_LONG_TYPE, .required = false,
      .help = PATH_COST_HELP},
    {.short_name = "-i", .long_name = INPUT_ARG, .has_arg = true,
      .arg_type = IFT_STR_TYPE, .required = true, 
      .help = "Input image path."},
    {.short_name = "-o", .long_name = OUTPUT_ARG, .has_arg = true,
      .arg_type = IFT_STR_TYPE, .required = true,
      .help = "Output label map image."},
    {.short_name = "-n", .long_name = NSP_ARG, .has_arg = true,
      .arg_type = IFT_LONG_TYPE, .required = true,
      .help = "Target number of superpixels."},
    {.short_name = "-a", .long_name = ALPHA_ARG, .has_arg = true,
      .arg_type = IFT_DBL_TYPE, .required = false,
      .help = ALPHA_HELP},
    {.short_name = "-b", .long_name = BETA_ARG, .has_arg = true,
      .arg_type = IFT_DBL_TYPE, .required = false,
      .help = BETA_HELP},
    {.short_name = "-t", .long_name = ITER_ARG, .has_arg = true,
      .arg_type = IFT_LONG_TYPE, .required = false,
      .help = "Number of IFT iterations."},
    {.short_name = "-p", .long_name = PRINT_ARG, .has_arg = true,
      .arg_type = IFT_LONG_TYPE, .required = false,
      .help = PRINT_HELP},
    {.short_name = "-l", .long_name = PREVSEG_ARG, .has_arg = true,
      .arg_type = IFT_STR_TYPE, .required = true,
      .help = "Previous segmentation label map.\n    See README for pixel level operation.\n"},
    {.short_name = "-x", .long_name = FORMAT_ARG, .has_arg = true,
      .arg_type = IFT_LONG_TYPE, .required = false,
      .help = OUTPUT_FORMAT_HELP}
  };
  int n_opts = sizeof(cmd_line_opts) / sizeof (iftCmdLineOpt);

  char program_description[IFT_STR_DEFAULT_SIZE] = 
    "This is a general demo to run RISF methods.";

  // Parser Setup
  iftCmdLineParser *parser =
    iftCreateCmdLineParser(program_description, n_opts, cmd_line_opts);
  iftDict *args = iftParseCmdLine(argc, argv, parser);
  iftDestroyCmdLineParser(&parser);

  return args;
}

