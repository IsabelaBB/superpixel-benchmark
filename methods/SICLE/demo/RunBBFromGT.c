/*****************************************************************************\
* RunBBFromGT.c
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-07-20
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#include "ift.h"
#include "iftArgs.h"
#include "iftMetrics.h"

void usage();
void readImgInputs
(const iftArgs *args, iftImage **gt_img, const char **path);
iftImage *calcMinBoundBoxImage
(const iftImage *gt_img);

int main(int argc, char const *argv[])
{
  //-------------------------------------------------------------------------//
  bool has_req, has_help;
  iftArgs *args;

  args = iftCreateArgs(argc, argv);

  has_req = iftExistArg(args, "gt") && iftExistArg(args, "out");
  has_help = iftExistArg(args, "help");

  if(has_req == false || has_help == true)
  {
    usage(); 
    iftDestroyArgs(&args);
    return EXIT_FAILURE;
  }
  //-------------------------------------------------------------------------//
  const char *OUT_PATH;
  int min_val;
  iftBoundingBox bb;
  iftImage *gt_img, *out_img;

  readImgInputs(args, &gt_img, &OUT_PATH);
  iftDestroyArgs(&args);
  
  min_val = iftMinimumValue(gt_img);
  #if IFT_OMP //-------------------------------------------------------------//
  #pragma omp parallel for
  #endif //------------------------------------------------------------------//
  for(int p = 0; p < gt_img->n; ++p)
    if(gt_img->val[p] == min_val) gt_img->val[p] = 0;

  bb = iftMinBoundingBox(gt_img, NULL);
  out_img = iftCreateImage(gt_img->xsize, gt_img->ysize, gt_img->zsize);
  iftDestroyImage(&gt_img);
  
  #if IFT_OMP //-------------------------------------------------------------//
  #pragma omp parallel for
  #endif //------------------------------------------------------------------//
  for(int z = bb.begin.z; z <= bb.end.z; ++z)
  {
    for(int y = bb.begin.y; y <= bb.end.y; ++y)
    {
      for(int x = bb.begin.x; x <= bb.end.x; ++x)
      {
        int p_idx;
        iftVoxel p_vxl;
        
        p_vxl.x = x; p_vxl.y = y; p_vxl.z = z;
        p_idx = iftGetVoxelIndex(out_img, p_vxl);

        out_img->val[p_idx] = 255;
      }
    }
  }

  iftWriteImageByExt(out_img, OUT_PATH);
  iftDestroyImage(&out_img);

  return EXIT_SUCCESS;
}

void usage()
{
  const int SKIP_IND = 15; // For indentation purposes
  printf("\nThe required parameters are:\n");
  printf("%-*s %s\n", SKIP_IND, "--gt", 
         "Input groundtruth image");
  printf("%-*s %s\n", SKIP_IND, "--out", 
         "Output groundtruth bounding-box image");

  printf("\nThe optional parameters are:\n");
  printf("%-*s %s\n", SKIP_IND, "--help", 
         "Prints this message");

  printf("\n");
}

void readImgInputs
(const iftArgs *args, iftImage **gt_img, const char **path)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(args != NULL);
  assert(gt_img != NULL);
  #endif //------------------------------------------------------------------//
  const char *PATH;

  if(iftHasArgVal(args, "gt") == true) 
  {
    PATH = iftGetArg(args, "gt");

    (*gt_img) = iftReadImageByExt(PATH);
  }
  else iftError("No ground-truth image path was given", "readImgInputs");

  if(iftHasArgVal(args, "out") == true)
  {
    PATH = iftGetArg(args, "out");

    if(iftIsImagePathnameValid(PATH) == true) (*path) = PATH;
    else iftError("Unknown image type", "readImgInputs");
  }
  else iftError("No output image path was given", "readImgInputs");
}

