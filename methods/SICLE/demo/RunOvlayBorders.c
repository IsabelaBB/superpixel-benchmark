/*****************************************************************************\
* RunOvlayBorders.c
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-02-28
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#include "ift.h"
#include "iftArgs.h"

void usage();
iftImage *ovlayBorders
(const iftImage *orig_img, const iftImage *label_img, const float thick, 
 const iftFColor rgb);
void readImgInputs
(const iftArgs *args, iftImage **img, iftImage **labels, const char **path);
void readOptArgs
(const iftArgs *args, float *thick, iftFColor *rgb);

int main(int argc, char const *argv[])
{
  //-------------------------------------------------------------------------//
  bool has_req, has_help;
  iftArgs *args;

  args = iftCreateArgs(argc, argv);

  has_req = iftExistArg(args, "img")  && iftExistArg(args, "labels") &&
            iftExistArg(args, "out");
  has_help = iftExistArg(args, "help");

  if(has_req == false || has_help == true)
  {
    usage(); 
    iftDestroyArgs(&args);
    return EXIT_FAILURE;
  }
  //-------------------------------------------------------------------------//
  const char *OVLAY_PATH;
  float thick;
  iftImage *img, *label_img, *ovlay_img;
  iftFColor rgb;

  readImgInputs(args, &img, &label_img, &OVLAY_PATH);
  readOptArgs(args, &thick, &rgb);
  iftDestroyArgs(&args);

  ovlay_img = ovlayBorders(img, label_img, thick, rgb);
  iftDestroyImage(&img);
  iftDestroyImage(&label_img);

  iftWriteImageByExt(ovlay_img, OVLAY_PATH);
  iftDestroyImage(&ovlay_img);

  return EXIT_SUCCESS;
}

void usage()
{
  const int SKIP_IND = 15; // For indentation purposes
  printf("\nThe required parameters are:\n");
  printf("%-*s %s\n", SKIP_IND, "--img", 
         "Input image");
  printf("%-*s %s\n", SKIP_IND, "--labels", 
         "Input label image");
  printf("%-*s %s\n", SKIP_IND, "--out", 
         "Output border overlayed image");

  printf("\nThe optional parameters are:\n");
  printf("%-*s %s\n", SKIP_IND, "--rgb", 
         "Comma-separated normalized RGB border color. Default: 0,0,0");
  printf("%-*s %s\n", SKIP_IND, "--thick", 
         "Border thickness. Default: 1.0");
  printf("%-*s %s\n", SKIP_IND, "--help", 
         "Prints this message");

  printf("\n");
}

iftImage *ovlayBorders
(const iftImage *orig_img, const iftImage *label_img, const float thick, 
 const iftFColor rgb)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(orig_img != NULL);
  assert(label_img != NULL);
  iftVerifyImageDomains(orig_img, label_img, "ovlayBorders");
  assert(thick > 0);
  #endif //------------------------------------------------------------------//
  int depth, norm_val;
  iftImage *ovlay_img;
  iftAdjRel *A;
  iftColor RGB, YCbCr;

  A = iftCircular(thick);

  depth = iftImageDepth(orig_img);
  norm_val = iftMaxImageRange(depth);
  ovlay_img = iftCreateColorImage(orig_img->xsize, orig_img->ysize, 1, depth);

  RGB.val[0] = rgb.val[0] * norm_val;
  RGB.val[1] = rgb.val[1] * norm_val;
  RGB.val[2] = rgb.val[2] * norm_val;
  YCbCr = iftRGBtoYCbCr(RGB, norm_val);

  #if IFT_OMP //-------------------------------------------------------------//
  #pragma omp parallel for
  #endif //------------------------------------------------------------------//
  for(int p = 0; p < ovlay_img->n; ++p)
  {
    bool is_border;
    int i;
    iftVoxel p_vxl;

    is_border = false;
    p_vxl = iftGetVoxelCoord(ovlay_img, p);

    i = 0;
    while(is_border == false && i < A->n)
    {
      iftVoxel adj_vxl;

      adj_vxl = iftGetAdjacentVoxel(A, p_vxl, i);

      if(iftValidVoxel(ovlay_img, adj_vxl) == true)
      {
        int adj_idx;

        adj_idx = iftGetVoxelIndex(ovlay_img, adj_vxl);

        if(label_img->val[p] != label_img->val[adj_idx]) 
          is_border = true;
      }

      ++i;
    }

    if(is_border == true)
    {
      ovlay_img->val[p] = YCbCr.val[0];
      ovlay_img->Cb[p] = YCbCr.val[1];
      ovlay_img->Cr[p] = YCbCr.val[2];
    }
    else
    {
      ovlay_img->val[p] = orig_img->val[p];
      if(iftIsColorImage(orig_img) == true)
      {
        ovlay_img->Cb[p] = orig_img->Cb[p];
        ovlay_img->Cr[p] = orig_img->Cr[p];
      }

    }
  }
  iftDestroyAdjRel(&A);

  if(depth != 8) iftConvertNewBitDepth(&ovlay_img, 8);

  return ovlay_img;
}

void readImgInputs
(const iftArgs *args, iftImage **img, iftImage **labels, const char **path)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(args != NULL);
  assert(img != NULL);
  assert(labels != NULL);
  assert(path != NULL);
  #endif //------------------------------------------------------------------//
  const char *PATH;

  if(iftHasArgVal(args, "img") == true) 
  {
    PATH = iftGetArg(args, "img");

    (*img) = iftReadImageByExt(PATH);
  }
  else iftError("No image path was given", "readImgInputs");

  if(iftHasArgVal(args, "labels") == true) 
  {
    PATH = iftGetArg(args, "labels");

    (*labels) = iftReadImageByExt(PATH);
  }
  else iftError("No label image path was given", "readImgInputs");

  iftVerifyImageDomains(*img, *labels, "readImgInputs");

  if(iftHasArgVal(args, "out") == true)
  {
    PATH = iftGetArg(args, "out");

    if(iftIsImagePathnameValid(PATH) == true) (*path) = PATH;
    else iftError("Unknown image type", "readImgInputs");
  }
  else iftError("No output image path was given", "readImgInputs");
}

void readOptArgs
(const iftArgs *args, float *thick, iftFColor *rgb)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(args != NULL);
  assert(thick != NULL);
  assert(rgb != NULL);
  #endif //------------------------------------------------------------------//

  if(iftExistArg(args, "thick") == true)
  {
    if(iftHasArgVal(args, "thick") == true) 
      (*thick) = atoi(iftGetArg(args, "thick"));
    else iftError("No border thickness was given", "readOptArgs");
  }
  else (*thick) = 1.0;

  if(iftExistArg(args, "rgb") == true)
  {
    if(iftHasArgVal(args, "rgb") == true)
    {
      const char *VAL;
      int i;
      char *tmp, *tok;


      VAL = iftGetArg(args, "rgb");
      tmp = iftCopyString(VAL);
      tok = strtok(tmp, ",");

      i = 0;
      while(tok != NULL && i < 3)
      {
        float c;

        c = atof(tok);
        if(c >= 0 && c <= 1) (*rgb).val[i] = c;
        else iftError("The color should be within [0,1]", "readOptArgs");

        tok = strtok(NULL, ",");
        i++;
      }
      if((tok != NULL && i == 3) || (tok == NULL && i < 2))
        iftError("Three colors are required for the RGB", "readOptArgs");

      free(tmp);
    }
    else iftError("No normalized RGB color was given", "readOptArgs");
  }
  else (*rgb).val[0] = (*rgb).val[1] = (*rgb).val[2] = 0.0;
}