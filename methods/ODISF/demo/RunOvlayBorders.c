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

typedef struct pointCut
{
  float val[2];
} pointCut;

void usage();
iftImage *ovlayBorders(const iftImage *orig_img, const iftImage *label_img, const float thick,
                       const iftFColor rgb);
void readImgInputs(const iftArgs *args, iftImage **img, iftImage **labels, const char **path);
void readOptArgs(const iftArgs *args, float *thick, iftFColor *rgb);
void writeOvlayImage(const iftImage *ovlay_img, const char *path);

iftImage *merge(const iftImage *ovlay_img1, const iftImage *ovlay_img2, const float thick, pointCut point);
void readImgInputsMult(const iftArgs *args, iftImage **img, iftImage **labels1, iftImage **labels2, pointCut *point, const char **path);

int main(int argc, char const *argv[])
{
  //-------------------------------------------------------------------------//
  bool has_req, has_help;
  iftArgs *args;

  args = iftCreateArgs(argc, argv);

  has_req = iftExistArg(args, "img") && iftExistArg(args, "labels") &&
            iftExistArg(args, "out");
  has_help = iftExistArg(args, "help");

  if (has_req == false || has_help == true)
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

  if (!iftExistArg(args, "mult"))
  {
    readImgInputs(args, &img, &label_img, &OVLAY_PATH);
    readOptArgs(args, &thick, &rgb);
    iftDestroyArgs(&args);

    ovlay_img = ovlayBorders(img, label_img, thick, rgb);
    iftDestroyImage(&img);
    iftDestroyImage(&label_img);

    iftWriteImageByExt(ovlay_img, OVLAY_PATH);
    iftDestroyImage(&ovlay_img);
  }
  else
  {
    iftImage *label_img2, *ovlay_img2, *result_img;
    pointCut point;

    readImgInputsMult(args, &img, &label_img, &label_img2, &point, &OVLAY_PATH);
    readOptArgs(args, &thick, &rgb);
    iftDestroyArgs(&args);

    ovlay_img = ovlayBorders(img, label_img, thick, rgb);
    ovlay_img2 = ovlayBorders(img, label_img2, thick, rgb);
    iftDestroyImage(&img);
    iftDestroyImage(&label_img);
    iftDestroyImage(&label_img2);
    result_img = merge(ovlay_img, ovlay_img2, thick, point);
    iftWriteImageByExt(result_img, OVLAY_PATH);
    iftDestroyImage(&ovlay_img);
    iftDestroyImage(&ovlay_img2);
    iftDestroyImage(&result_img);
  }

  return EXIT_SUCCESS;
}

void usage()
{
  const int SKIP_IND = 15; // For indentation purposes
  printf("\nThe required parameters are:\n");
  printf("%-*s %s\n", SKIP_IND, "--img",
         "Input image");
  printf("%-*s %s\n", SKIP_IND, "--labels",
         "Input label");
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

iftImage *ovlayBorders(const iftImage *orig_img, const iftImage *label_img, const float thick,
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
  ovlay_img = iftCreateColorImage(orig_img->xsize, orig_img->ysize,
                                  orig_img->zsize, depth);

  RGB.val[0] = rgb.val[0] * norm_val;
  RGB.val[1] = rgb.val[1] * norm_val;
  RGB.val[2] = rgb.val[2] * norm_val;
  YCbCr = iftRGBtoYCbCr(RGB, norm_val);

#if IFT_OMP //-------------------------------------------------------------//
#pragma omp parallel for
#endif //------------------------------------------------------------------//
  for (int p = 0; p < ovlay_img->n; ++p)
  {
    bool is_border;
    int i;
    iftVoxel p_vxl;

    is_border = false;
    p_vxl = iftGetVoxelCoord(ovlay_img, p);

    i = 0;
    while (is_border == false && i < A->n)
    {
      iftVoxel adj_vxl;

      adj_vxl = iftGetAdjacentVoxel(A, p_vxl, i);

      if (iftValidVoxel(ovlay_img, adj_vxl) == true)
      {
        int adj_idx;

        adj_idx = iftGetVoxelIndex(ovlay_img, adj_vxl);

        if (label_img->val[p] != label_img->val[adj_idx])
          is_border = true;
      }

      ++i;
    }

    if (is_border == true)
    {
      ovlay_img->val[p] = YCbCr.val[0];
      ovlay_img->Cb[p] = YCbCr.val[1];
      ovlay_img->Cr[p] = YCbCr.val[2];
    }
    else
    {
      ovlay_img->val[p] = orig_img->val[p];
      if (iftIsColorImage(orig_img) == true)
      {
        ovlay_img->Cb[p] = orig_img->Cb[p];
        ovlay_img->Cr[p] = orig_img->Cr[p];
      }
    }
  }
  iftDestroyAdjRel(&A);

  if (depth != 8)
    iftConvertNewBitDepth(&ovlay_img, 8);

  return ovlay_img;
}

void readImgInputs(const iftArgs *args, iftImage **img, iftImage **labels, const char **path)
{
#if IFT_DEBUG //-----------------------------------------------------------//
  assert(args != NULL);
  assert(img != NULL);
  assert(labels != NULL);
  assert(path != NULL);
#endif //------------------------------------------------------------------//
  const char *PATH;

  if (iftHasArgVal(args, "img") == true)
  {
    PATH = iftGetArg(args, "img");

    (*img) = iftReadImageByExt(PATH);
  }
  else
    iftError("No image path was given", "readImgInputs");

  if (iftHasArgVal(args, "labels") == true)
  {
    PATH = iftGetArg(args, "labels");

    (*labels) = iftReadImageByExt(PATH);
  }
  else
    iftError("No label image path was given", "readImgInputs");

  iftVerifyImageDomains(*img, *labels, "readImgInputs");

  if (iftHasArgVal(args, "out") == true)
    (*path) = iftGetArg(args, "out");
  else
    iftError("No output image path was given", "readImgInputs");
}

void readOptArgs(const iftArgs *args, float *thick, iftFColor *rgb)
{
#if IFT_DEBUG //-----------------------------------------------------------//
  assert(args != NULL);
  assert(thick != NULL);
  assert(rgb != NULL);
#endif //------------------------------------------------------------------//

  if (iftExistArg(args, "thick") == true)
  {
    if (iftHasArgVal(args, "thick") == true)
      (*thick) = atoi(iftGetArg(args, "thick"));
    else
      iftError("No border thickness was given", "readOptArgs");
  }
  else
    (*thick) = 1.0;

  if (iftExistArg(args, "rgb") == true)
  {
    if (iftHasArgVal(args, "rgb") == true)
    {
      const char *VAL;
      int i;
      char *tmp, *tok;

      VAL = iftGetArg(args, "rgb");
      tmp = iftCopyString(VAL);
      tok = strtok(tmp, ",");

      i = 0;
      while (tok != NULL && i < 3)
      {
        float c;

        c = atof(tok);
        if (c >= 0 && c <= 1)
          (*rgb).val[i] = c;
        else
          iftError("The color should be within [0,1]", "readOptArgs");

        tok = strtok(NULL, ",");
        i++;
      }
      if ((tok != NULL && i == 3) || (tok == NULL && i < 2))
        iftError("Three colors are required for the RGB", "readOptArgs");

      free(tmp);
    }
    else
      iftError("No normalized RGB color was given", "readOptArgs");
  }
  else
    (*rgb).val[0] = (*rgb).val[1] = (*rgb).val[2] = 0.0;
}

void readImgInputsMult(const iftArgs *args, iftImage **img, iftImage **labels1, iftImage **labels2, pointCut *point, const char **path)
{
#if IFT_DEBUG //-----------------------------------------------------------//
  assert(args != NULL);
  assert(img != NULL);
  assert(labels != NULL);
  assert(path != NULL);
#endif //------------------------------------------------------------------//
  const char *PATH;

  if (iftHasArgVal(args, "img") == true)
  {
    PATH = iftGetArg(args, "img");

    (*img) = iftReadImageByExt(PATH);
  }
  else
    iftError("No image path was given", "readImgInputs");

  if (iftHasArgVal(args, "labels") == true)
  {
    PATH = iftGetArg(args, "labels");

    (*labels1) = iftReadImageByExt(PATH);
  }
  else
    iftError("No label image path was given", "readImgInputs");

  if (iftHasArgVal(args, "labels2") == true)
  {
    PATH = iftGetArg(args, "labels2");

    (*labels2) = iftReadImageByExt(PATH);
  }
  else
    iftError("No label2 image path was given", "readImgInputs");

  if (iftHasArgVal(args, "point") == true)
  {

    const char *VAL;
    int i;
    char *tmp, *tok;

    VAL = iftGetArg(args, "point");
    tmp = iftCopyString(VAL);
    tok = strtok(tmp, ",");

    i = 0;
    while (tok != NULL && i < 3)
    {
      float c;

      c = atof(tok);
      (*point).val[i] = c;

      tok = strtok(NULL, ",");
      i++;
    }
    if ((tok == NULL && i != 2))
      iftError("Two values are required to point cut", "readOptArgs");

    free(tmp);
  }
  else
    (*point).val[0] = (*point).val[1] = 50.0;

  iftVerifyImageDomains(*img, *labels1, "readImgInputs");
  iftVerifyImageDomains(*img, *labels2, "readImgInputs");

  if (iftHasArgVal(args, "out") == true)
    (*path) = iftGetArg(args, "out");
  else
    iftError("No output image path was given", "readImgInputs");
}



int check (int *v, int n){
int i;
for (i=0; i<n; i++){
    if (v[i] != -1){
        return -1;
        break;
    }
    else return 1;
}
}


void sortVector (int *v, int n)
{
    int i, k, j=0, vp[n];
    int maxIndex=0;//you need to add this variable in order to keep track of the maximum value in each iteration

    while (check(v,n) == -1)
    {
        for (i=0; i<n; i++)
        {
            maxIndex=i; //you suppose that the maximum is the first element in each loop 
            for (k=i+1; k<n; k++)
            {
                if (v[k] < v[maxIndex])
                  maxIndex=k; // if there is another element greater you preserve its index in the variable 
            }
           //after finishing the loop above you have the greatest variable in the array  which has the index stored in maxIndex
                    vp[i] = v[maxIndex]; // put it in vp array
                    v[maxIndex]=v[i];//put it in treated elements zone
                    v[i]=-1;// make it -1
                    j++;
        }
    }

    for (i=0; i<n; i++)
        v[i] = vp[i];
}


iftImage *merge(const iftImage *ovlay_img1, const iftImage *ovlay_img2, const float thick, pointCut point)
{
#if IFT_DEBUG //-----------------------------------------------------------//
  assert(orig_img != NULL);
  assert(label_img != NULL);
  iftVerifyImageDomains(orig_img, label_img, "ovlayBorders");
  assert(thick > 0);
#endif //------------------------------------------------------------------//
  int depth, norm_val;
  iftImage *result_img;
  iftAdjRel *A;
  iftColor RGB, YCbCr;

  A = iftCircular(4.0);

  depth = iftImageDepth(ovlay_img1);
  norm_val = iftMaxImageRange(depth);
  result_img = iftCreateColorImage(ovlay_img1->xsize, ovlay_img1->ysize,
                                   ovlay_img1->zsize, depth);

  RGB.val[0] = RGB.val[1] = RGB.val[2] = 255;
  YCbCr = iftRGBtoYCbCr(RGB, norm_val);

  // draw line: https://stackoverflow.com/questions/64561090/algorithm-for-drawing-diagonals-in-a-picture
  /*
    aspect = (float)c/r
    for (int i = 0; i < r; i++)
      matrix[i][(int)(aspect*i)] = color_of_diagonal
  */

  bool usingImg1 = true;

  int x0 = point.val[0];
  int y0 = 0;
  int x1 = 0;
  int y1 = point.val[1];

  int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
  int dy = abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
  int err = (dx > dy ? dx : -dy) / 2, e2;

  int pointsLine[ovlay_img1->xsize*ovlay_img1->xsize];

  int numPoints = 0;
  for (;;)
  {
    int p = x0 * ovlay_img1->xsize + y0;

    pointsLine[numPoints] = p;
    numPoints++;

    if (x0 == x1 && y0 == y1)
      break;

    e2 = err;
    if (e2 > -dx)
    {
      err -= dy;
      x0 += sx;
    }
    if (e2 < dy)
    {
      err += dx;
      y0 += sy;
    }

    
  }

  sortVector(pointsLine, numPoints);

  int pointId = 0;
//#if IFT_OMP //-------------------------------------------------------------//
//#pragma omp parallel for
//#endif //------------------------------------------------------------------//
  for (int p = 0; p < result_img->n; ++p)
  {
    int x = p % ovlay_img1->xsize;
    //int y = p / ovlay_img1->xsize;

    if (x == 0)
      usingImg1 = true;

    if(pointId == numPoints - 1)
      usingImg1 = false;
    
    if (pointsLine[pointId] == p && (pointId < numPoints - 1))
    {
      //result_img->val[p] = YCbCr.val[0];
      //result_img->Cb[p] = YCbCr.val[1];
      //result_img->Cr[p] = YCbCr.val[2];
      usingImg1 = false;
      pointId++;
    }
    else
    {

      if (usingImg1)
      {
        result_img->val[p] = ovlay_img1->val[p];
        result_img->Cb[p] = ovlay_img1->Cb[p];
        result_img->Cr[p] = ovlay_img1->Cr[p];
      }
      else
      {
        result_img->val[p] = ovlay_img2->val[p];
        result_img->Cb[p] = ovlay_img2->Cb[p];
        result_img->Cr[p] = ovlay_img2->Cr[p];
      }
    }
  }

  
  for (pointId=0; pointId < numPoints; pointId++)
  {

    int p = pointsLine[pointId];
    result_img->val[p] = YCbCr.val[0];
    result_img->Cb[p] = YCbCr.val[1];
    result_img->Cr[p] = YCbCr.val[2];

    int i;
    iftVoxel p_vxl;
    p_vxl = iftGetVoxelCoord(ovlay_img1, p);

    i = 0;
    while (i < A->n)
    {
      iftVoxel adj_vxl;

      adj_vxl = iftGetAdjacentVoxel(A, p_vxl, i);

      if (iftValidVoxel(ovlay_img1, adj_vxl) == true)
      {
        int adj_idx;

        adj_idx = iftGetVoxelIndex(ovlay_img1, adj_vxl);

        result_img->val[adj_idx] = YCbCr.val[0];
        result_img->Cb[adj_idx] = YCbCr.val[1];
        result_img->Cr[adj_idx] = YCbCr.val[2];
      }

      ++i;
    }
  }

  iftDestroyAdjRel(&A);

  if (depth != 8)
    iftConvertNewBitDepth(&result_img, 8);

  return result_img;
}