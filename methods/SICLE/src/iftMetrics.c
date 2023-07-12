/*****************************************************************************\
* iftMetrics.c
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-03-08
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#include "iftMetrics.h"

/*****************************************************************************\
*
*                              PRIVATE FUNCTIONS
*
\*****************************************************************************/
//===========================================================================//
// GENERAL & AUXILIARY
//===========================================================================//
/*
  Counts the number of spels in which each label intersects with an object 
  defined in the ground-truth
*/
int **_iftCalcLabelGTIntersec
(const iftImage *label_img, const iftImage *gt_img)
{
  #ifdef IFT_DEBUG //--------------------------------------------------------//
  assert(label_img != NULL);
  assert(gt_img != NULL);
  iftVerifyImageDomains(label_img, gt_img, "_iftCalcLabelGTIntersec");
  #endif //------------------------------------------------------------------//
  int min_label, max_label, num_labels, min_gt, max_gt, num_gt;
  int **inter;

  iftMinMaxValues(label_img, &min_label, &max_label);
  num_labels = max_label - min_label + 1;
  iftMinMaxValues(gt_img, &min_gt, &max_gt);
  num_gt = max_gt - min_gt + 1;

  inter = calloc(num_labels, sizeof(int*));
  assert(inter != NULL);

  for(int i = 0; i < num_labels; ++i) // For each label
  {
    inter[i] = calloc(num_gt, sizeof(int));
    assert(inter[i] != NULL);
  }

  for(int p = 0; p < label_img->n; ++p) // For each spel
  {
    int label, gt;

    label = label_img->val[p] - min_label;
    gt = gt_img->val[p] - min_gt;

    inter[label][gt]++; // Add one more spel to the intersection
  }

  return inter;
}

/*****************************************************************************\
*
*                               PUBLIC FUNCTIONS
*
\*****************************************************************************/
float iftBoundRecall
(const iftImage *label_img, const iftImage *gt_img)
{
  #ifdef IFT_DEBUG //--------------------------------------------------------//
  assert(label_img != NULL);
  assert(gt_img != NULL);
  iftVerifyImageDomains(label_img, gt_img, "iftBoundRecall");
  #endif //------------------------------------------------------------------//
  int tp, fn, r;
  float bound_rec; 
  iftAdjRel *A;

  r = ceil(0.0025 * iftDiagonalSize(label_img)); // See [1] for more details
  // The adjacency below is equivalent to a square neighborhood of 2r+1 x 2r+1
  A = iftCircular(r * sqrtf(2.0));

  tp = fn = 0;
  #ifdef IFT_OMP //----------------------------------------------------------//
  #pragma omp parallel for reduction(+:tp, fn)
  #endif //------------------------------------------------------------------//
  for(int p = 0; p < label_img->n; ++p) // For each pixel
  {
    bool is_label_border, is_gt_border;
    iftVoxel p_vxl;
    
    p_vxl = iftGetVoxelCoord(label_img, p);

    is_gt_border = is_label_border = false;//Set as a non-border pixel(for now)
    for(int i = 1; i < A->n; ++i) // For each adjacent pixel (i = 0 == p)
    {
      iftVoxel adj_vxl;

      adj_vxl = iftGetAdjacentVoxel(A, p_vxl, i);

      // If it is within the image domain
      if(iftValidVoxel(label_img, adj_vxl) == true) 
      {
        int adj_idx;

        adj_idx = iftGetVoxelIndex(label_img, adj_vxl);

        // If in the GT image, p and adj have different labels
        if(gt_img->val[p] != gt_img->val[adj_idx]) is_gt_border = true;
        // If in the label image, p and adj have different labels
        if(label_img->val[p] != label_img->val[adj_idx]) is_label_border = true;
      }
    }

    if(is_gt_border == true) // If p is a GT border pixel
    {
      // If it is also a border spel in the label image
      if(is_label_border == true) tp++; // It was correctly assigned as border
      else fn++; // It was incorrectly assigned as non-border
    } // Then, don't care for true negatives or false positives 
  }

  bound_rec = tp/(float)(tp + fn); // Compute recall

  iftDestroyAdjRel(&A);

  return bound_rec;
}

float iftUnderSegmError
(const iftImage *label_img, const iftImage *gt_img)
{
  #ifdef IFT_DEBUG //--------------------------------------------------------//
  assert(label_img != NULL);
  assert(gt_img != NULL);
  iftVerifyImageDomains(label_img, gt_img, "iftUnderSegmError");
  #endif //------------------------------------------------------------------//
  int min_label, max_label, num_labels, min_gt, max_gt, num_gt;
  float under_segm;
  int **inter;

  iftMinMaxValues(label_img, &min_label, &max_label);
  num_labels = max_label - min_label + 1;
  iftMinMaxValues(gt_img, &min_gt, &max_gt);
  num_gt = max_gt - min_gt + 1;

  // Calculate the intersection between the superpixels and the GT objects
  inter = _iftCalcLabelGTIntersec(label_img, gt_img);

  under_segm = 0.0;
  #ifdef IFT_OMP //----------------------------------------------------------//
  #pragma omp parallel for reduction(+:under_segm)
  #endif //------------------------------------------------------------------//
  for(int i = 0; i < num_labels; ++i) // For each label
  {
    int error, size;

    size = 0;
    for(int j = 0; j < num_gt; ++j) // Get the label's size (in pixels)
      size += inter[i][j]; // Accumulate its intersection with the GT object

    error = 0;
    for(int j = 0; j < num_gt; ++j) // For each GT object
    {
      if(inter[i][j] > 0)//If exists at least one spel in their intersection
        // Gets the lowest error between inner and outer leakings
        // |Sj - Gi| = |Sj| - |Sj intersec Gi| in this case
        error += iftMin(inter[i][j], size - inter[i][j]);
    }

    under_segm += error; // Accumulate the label's leaking errors
  }
  under_segm /= (float)label_img->n; // Normalize to [0,1]

  for(int i = 0; i < num_labels; ++i)
    free(inter[i]);
  free(inter);

  return under_segm;
}
