/*****************************************************************************\
* iftSICLE.c
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-07-08
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#include "iftSICLE.h"

/*****************************************************************************\
*
*                               PRIVATE STRUCTS
*
\*****************************************************************************/
/*
  Contains all IFT's auxiliary data
*/
typedef struct _ift_iftdata
{
  int num_pixels; // Number of pixels within the image
  int *id_map; // Seed id map (a.k.a. root map)
  int *pred_map; // Predecessor map
  double *cost_map; // Accumulated pathcost map
} _iftIFTData;

/*
  Contains the forest's data through its trees (and their relation)
*/
typedef struct _ift_foreststats
{
  int num_trees, num_feats; // Number of trees and the number of its features
  int *tree_size; // Size of each tree (in pixels)
  float *tree_sal; // Mean saliency value of each tree
  float **tree_feats; // Mean feature vector of each tree
  iftBMap **tree_adj; // Tree adjacency relation
} _iftForestStats;

/*****************************************************************************\
*
*                               PUBLIC STRUCTS
*
\*****************************************************************************/
struct ift_sicle_alg
{
  bool use_diag_adj; // Use diagonal adjacents (i.e., 8-adjacency)?
  bool enable_boost; // Boost delineation?
  int n0, nf; // Initial number of seeds and final quantity of superpixels
  int max_iters; // Maximum number of iterations for segmentation
  int num_scales; // Number of superpixel segmentation scales
  int *scales; // Statically decreasing order of superpixel segmentation scales
  float *saliency; // Normalized object saliency values (if they are provided)
  iftMImage *mimg; // Multiband image
  iftBMap *mask; // Mask indicating the ROI (if it exists)
  iftSICLESampl sampl_opt; // Seed sampling option
  iftSICLEArcCost arc_opt; // Arc-cost estimation option
  iftSICLERem rem_opt; // Seed removal criterion option
};

/*****************************************************************************\
*
*                               PRIVATE FUNCTIONS
*
\*****************************************************************************/
//===========================================================================//
// GENERAL
//===========================================================================//
/*
  Computes the L2-norm (a.k.a. Euclidean distance) between two feature vectors
  in the same space (i.e., have the same number of features)
*/
inline double _iftEuclDist
(const float *feat_a, const float *feat_b, const int num_feats)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(feat_a != NULL);
  assert(feat_b != NULL);
  assert(num_feats > 0);
  #endif //------------------------------------------------------------------//
  double dist;

  dist = 0;
  for(int j = 0; j < num_feats; ++j) // L2-norm
    dist += (feat_a[j] - feat_b[j]) * (feat_a[j] - feat_b[j]);
  dist = sqrtf(dist);

  return dist;
}

/*
  Creates a label image given an IFT id map generated during
  segmentation.
*/
iftImage *_iftCreateLabelImageFromIdMap
(const iftSICLE *sicle, const iftIntArray *seeds, const _iftIFTData *iftdata)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(seeds != NULL);
  assert(iftdata != NULL);
  #endif //------------------------------------------------------------------//
  int label;
  iftImage *label_img;

  label_img = iftCreateImage(sicle->mimg->xsize, sicle->mimg->ysize,
                             sicle->mimg->zsize);

  label = 1; // 0 is for the background
  for(long i = 0; i < seeds->n; ++i)
    label_img->val[seeds->val[i]] = label++; // Assign unique label

  #if IFT_OMP //-------------------------------------------------------------//
  #pragma omp parallel for
  #endif //------------------------------------------------------------------//
  for(int p = 0; p < label_img->n; ++p) // For each pixel
  {
    int seed_id;

    seed_id = iftdata->id_map[p];
    if(seed_id != IFT_NIL) // If it is not in the background
      label_img->val[p] = label_img->val[seeds->val[seed_id]];// Copy root label
  }

  return label_img;
}

//===========================================================================//
// IFT DATA
//===========================================================================//
/*
  Creates a new empty IFT data container whose mappings' sizes are defined
  by the number of pixels given
*/
_iftIFTData *_iftCreateIFTData
(const int num_pixels)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(num_pixels > 0);
  #endif //------------------------------------------------------------------//
  _iftIFTData *iftdata;

  iftdata = malloc(sizeof(_iftIFTData));
  assert(iftdata != NULL);

  iftdata->num_pixels = num_pixels;
  iftdata->id_map = calloc(num_pixels, sizeof(int));
  assert(iftdata->id_map != NULL);
  iftdata->pred_map = calloc(num_pixels, sizeof(int));
  assert(iftdata->pred_map != NULL);
  iftdata->cost_map = calloc(num_pixels, sizeof(double));
  assert(iftdata->cost_map != NULL);

  return iftdata;
}

/*
  Deallocates the memory of the given object and sets it to NULL
*/
void _iftDestroyIFTData
(_iftIFTData **iftdata)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(iftdata != NULL);
  assert(*iftdata != NULL);
  #endif //------------------------------------------------------------------//
  free((*iftdata)->id_map);
  free((*iftdata)->pred_map);
  free((*iftdata)->cost_map);

  free(*iftdata);
  (*iftdata) = NULL;
}

/*
  Resets the IFT mappings to their initial (pre-IFT) values
*/
void _iftResetIFTData
(const iftSICLE *sicle, _iftIFTData **iftdata)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(iftdata != NULL);
  assert(*iftdata != NULL);
  #endif //------------------------------------------------------------------//
  #if IFT_OMP //-------------------------------------------------------------//
  #pragma omp parallel for
  #endif //------------------------------------------------------------------//
  for(int p = 0; p < sicle->mimg->n; ++p) // For each pixel
  {
    (*iftdata)->id_map[p] = IFT_NIL;   // If p is not in the background, such
    (*iftdata)->pred_map[p] = IFT_NIL; // values will surely be changed

    // If a mask was not provided OR it is not a forbidden pixel
    if(sicle->mask == NULL || iftBMapValue(sicle->mask, p) == true)
      (*iftdata)->cost_map[p] = IFT_INFINITY_DBL; // Exposed to conquering
    else // Then, it is forbidden
      (*iftdata)->cost_map[p] = IFT_INFINITY_DBL_NEG; // Safe from conquering
  }
}

//===========================================================================//
// FOREST STATISTICS
//===========================================================================//
/*
  Creates an empty IFT forest statistics object based on the IFT data
  generated in the segmentation and its inputs.
*/
_iftForestStats *_iftCreateForestStats
(const iftSICLE *sicle, const _iftIFTData *iftdata, const iftIntArray *seeds)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(iftdata != NULL);
  assert(seeds != NULL);
  #endif //------------------------------------------------------------------//
  _iftForestStats *forstats;

  forstats = malloc(sizeof(_iftForestStats));
  assert(forstats != NULL);

  forstats->num_trees = (int)seeds->n;
  forstats->num_feats = sicle->mimg->m;

  forstats->tree_size = calloc(forstats->num_trees, sizeof(int));
  assert(forstats->tree_size != NULL);

  forstats->tree_sal = calloc(forstats->num_trees, sizeof(float));
  assert(forstats->tree_sal != NULL);

  forstats->tree_feats = calloc(forstats->num_trees, sizeof(float*));
  assert(forstats->tree_feats != NULL);

  forstats->tree_adj = calloc(forstats->num_trees, sizeof(iftBMap*));
  assert(forstats->tree_adj != NULL);

  #if IFT_OMP //-------------------------------------------------------------//
  #pragma omp parallel for
  #endif //------------------------------------------------------------------//
  for(long i = 0; i < seeds->n; ++i) // For each seed/tree
  {
    forstats->tree_feats[i] = calloc(forstats->num_feats, sizeof(float));
    assert(forstats->tree_feats[i] != NULL);

    forstats->tree_adj[i] = iftCreateBMap(forstats->num_trees);
    assert(forstats->tree_adj[i] != NULL);
  }

  return forstats;
}

/*
  Deallocates the memory of the given object and sets it to NULL
*/
void _iftDestroyForestStats
(_iftForestStats **forstats)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(forstats != NULL);
  assert(*forstats != NULL);
  #endif //------------------------------------------------------------------//
  free((*forstats)->tree_size);
  free((*forstats)->tree_sal);
  for(int i = 0; i < (*forstats)->num_trees; ++i) // For each tree
  {
    if((*forstats)->tree_feats[i] != NULL)
      free((*forstats)->tree_feats[i]);
    if((*forstats)->tree_adj[i] != NULL)
      iftDestroyBMap(&((*forstats)->tree_adj[i]));
  }
  free((*forstats)->tree_feats);
  free((*forstats)->tree_adj);

  free(*forstats);
  (*forstats) = NULL;
}

/*
  Calculates the forest statistics object based on the IFT data generated in
  the segmentation and its inputs.
*/
_iftForestStats *_iftCalcForestStats
(const iftSICLE *sicle, const _iftIFTData *iftdata, const iftIntArray *seeds)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(iftdata != NULL);
  assert(seeds != NULL);
  #endif //------------------------------------------------------------------//
  _iftForestStats *forstats;
  iftAdjRel *A;

  forstats = _iftCreateForestStats(sicle, iftdata, seeds);

  // Use the same adjacency as the one used in segmentation
  if(sicle->use_diag_adj == true) A = iftCircular(sqrtf(2.0));
  else A = iftCircular(1.0);

  for(int p = 0; p < sicle->mimg->n; ++p) // For each pixel
  {
    int seed_id;

    seed_id = iftdata->id_map[p];
    if(seed_id != IFT_NIL) // If it is not in the background
    {
      iftVoxel p_vxl;

      p_vxl = iftMGetVoxelCoord(sicle->mimg, p);

      forstats->tree_size[seed_id]++;

      if(sicle->saliency != NULL) // If a saliency map was provided
        forstats->tree_sal[seed_id] += sicle->saliency[p];

      for(int j = 0; j < sicle->mimg->m; ++j) // For each feature
        forstats->tree_feats[seed_id][j] += sicle->mimg->val[p][j];

      for(int i = 1; i < A->n; ++i) // For each adjacent
      {
        iftVoxel adj_vxl;

        adj_vxl = iftGetAdjacentVoxel(A, p_vxl, i);

        // If it is within the image domain
        if(iftMValidVoxel(sicle->mimg, adj_vxl) == true)
        {
          int adj_idx, adj_seed_id;

          adj_idx = iftMGetVoxelIndex(sicle->mimg, adj_vxl);
          adj_seed_id = iftdata->id_map[adj_idx];

          // If the adjacent is not in the background AND their labels differ
          // AND their trees were not already marked as adjacents
          if(adj_seed_id != IFT_NIL && seed_id != adj_seed_id &&
             iftBMapValue(forstats->tree_adj[seed_id], adj_seed_id) == false)
            // They are adjacents
            iftBMapSet1(forstats->tree_adj[seed_id], adj_seed_id);
        }
      }
    }
  }

  for(long i = 0; i < seeds->n; ++i) // For each seed/tree
  {
    if(sicle->saliency != NULL) // If a saliency map was provided
      // Compute mean saliency
      forstats->tree_sal[i] /= (float)forstats->tree_size[i];
    else forstats->tree_sal[i] = 1.0; // Then, all trees are relevant

    for(int j = 0; j < sicle->mimg->m; ++j) // For each feature
      // Compute mean feature
      forstats->tree_feats[i][j] /= (float)forstats->tree_size[i];
  }
  iftDestroyAdjRel(&A);

  return forstats;
}

//===========================================================================//
// SEED SAMPLING
//===========================================================================//
/*
  Calculates the stride within each axis in order to sample the respective
  proportional (to its length) quantity of seeds.
*/
void _iftCalcAxisStride
(const iftSICLE *sicle, float *xstride, float *ystride)
{
  #ifdef IFT_DEBUG //--------------------------------------------------------//
  assert(sicle != NULL);
  assert(xstride != NULL);
  assert(ystride != NULL);
  #endif //------------------------------------------------------------------//
  int total;
  float perc_x, perc_y, c;

  total = sicle->mimg->xsize + sicle->mimg->ysize;
  perc_x = sicle->mimg->xsize/(float)total;
  perc_y = sicle->mimg->ysize/(float)total;

  // N0 = Kx * Ky == (Px * C)(Py * C) == (PxPy)C^2 for 2D
  c = (float)sqrtf(sicle->n0/(float)(perc_x * perc_y));

  // Distribute seeds proportionally to the axis' size
  (*xstride) = sicle->mimg->xsize/(float)(perc_x * c);
  (*ystride) = sicle->mimg->ysize/(float)(perc_y * c);
}

/*
  Samples the desired initial number of seeds by a grid scheme selection. If a
  candidate is within a forbidden location, it is not resampled.
*/
iftIntArray *_iftRunGridSampl
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  int x0, xf, y0, yf;
  float xstride, ystride;
  iftSet *tmp_seeds;
  iftIntArray *seeds;

  // Calculate each axis' stride for equal distribution of seeds
  _iftCalcAxisStride(sicle, &xstride, &ystride);

  if(xstride < 1.0 || ystride < 1.0) // If jump size may fall in the same spel
    iftError("Excessive number of seeds!", "_iftRunGridSampl");

  x0 = (int)(xstride/2.0); xf = sicle->mimg->xsize - 1;
  y0 = (int)(ystride/2.0); yf = sicle->mimg->ysize - 1;

  tmp_seeds = NULL; // Temporary storage for agglomerating the seed indexes
  for(int y = y0; y <= yf; y = (int)(y + ystride)) // For each y-index
  {
    for(int x = x0; x <= xf; x = (int)(x + xstride)) // For each x-index
    {
      int curr_idx;
      iftVoxel curr_vxl;

      curr_vxl.x = x; curr_vxl.y = y; curr_vxl.z = 0; // Create the seed pixel

      curr_idx = iftMGetVoxelIndex(sicle->mimg, curr_vxl); // Get its index

      // If a mask was not provided OR it is not a forbidden spel
      if(sicle->mask == NULL || iftBMapValue(sicle->mask, curr_idx) == true)
        iftInsertSet(&tmp_seeds, curr_idx); // Add it as a seed
    }
  }

  seeds = iftSetToArray(tmp_seeds); // Convert to an array for simplicity
  iftDestroySet(&tmp_seeds);

  return seeds;
}

/*
  Samples the desired initial number of seeds by random selection. If a
  candidate is within a forbidden location, it is resampled.
*/
iftIntArray *_iftRunRndSampl
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  int num_sampled;
  iftBMap *marked;
  iftIntArray *seeds;

  seeds = iftCreateIntArray(sicle->n0); // Assured N0 seeds to be sampled
  marked = iftCreateBMap(sicle->mimg->n); // Marked indexes

  num_sampled = 0;
  while(num_sampled < sicle->n0) // Until the desired quantity is reached
  {
    int idx;

    idx = iftRandomInteger(0, sicle->mimg->n - 1);

    // If (a mask was not provided OR it is not a forbidden spel) AND
    //    it was not selected as seed in a previous iteration
    if((sicle->mask == NULL || iftBMapValue(sicle->mask, idx) == true) &&
       iftBMapValue(marked, idx) == false)
    {
      seeds->val[num_sampled] = idx;
      iftBMapSet1(marked, idx);
      num_sampled++;
    }
    // Then, such seed must be resampled to another position
  }

  iftDestroyBMap(&marked);

  return seeds;
}

/*
  Samples the initial seed set in accordance to the configuration given in
  parameter.
*/
iftIntArray* _iftSampleSeeds
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  iftIntArray *seeds;

  seeds = NULL;
  // Sample the initial seed set
  if(sicle->sampl_opt == IFT_SICLE_SAMPL_GRID)
    seeds = _iftRunGridSampl(sicle);
  else if(sicle->sampl_opt == IFT_SICLE_SAMPL_RND)
    seeds = _iftRunRndSampl(sicle);
  else
    iftError("Unknown seed sampling option", "_iftSampleSeeds");

  return seeds;
}

//===========================================================================//
// IMAGE FORESTING TRANSFORM
//===========================================================================//
/*
  Gets and/or initializes the seeds' feature vector with respect to the
  arc-cost estimation defined by the user.
*/
float **_iftGetSeedsFeats
(const iftSICLE *sicle, const iftIntArray *seeds, const iftIntArray *old_seeds,
  const _iftIFTData *iftdata)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(seeds != NULL);
  #endif //------------------------------------------------------------------//
  float **seed_feats;

  seed_feats = calloc(seeds->n, sizeof(float*));
  assert(seed_feats != NULL);

  for(long i = 0; i < seeds->n; ++i) // For each seed/tree
  {
    seed_feats[i] = calloc(sicle->mimg->m, sizeof(float));
    assert(seed_feats[i] != NULL);

    for(int j = 0; j < sicle->mimg->m; ++j) // For each feature
      if(sicle->arc_opt == IFT_SICLE_ARCCOST_ROOT)
        seed_feats[i][j] = sicle->mimg->val[seeds->val[i]][j]; //Copy its values
      // for dynamic estimation, there's nothing to do here
  }

  return seed_feats;
}

/*
  Runs the seed-restricted IFT considering the options defined in the
  parameter. Note that the IFT data is updated inplace.
*/
void _iftRunIFT
(const iftSICLE *sicle, const iftIntArray *seeds, const iftIntArray *old_seeds,
 _iftIFTData **iftdata)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(seeds != NULL);
  assert(iftdata != NULL);
  assert(*iftdata != NULL);
  #endif //------------------------------------------------------------------//
  int alpha;
  int *tree_size;
  float **seed_feats;
  iftAdjRel *A;
  iftDHeap *heap;
	
  // If the user desires to boost delineation for passive environments 
  alpha = (sicle->enable_boost) ? 1 : 0;

  // If it is permitted to consider the diagonal neighbors
  if(sicle->use_diag_adj == true) A = iftCircular(sqrtf(2.0));
  else A = iftCircular(1.0);

  tree_size = calloc(seeds->n, sizeof(int));
  assert(tree_size != NULL);

  // Get the seed features based on the user's selection of arc estimation
  seed_feats = _iftGetSeedsFeats(sicle, seeds, old_seeds, *iftdata);

  _iftResetIFTData(sicle, iftdata); // Clear the IFT data for a new iteration

  heap = iftCreateDHeap(sicle->mimg->n, (*iftdata)->cost_map);
  iftSetRemovalPolicyDHeap(heap, MINVALUE);

  // Place the avaliable seeds in the mappings
  for(long i = 0; i < seeds->n; ++i) // For each seed
  {
    int idx;

    idx = seeds->val[i];
    (*iftdata)->cost_map[idx] = 0; // Assuring that seeds are 1st in line
    (*iftdata)->id_map[idx] = i; // A seed maps to its own id
    iftInsertDHeap(heap, idx);
  }

  while(iftEmptyDHeap(heap) == false) // While exists a pixel to be evaluated
  {
    int p_idx, seed_id;
    float seed_saliency;
    float *feats;
    iftVoxel p_vxl;

    p_idx = iftRemoveDHeap(heap);
    p_vxl = iftMGetVoxelCoord(sicle->mimg, p_idx);

    // Get seed id and saliency (when provided)
    seed_id = (*iftdata)->id_map[p_idx];
	if(sicle->saliency == NULL) seed_saliency = 0;
	else seed_saliency = sicle->saliency[seeds->val[seed_id]];

    // Update the seed features
    if(sicle->arc_opt == IFT_SICLE_ARCCOST_DYN) // If it is an dynamic arc cost
    {
      for(int j = 0; j < sicle->mimg->m; ++j) // For each feature
      {
        seed_feats[seed_id][j] *= tree_size[seed_id]; // Get accumulated value
        seed_feats[seed_id][j] += sicle->mimg->val[p_idx][j]; // Add features
        seed_feats[seed_id][j] /= (float)(tree_size[seed_id] + 1); // Get mean
      }
    } // Then, nothing to do for root-based arc costs
    tree_size[seed_id]++; // Increase tree's size

    feats = seed_feats[seed_id]; // For readability

    for(int i = 1; i < A->n; ++i) // For each adjacent pixel
    {
      iftVoxel adj_vxl;

      adj_vxl = iftGetAdjacentVoxel(A, p_vxl, i);

      // If it is within the image domain
      if(iftMValidVoxel(sicle->mimg, adj_vxl) == true)
      {
        int adj_idx;

        adj_idx = iftMGetVoxelIndex(sicle->mimg, adj_vxl);

        if(heap->color[adj_idx] != IFT_BLACK) // If it wasn't ordely removed
        {
          double arccost, pathcost;
		  float adj_saliency,tmp;
          float *adj_feats;
		  
		  // Get adjacent features and saliency (when provided)
          adj_feats = sicle->mimg->val[adj_idx];
		  if(sicle->saliency == NULL) adj_saliency = 0;
		  else adj_saliency = sicle->saliency[adj_idx];
		  
          arccost = _iftEuclDist(feats, adj_feats, sicle->mimg->m); // L2-Norm
		  tmp = fabs(seed_saliency - adj_saliency);
          arccost = pow(arccost,1 + alpha * tmp);

          pathcost = iftMax((*iftdata)->cost_map[p_idx], arccost); // Fmax
          // Does the new path offer a lesser accumulated cost?
          if(pathcost < (*iftdata)->cost_map[adj_idx])
          {
            if (heap->color[adj_idx] == IFT_GRAY) // If it is within the heap
              iftRemoveDHeapElem(heap, adj_idx); // Remove for updating

            (*iftdata)->id_map[adj_idx] = seed_id;    // Mark as conquered by
            (*iftdata)->pred_map[adj_idx] = p_idx;    // the current pixel
            (*iftdata)->cost_map[adj_idx] = pathcost; // evaluated

            // Insert to be evaluated in the future
            iftInsertDHeap(heap, adj_idx);
          }
        }
      }
    }
  }

  for(long i = 0; i < seeds->n; ++i) // For each seed
    free(seed_feats[i]);
  free(seed_feats);
  free(tree_size);

  iftDestroyDHeap(&heap);
  iftDestroyAdjRel(&A);
}

//===========================================================================//
// SEED PRESERVATION CURVE
//===========================================================================//
/*
  Calculates the necessary number of iterations for reaching the desired
  number of superpixels, with respect to the maximum number of iterations
  given, and assuming that the final number of superpixels is always greater
  than 1.
*/
inline int _iftCalcNumIters
(const iftSICLE *sicle, const int real_n0)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  const float MIN_NF = 1.0; // Nf should always be > 1
  float n0, nf, max_iters;  // For readability
  float base, exp; // Base and its exponent
  float approx;
  int num_iters;

  n0 = real_n0;
  nf = (float)sicle->nf;
  max_iters = (float)sicle->max_iters;

  approx = 0;

  // f(i) = N0 * base^(-i)
  if(sicle->num_scales == 0)
  {
    exp = 1.0/(max_iters - 1);
  	base = (float)pow(n0/MIN_NF, exp);
  	approx = iftLog(n0/nf, base);
  	num_iters = ceil(approx) + 1; // Approximate number + the final one (Ni = Nf)
  }
  else num_iters = sicle->num_scales + 1; // Scales + Nf

  return num_iters;
}

/*
  Calculates the quantity of irrelevant seeds to be removed in the given
  iteration, considering the configuration given in parameter.
*/
inline int _iftCalcNumToRem
(const iftSICLE *sicle, const int real_n0, const int iter)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(iter >= 0);
  #endif //------------------------------------------------------------------//
  const float MIN_NF = 1.0; // Nf should always be > 1
  float n0, nf, max_iters; // For readability
  float base, exp; // Base and its exponent
  float perc; // Percentage of N0 seeds to be selected
  int ni;

  n0 = real_n0;
  nf = (float)sicle->nf;
  max_iters = (float)sicle->max_iters;

  perc = 0;

  // f(i) = N0 * base^(-i)
  if(sicle->num_scales == 0)
  {
    exp = 1.0/(max_iters - 1);
    base = (float)pow(n0/MIN_NF, exp);
    perc = pow(base, -iter);
    ni = iftMax(iftRound(n0 * perc), nf); // Curve "cutting"
  }
  else ni = sicle->scales[iter - 1];

  return ni;
}


//===========================================================================//
// SEED REMOVAL
//===========================================================================//
/*
  Assigns a random priority to each seed existing in the given set.
*/
double *_iftCalcRndSeedPrio
(const iftIntArray *seeds)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(seeds != NULL);
  #endif //------------------------------------------------------------------//
  double *prio;

  prio = calloc(seeds->n, sizeof(double));

  for(long i = 0; i < seeds->n; ++i) // For each seed
    prio[i] = iftRandomUniform(0.0, 1.0); // Set a random priority

  return prio;
}

/*
  Assigns a priority value to each seed based on its properties and based on
  its contrast to its adjacents (considering either saliency or color). If a
  saliency map is provided, then all options become object-based.
*/
double *_iftCalcSeedPrio
(const iftSICLE *sicle, const _iftIFTData *iftdata, const iftIntArray *seeds)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(seeds != NULL);
  assert(iftdata != NULL);
  #endif //------------------------------------------------------------------//
  int alpha;
  double *prio;
  _iftForestStats *forstats;

  prio = calloc(seeds->n, sizeof(double));
  // Calculate the statistics of the generated forest
  forstats = _iftCalcForestStats(sicle, iftdata, seeds);

  alpha = (sicle->enable_boost) ? 1 : 0;

  #if IFT_OMP //-------------------------------------------------------------//
  #pragma omp parallel for
  #endif //------------------------------------------------------------------//
  for(long i = 0; i < seeds->n; ++i)
  {
    double size_perc, min_color_grad, max_color_grad, max_sal_grad;

    // Size percentage with respect to the image's size (in pixels)
    size_perc = forstats->tree_size[i]/(float)sicle->mimg->n;

    max_color_grad = max_sal_grad = 0.0;
    min_color_grad = IFT_INFINITY_DBL;
    for(int j = 0; j < forstats->num_trees; ++j) // For each seed/tree
    {
      double sal_grad, grad;

      if(iftBMapValue(forstats->tree_adj[i], j) == true)//If they are adjacents
      {
        // L2-Norm
        grad = _iftEuclDist(forstats->tree_feats[i], forstats->tree_feats[j],
                            forstats->num_feats);
        // L1-Norm
        sal_grad = fabs(forstats->tree_sal[i] - forstats->tree_sal[j]);

        // Update values
        if(grad < min_color_grad) min_color_grad = grad;
        if(grad > max_color_grad) max_color_grad = grad;
        if(sal_grad > max_sal_grad) max_sal_grad = sal_grad;
      } // Then, do not consider in the computation
    }

    if(sicle->rem_opt == IFT_SICLE_REM_MINCONTR)
      prio[i] = min_color_grad;
    else if(sicle->rem_opt == IFT_SICLE_REM_MAXCONTR)
      prio[i] = max_color_grad;
    else if(sicle->rem_opt == IFT_SICLE_REM_SIZE)
      prio[i] = size_perc;
    else if(sicle->rem_opt == IFT_SICLE_REM_MAXSC)
      prio[i] = size_perc * max_color_grad;
    else if(sicle->rem_opt == IFT_SICLE_REM_MINSC)
      prio[i] = size_perc * min_color_grad;

    prio[i] *= iftMax((1-alpha) * forstats->tree_sal[i], max_sal_grad);
  }
  _iftDestroyForestStats(&forstats);

  return prio;
}

/*
  Selects the desired quantity of seeds with highest relevance --- with
  respect to the given criterion --- and returns them.
*/
iftIntArray *_iftRemSeeds
(const iftSICLE *sicle, const int ni, const _iftIFTData *iftdata,
 const iftIntArray *seeds)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(seeds != NULL);
  assert(ni > 1 && ni <= seeds->n);
  #endif //------------------------------------------------------------------//
  double *prio;
  iftDHeap *heap;
  iftIntArray *new_seeds;

  prio = NULL;
  // Compute the seed priorities
  if(sicle->rem_opt == IFT_SICLE_REM_RND)
    prio = _iftCalcRndSeedPrio(seeds);
  else if(sicle->rem_opt == IFT_SICLE_REM_MINCONTR ||
          sicle->rem_opt == IFT_SICLE_REM_MAXCONTR ||
          sicle->rem_opt == IFT_SICLE_REM_SIZE ||
          sicle->rem_opt == IFT_SICLE_REM_MAXSC ||
          sicle->rem_opt == IFT_SICLE_REM_MINSC )
    prio = _iftCalcSeedPrio(sicle, iftdata, seeds);
  else
    iftError("Unknown seed removal criterion", "_iftRemSeeds");

  heap = iftCreateDHeap(seeds->n, prio);
  iftSetRemovalPolicyDHeap(heap, MAXVALUE); // Select the most relevants

  for(long i = 0; i < seeds->n; ++i) // For each seed
    iftInsertDHeap(heap, i); // Insert with its respective priority

  new_seeds = iftCreateIntArray(ni); // Assured to select ni seeds
  for(long i = 0; i < ni; ++i) // For each seed selected
  {
    int seed_id;

    seed_id = iftRemoveDHeap(heap); // Remove the one with highest priority
    new_seeds->val[i] = seeds->val[seed_id]; // Add its index to the array
  } // The rest are discarted

  free(prio);
  iftDestroyDHeap(&heap);

  return new_seeds;
}

/*****************************************************************************\
*
*                               PUBLIC FUNCTIONS
*
\*****************************************************************************/
//===========================================================================//
// CONSTRUCTOR & DESTRUCTOR
//===========================================================================//
iftSICLE *iftCreateSICLE
(const iftImage *img, const iftImage *mask, const iftImage *objsm)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(img != NULL);
  if(mask != NULL) iftVerifyImageDomains(img, mask, "iftCreateSICLE");
  if(objsm != NULL) iftVerifyImageDomains(img, objsm, "iftCreateSICLE");
  #endif //------------------------------------------------------------------//
  iftSICLE *sicle;

  sicle = malloc(sizeof(*sicle));
  assert(sicle != NULL);

  if(iftIsColorImage(img) == true) // If it is a colored image or video
    sicle->mimg = iftImageToMImage(img, LAB_CSPACE); // Set to CIELAB
  else // Then, deal with the luminosity values
    sicle->mimg = iftImageToMImage(img, GRAY_CSPACE);

  if(mask != NULL) // If a mask was provided
    sicle->mask = iftBinImageToBMap(mask); // Saves a lot of memory
  else sicle->mask = NULL; // Then, all spels are to be conquered

  if(objsm != NULL) // If a saliency map was provided
  {
    int max_sal, min_sal;

    iftMinMaxValues(objsm, &min_sal, &max_sal); // For normalizing the values

    sicle->saliency = calloc(objsm->n, sizeof(float));
    assert(sicle->saliency != NULL);

    #if IFT_OMP //-----------------------------------------------------------//
    #pragma omp parallel for
    #endif //----------------------------------------------------------------//
    for(int p = 0; p < objsm->n; ++p) // For each spel
    {
      sicle->saliency[p] = (objsm->val[p] - min_sal);   // Normalize saliency
      sicle->saliency[p] /= (float)(max_sal - min_sal); // to [0,1]
    }
  }
  else sicle->saliency = NULL; // Then, no object information is considered

  // Default values
  sicle->n0 = 3000;
  sicle->nf = 200;
  sicle->max_iters = 5;
  sicle->use_diag_adj = true;
  sicle->num_scales = 0;
  sicle->scales = NULL;
  sicle->enable_boost = false;
  // Default configuration
  sicle->sampl_opt = IFT_SICLE_SAMPL_RND;
  sicle->arc_opt = IFT_SICLE_ARCCOST_ROOT;
  sicle->rem_opt = IFT_SICLE_REM_MINSC;

  return sicle;
}

void iftDestroySICLE
(iftSICLE **sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  #endif //------------------------------------------------------------------//
  iftDestroyMImage(&((*sicle)->mimg));
  if((*sicle)->mask != NULL) iftDestroyBMap(&((*sicle)->mask));
  if((*sicle)->saliency != NULL) free((*sicle)->saliency);
  if((*sicle)->num_scales > 0) free((*sicle)->scales);

  free(*sicle);
  (*sicle) = NULL;
}

//===========================================================================//
// GETTERS
//===========================================================================//
inline int iftSICLEGetMaxIters
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  return sicle->max_iters;
}

inline int iftSICLEGetN0
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  return sicle->n0;
}

inline int iftSICLEGetNf
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  return sicle->nf;
}

inline int iftSICLEGetNumScales
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  return sicle->num_scales;
}

inline int *iftSICLEGetScales
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  return sicle->scales;
}

//===========================================================================//
// SETTERS
//===========================================================================//
inline void iftSICLESetMaxIters
(iftSICLE **sicle, const int max_iters)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  assert(max_iters > 0);
  #endif //------------------------------------------------------------------//
  (*sicle)->max_iters = max_iters;
}

inline void iftSICLESetN0
(iftSICLE **sicle, const int n0)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  assert(n0 >= 2 && n0 >= (*sicle)->nf);
  #endif //------------------------------------------------------------------//
  (*sicle)->n0 = n0;
}

inline void iftSICLESetNf
(iftSICLE **sicle, const int nf)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  assert(nf >= 2 && nf <= (*sicle)->n0);
  #endif //------------------------------------------------------------------//
  (*sicle)->nf = nf;
}

void iftSICLESetScales
(iftSICLE **sicle, const int num_scales, const int* scales)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  assert(num_scales >= 0);
  #endif //------------------------------------------------------------------//
  (*sicle)->num_scales = num_scales + 1;
  if(num_scales > 0)
  {
	(*sicle)->scales = calloc(num_scales + 1, sizeof(int));
	for(int i = 0; i < num_scales; ++i) (*sicle)->scales[i] = scales[i];
	(*sicle)->scales[num_scales] = (*sicle)->nf;
  }
}

inline void iftSICLEUseDiagAdj
(iftSICLE **sicle, const bool use)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  #endif //------------------------------------------------------------------//
  (*sicle)->use_diag_adj = use;
}

inline void iftSICLEEnableBoost
(iftSICLE **sicle, const bool enable)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  #endif //------------------------------------------------------------------//
  (*sicle)->enable_boost = enable;
}

//===========================================================================//
// SEED SAMPLING
//===========================================================================//
inline iftSICLESampl iftSICLEGetSamplOpt
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  return sicle->sampl_opt;
}

inline void iftSICLESetSamplOpt
(iftSICLE **sicle, const iftSICLESampl sampl_opt)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  #endif //------------------------------------------------------------------//
  (*sicle)->sampl_opt = sampl_opt;
}

//===========================================================================//
// ARC COST FUNCTION
//===========================================================================//
inline iftSICLEArcCost iftSICLEGetArcCostOpt
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  return sicle->arc_opt;
}

inline void iftSICLESetArcCostOpt
(iftSICLE **sicle, const iftSICLEArcCost arc_opt)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  #endif //------------------------------------------------------------------//
  (*sicle)->arc_opt = arc_opt;
}

//===========================================================================//
// SEED REMOVAL
//===========================================================================//
inline iftSICLERem iftSICLEGetRemOpt
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  return sicle->rem_opt;
}

inline void iftSICLESetRemOpt
(iftSICLE **sicle, const iftSICLERem rem_opt)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  #endif //------------------------------------------------------------------//
  (*sicle)->rem_opt = rem_opt;
}

//===========================================================================//
// VERIFIERS
//===========================================================================//
inline bool iftSICLEUsingDiagAdj
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  return sicle->use_diag_adj;
}

inline bool iftSICLEBoostEnabled
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  return sicle->enable_boost;
}
//===========================================================================//
// RUNNER
//===========================================================================//
iftImage **iftRunSICLE
(const iftSICLE *sicle)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  #endif //------------------------------------------------------------------//
  int num_iters, real_n0;
  iftIntArray *seeds;
  iftImage **label_img;
  _iftIFTData *iftdata;

  seeds = _iftSampleSeeds(sicle); // Sample the initial seed set
  // Since the latter may not guarantee N0, we need to consider the REAL N0
  // over the one desired by the user
  real_n0 = (int)seeds->n;
  #if IFT_DEBUG //---------------------------------------------------------//
  fprintf(stderr, "iftRunSICLE: Number of seeds = %ld\n", seeds->n);
  #endif //----------------------------------------------------------------//

  iftdata = _iftCreateIFTData(sicle->mimg->n); // Create an empty IFT data

  // Get the necessary number of iterations for segmentation
  num_iters = _iftCalcNumIters(sicle, real_n0);
  if(num_iters <= 1) // At least two iterations MUST be performed
    iftError("The number of seeds is too low", "iftRunSICLE");

  label_img = calloc(num_iters, sizeof(iftImage*));

  if(label_img == NULL)
    iftError("Could not alloc array of Images","iftRunSICLE");

  for(int iter = 1; iter <= num_iters; ++iter) // For each iteration
  {
    int ni;

    #if IFT_DEBUG //---------------------------------------------------------//
    fprintf(stderr, "iftRunSICLE: Iteration %d of %d\n", iter, num_iters);
    #endif //----------------------------------------------------------------//
    // Segment image using the IFT
    _iftRunIFT(sicle, seeds, seeds, &iftdata);

    // Creates a label image based on the seed id's in the IFT's id map
	if(sicle->num_scales > 0 && iter < num_iters)
      label_img[num_iters - iter] = _iftCreateLabelImageFromIdMap(sicle, seeds, iftdata);

    if(iter < num_iters) // If it is not the last iteration
    {
	  iftIntArray *old_seeds;

      // Calculate the number of seeds to be selected (for the next iteration)
      ni = _iftCalcNumToRem(sicle, real_n0, iter);
      #if IFT_DEBUG //-------------------------------------------------------//
      fprintf(stderr, "iftRunSICLE: Ni = %d\n", ni);
      #endif //--------------------------------------------------------------//
      old_seeds = seeds;
      // Select the Ni relevant seeds (and remove the others)
      seeds = _iftRemSeeds(sicle, ni, iftdata, old_seeds);
  	  iftDestroyIntArray(&old_seeds);
    }
  }
  label_img[0] = _iftCreateLabelImageFromIdMap(sicle, seeds, iftdata);

  _iftDestroyIFTData(&iftdata);
  iftDestroyIntArray(&seeds);

  return label_img;
}
