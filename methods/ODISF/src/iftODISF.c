/*****************************************************************************\
* iftODISF.c
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-06-15
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#include "iftODISF.h"

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
  int num_spels; // Number of spels within the image/video
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
  int *tree_size; // Size of each tree (in spels)
  float *tree_sal; // Mean saliency value of each tree
  float **tree_feats; // Mean feature vector of each tree
  iftBMap **tree_adj; // Tree adjacency relation
} _iftForestStats;

/*****************************************************************************\
*
*                               PUBLIC STRUCTS
*
\*****************************************************************************/
struct ift_odisf_alg
{
  bool use_diag_adj; // Use diagonal adjacents (i.e., 8-adjacency)?
  int n0, nf; // Initial number of seeds and final quantity of superspels
  float *saliency; // Normalized object saliency values (if they are provided)
  iftMImage *mimg; // Multiband image/video
  iftBMap *mask; // Mask indicating the ROI (if it exists)
  iftODISFSampl sampl_opt; // Seed sampling option
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
(const iftODISF *odisf, const iftIntArray *seeds, const _iftIFTData *iftdata)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(seeds != NULL);
  assert(iftdata != NULL);
  #endif //------------------------------------------------------------------//
  int label;
  iftImage *label_img;

  label_img = iftCreateImage(odisf->mimg->xsize, odisf->mimg->ysize, 
                             odisf->mimg->zsize);
    
  label = 1; // 0 is for the background
  for(long i = 0; i < seeds->n; ++i)
    label_img->val[seeds->val[i]] = label++; // Assign unique label

  #if IFT_OMP //-------------------------------------------------------------//
  #pragma omp parallel for
  #endif //------------------------------------------------------------------//
  for(int p = 0; p < label_img->n; ++p) // For each spel
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
  by the number of spels given
*/
_iftIFTData *_iftCreateIFTData
(const int num_spels)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(num_spels > 0);
  #endif //------------------------------------------------------------------//
  _iftIFTData *iftdata;

  iftdata = malloc(sizeof(_iftIFTData));
  assert(iftdata != NULL);

  iftdata->num_spels = num_spels;
  iftdata->id_map = calloc(num_spels, sizeof(int));
  assert(iftdata->id_map != NULL);
  iftdata->pred_map = calloc(num_spels, sizeof(int));
  assert(iftdata->pred_map != NULL);
  iftdata->cost_map = calloc(num_spels, sizeof(double));
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
(const iftODISF *odisf, _iftIFTData **iftdata)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(iftdata != NULL);
  assert(*iftdata != NULL);
  #endif //------------------------------------------------------------------//
  #if IFT_OMP //-------------------------------------------------------------//
  #pragma omp parallel for
  #endif //------------------------------------------------------------------//
  for(int p = 0; p < odisf->mimg->n; ++p) // For each spel
  {
    (*iftdata)->id_map[p] = IFT_NIL;   // If p is not in the background, such
    (*iftdata)->pred_map[p] = IFT_NIL; // values will surely be changed

    // If a mask was not provided OR it is not a forbidden spel
    if(odisf->mask == NULL || iftBMapValue(odisf->mask, p) == true)
      (*iftdata)->cost_map[p] = IFT_INFINITY_DBL; // Exposed to conquering
    else // Then, it is forbidden
      (*iftdata)->cost_map[p] = IFT_INFINITY_DBL_NEG; // Safe from conquering
  }
}

//===========================================================================//
// FOREST STATISTICS
//===========================================================================//
/*
  Creates an empty IFT forest statistics object based on the IFT data generated
  in the segmentation and its inputs.
*/
_iftForestStats *_iftCreateForestStats
(const iftODISF *odisf, const _iftIFTData *iftdata, const iftIntArray *seeds)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(iftdata != NULL);
  assert(seeds != NULL);
  #endif //------------------------------------------------------------------//
  _iftForestStats *forstats;

  forstats = malloc(sizeof(_iftForestStats));
  assert(forstats != NULL);

  forstats->num_trees = (int)seeds->n;
  forstats->num_feats = odisf->mimg->m;
  
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
(const iftODISF *odisf, const _iftIFTData *iftdata, const iftIntArray *seeds)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(iftdata != NULL);
  assert(seeds != NULL);
  #endif //------------------------------------------------------------------//
  _iftForestStats *forstats;
  iftAdjRel *A;

  forstats = _iftCreateForestStats(odisf, iftdata, seeds);

  // Use the same adjacency as the one used in segmentation
  if(odisf->use_diag_adj == true) A = iftCircular(sqrtf(2.0));
  else A = iftCircular(1.0);

  for(int p = 0; p < odisf->mimg->n; ++p) // For each spel
  {
    int seed_id;

    seed_id = iftdata->id_map[p];
    if(seed_id != IFT_NIL) // If it is not in the background
    {
      iftVoxel p_vxl;

      p_vxl = iftMGetVoxelCoord(odisf->mimg, p);

      forstats->tree_size[seed_id]++;

      if(odisf->saliency != NULL) // If a saliency map was provided
        forstats->tree_sal[seed_id] += odisf->saliency[p];

      for(int j = 0; j < odisf->mimg->m; ++j) // For each feature
        forstats->tree_feats[seed_id][j] += odisf->mimg->val[p][j];

      for(int i = 1; i < A->n; ++i) // For each adjacent
      {
        iftVoxel adj_vxl;

        adj_vxl = iftGetAdjacentVoxel(A, p_vxl, i);

        // If it is within the image domain
        if(iftMValidVoxel(odisf->mimg, adj_vxl) == true)
        {
          int adj_idx, adj_seed_id;

          adj_idx = iftMGetVoxelIndex(odisf->mimg, adj_vxl);
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
    if(odisf->saliency != NULL) // If a saliency map was provided
      // Compute mean saliency
      forstats->tree_sal[i] /= (float)forstats->tree_size[i];
    else forstats->tree_sal[i] = 1.0; // Then, all trees are relevant
    
    for(int j = 0; j < odisf->mimg->m; ++j) // For each feature
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
  Samples the desired initial number of seeds by a grid scheme selection. If a 
  candidate is within a forbidden location, it is not resampled.
*/
iftIntArray *_iftRunGridSampl
(const iftODISF *odisf)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  #endif //------------------------------------------------------------------//
  int x0, xf, y0, yf;
  float xstride, ystride;
  iftSet *tmp_seeds;
  iftIntArray *seeds;

  // Calculate each axis' stride for equal distribution of seeds
  xstride = ystride = sqrtf(odisf->mimg->n/(float)odisf->n0);
  
  if(xstride < 1.0 || ystride < 1.0) // If jump size may fall in the same spel
    iftError("Excessive number of seeds!", "_iftRunGridSampl");

  x0 = (int)(xstride/2.0); xf = odisf->mimg->xsize - 1;
  y0 = (int)(ystride/2.0); yf = odisf->mimg->ysize - 1;

  tmp_seeds = NULL; // Temporary storage for agglomerating the seed indexes
  for(int y = y0; y <= yf; y = (int)(y + ystride)) // For each y-index
  {
    for(int x = x0; x <= xf; x = (int)(x + xstride)) // For each x-index
    {
      int curr_idx;
      iftVoxel curr_vxl;

      curr_vxl.x = x; curr_vxl.y = y; curr_vxl.z = 0; // Create the seed voxel

      curr_idx = iftMGetVoxelIndex(odisf->mimg, curr_vxl); // Get its index

      // If a mask was not provided OR it is not a forbidden spel
      if(odisf->mask == NULL || iftBMapValue(odisf->mask, curr_idx) == true)
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
(const iftODISF *odisf)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  #endif //------------------------------------------------------------------//
  int num_sampled;
  iftBMap *marked;
  iftIntArray *seeds;

  seeds = iftCreateIntArray(odisf->n0); // Assured N0 seeds to be sampled
  marked = iftCreateBMap(odisf->mimg->n); // Marked indexes

  num_sampled = 0;
  while(num_sampled < odisf->n0) // Until the desired quantity is reached
  {
    int idx;

    idx = iftRandomInteger(0, odisf->mimg->n - 1);

    // If (a mask was not provided OR it is not a forbidden spel) AND
    //    it was not selected as seed in a previous iteration
    if((odisf->mask == NULL || iftBMapValue(odisf->mask, idx) == true) &&
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
(const iftODISF *odisf)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  #endif //------------------------------------------------------------------//
  iftIntArray *seeds;

  seeds = NULL;
  // Sample the initial seed set
  if(odisf->sampl_opt == IFT_ODISF_SAMPL_GRID)
    seeds = _iftRunGridSampl(odisf); 
  else if(odisf->sampl_opt == IFT_ODISF_SAMPL_RND)
    seeds = _iftRunRndSampl(odisf);
  else
    iftError("Unknown seed sampling option", "_iftSampleSeeds");

  return seeds;
}

//===========================================================================//
// IMAGE FORESTING TRANSFORM
//===========================================================================//
/*
  Runs the seed-restricted IFT considering the options defined in the 
  parameter. Note that the IFT data is updated inplace.
*/
void _iftRunIFT
(const iftODISF *odisf, const iftIntArray *seeds, _iftIFTData **iftdata)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(seeds != NULL);
  assert(iftdata != NULL);
  assert(*iftdata != NULL);
  #endif //------------------------------------------------------------------//
  int *tree_size;
  float **seed_feats;
  iftAdjRel *A;
  iftDHeap *heap;

  // If it is permitted to consider the diagonal neighbors
  if(odisf->use_diag_adj == true) A = iftCircular(sqrtf(2.0));
  else A = iftCircular(1.0);

  tree_size = calloc(seeds->n, sizeof(int));
  assert(tree_size != NULL);

  // Get the seed features based on the user's selection of arc estimation
  seed_feats = calloc(seeds->n, sizeof(float*));
  assert(seed_feats != NULL);

  for(long i = 0; i < seeds->n; ++i) // For each seed/tree
  {
    seed_feats[i] = calloc(odisf->mimg->m, sizeof(float));
    assert(seed_feats[i] != NULL);
  }
  
  _iftResetIFTData(odisf, iftdata); // Clear the IFT data for a new iteration

  heap = iftCreateDHeap(odisf->mimg->n, (*iftdata)->cost_map);
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

  while(iftEmptyDHeap(heap) == false) // While exists a spel to be evaluated
  {
    int p_idx, seed_id;
    float *feats;
    iftVoxel p_vxl;

    p_idx = iftRemoveDHeap(heap);
    p_vxl = iftMGetVoxelCoord(odisf->mimg, p_idx);
    
    seed_id = (*iftdata)->id_map[p_idx];

    // Update the seed features 
    for(int j = 0; j < odisf->mimg->m; ++j) // For each feature
    {
      seed_feats[seed_id][j] *= tree_size[seed_id]; // Get accumulated value
      seed_feats[seed_id][j] += odisf->mimg->val[p_idx][j]; // Add features
      seed_feats[seed_id][j] /= (float)(tree_size[seed_id] + 1); // Get mean
    }
    tree_size[seed_id]++; // Increase tree's size
    
    feats = seed_feats[seed_id]; // For readability

    for(int i = 1; i < A->n; ++i) // For each adjacent spel
    {
      iftVoxel adj_vxl;

      adj_vxl = iftGetAdjacentVoxel(A, p_vxl, i);

      // If it is within the image domain
      if(iftMValidVoxel(odisf->mimg, adj_vxl) == true)
      {
        int adj_idx;

        adj_idx = iftMGetVoxelIndex(odisf->mimg, adj_vxl);

        if(heap->color[adj_idx] != IFT_BLACK) // If it wasn't ordely removed
        {
          double arccost, pathcost;
          float *adj_feats;

          adj_feats = odisf->mimg->val[adj_idx];
          
          arccost = _iftEuclDist(feats, adj_feats, odisf->mimg->m); // L2-Norm
          pathcost = iftMax((*iftdata)->cost_map[p_idx], arccost); // Fmax

          // Does the new path offer a lesser accumulated cost?
          if(pathcost < (*iftdata)->cost_map[adj_idx])
          {
            if (heap->color[adj_idx] == IFT_GRAY) // If it is within the heap
              iftRemoveDHeapElem(heap, adj_idx); // Remove for updating

            (*iftdata)->id_map[adj_idx] = seed_id;    // Mark as conquered by
            (*iftdata)->pred_map[adj_idx] = p_idx;    // the current spel
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
// SEED REMOVAL
//===========================================================================//
/*
  Assigns a priority value to each seed based on its properties and based on
  its contrast to its adjacents (considering either saliency or color). If a
  saliency map is provided, then it becomes object-based.
*/
double *_iftCalcSeedPrio
(const iftODISF *odisf, const _iftIFTData *iftdata, const iftIntArray *seeds)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(seeds != NULL);
  assert(iftdata != NULL);
  #endif //------------------------------------------------------------------//
  double *prio;
  _iftForestStats *forstats;

  prio = calloc(seeds->n, sizeof(double));
  // Calculate the statistics of the generated forest
  forstats = _iftCalcForestStats(odisf, iftdata, seeds);

  #if IFT_OMP //-------------------------------------------------------------//
  #pragma omp parallel for
  #endif //------------------------------------------------------------------//
  for(long i = 0; i < seeds->n; ++i)
  {
    double size_perc, min_color_grad, max_sal_grad;

    // Size percentage with respect to the image's size (in spels)
    size_perc = forstats->tree_size[i]/(float)odisf->mimg->n;

    max_sal_grad = 0.0;
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
        if(sal_grad > max_sal_grad) max_sal_grad = sal_grad;
      } // Then, do not consider in the computation
    }

    prio[i] = size_perc * min_color_grad;
    prio[i] *= iftMax(forstats->tree_sal[i], max_sal_grad);
  }
  _iftDestroyForestStats(&forstats);

  return prio;
}

/*
  Selects the desired quantity of seeds with highest relevance and returns them.
*/
iftIntArray *_iftRemSeeds
(const iftODISF *odisf, const int ni, const _iftIFTData *iftdata, 
 const iftIntArray *seeds)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(seeds != NULL);
  assert(ni > 1 && ni <= seeds->n);
  #endif //------------------------------------------------------------------//
  double *prio;
  iftDHeap *heap;
  iftIntArray *new_seeds;

  // Compute the seed priority
  prio = _iftCalcSeedPrio(odisf, iftdata, seeds);

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
iftODISF *iftCreateODISF
(const iftImage *img, const iftImage *mask, const iftImage *objsm)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(img != NULL);
  if(mask != NULL) iftVerifyImageDomains(img, mask, "iftCreateODISF");
  if(objsm != NULL) iftVerifyImageDomains(img, objsm, "iftCreateODISF");
  #endif //------------------------------------------------------------------//
  iftODISF *odisf;

  odisf = malloc(sizeof(*odisf));
  assert(odisf != NULL);

  if(iftIsColorImage(img) == true) // If it is a colored image or video
    odisf->mimg = iftImageToMImage(img, LAB_CSPACE); // Set to CIELAB
  else // Then, deal with the luminosity values
    odisf->mimg = iftImageToMImage(img, GRAY_CSPACE);
  
  if(mask != NULL) // If a mask was provided
    odisf->mask = iftBinImageToBMap(mask); // Saves a lot of memory
  else odisf->mask = NULL; // Then, all spels are to be conquered

  if(objsm != NULL) // If a saliency map was provided
  {
    int max_sal, min_sal;
    
    iftMinMaxValues(objsm, &min_sal, &max_sal); // For normalizing the values

    odisf->saliency = calloc(objsm->n, sizeof(float));
    assert(odisf->saliency != NULL);

    #if IFT_OMP //-----------------------------------------------------------//
    #pragma omp parallel for
    #endif //----------------------------------------------------------------//
    for(int p = 0; p < objsm->n; ++p) // For each spel
    {
      odisf->saliency[p] = (objsm->val[p] - min_sal);   // Normalize saliency
      odisf->saliency[p] /= (float)(max_sal - min_sal); // to [0,1]
    }
  } 
  else odisf->saliency = NULL; // Then, no object information is considered

  // Default values
  odisf->n0 = 8000;
  odisf->nf = 200;
  odisf->use_diag_adj = true;
  // Default configuration
  odisf->sampl_opt = IFT_ODISF_SAMPL_GRID;

  return odisf;
}

void iftDestroyODISF
(iftODISF **odisf)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(*odisf != NULL);
  #endif //------------------------------------------------------------------//
  iftDestroyMImage(&((*odisf)->mimg));
  if((*odisf)->mask != NULL) iftDestroyBMap(&((*odisf)->mask));
  if((*odisf)->saliency != NULL) free((*odisf)->mask);
  
  free(*odisf);
  (*odisf) = NULL;
}

//===========================================================================//
// GETTERS
//===========================================================================//
inline int iftODISFGetN0
(const iftODISF *odisf)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  #endif //------------------------------------------------------------------//
  return odisf->n0;
}

inline int iftODISFGetNf
(const iftODISF *odisf)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  #endif //------------------------------------------------------------------//
  return odisf->nf;
}

//===========================================================================//
// SETTERS
//===========================================================================//
inline void iftODISFSetN0
(iftODISF **odisf, const int n0)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(*odisf != NULL);
  assert(n0 >= 2 && n0 >= (*odisf)->nf);
  #endif //------------------------------------------------------------------//
  (*odisf)->n0 = n0;
}

inline void iftODISFSetNf
(iftODISF **odisf, const int nf)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(*odisf != NULL);
  assert(nf >= 2 && nf <= (*odisf)->n0);
  #endif //------------------------------------------------------------------//
  (*odisf)->nf = nf;
}

inline void iftODISFUseDiagAdj
(iftODISF **odisf, const bool use)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(*odisf != NULL);
  #endif //------------------------------------------------------------------//
  (*odisf)->use_diag_adj = use;
}

//===========================================================================//
// SEED SAMPLING
//===========================================================================//
inline iftODISFSampl iftODISFGetSamplOpt
(const iftODISF *odisf)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  #endif //------------------------------------------------------------------//
  return odisf->sampl_opt;
}

inline void iftODISFSetSamplOpt
(iftODISF **odisf, const iftODISFSampl sampl_opt)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(*odisf != NULL);
  #endif //------------------------------------------------------------------//
  (*odisf)->sampl_opt = sampl_opt;
}

//===========================================================================//
// VERIFIERS
//===========================================================================//
inline bool iftODISFUsingDiagAdj
(const iftODISF *odisf)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  #endif //------------------------------------------------------------------//
  return odisf->use_diag_adj; 
}

//===========================================================================//
// RUNNER
//===========================================================================//
iftImage *iftRunODISF
(const iftODISF *odisf)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  #endif //------------------------------------------------------------------//
  int iter, real_n0, ni, ni_prev;
  iftIntArray *seeds, *old_seeds;
  iftImage *label_img;
  _iftIFTData *iftdata;
  
  old_seeds = NULL;
  seeds = _iftSampleSeeds(odisf); // Sample the initial seed set

  // Since the latter may not guarantee N0, we need to consider the REAL N0 
  // over the one desired by the user
  real_n0 = ni = (int)seeds->n; 
  #if IFT_DEBUG //---------------------------------------------------------//
  fprintf(stderr, "iftRunODISF: Number of seeds = %ld\n", seeds->n);
  #endif //----------------------------------------------------------------//

  iftdata = _iftCreateIFTData(odisf->mimg->n); // Create an empty IFT data
  
  iter = 1;
  do
  {
    #if IFT_DEBUG //---------------------------------------------------------//
    fprintf(stderr, "iftRunODISF: Iteration %d\n", iter);
    #endif //----------------------------------------------------------------//

    // Segment image using the IFT
    _iftRunIFT(odisf, seeds, &iftdata);

    // Calculate the number of seeds to be selected (for the next iteration)
    ni_prev = ni;
    ni = iftMax(iftRound(real_n0 * exp(-iter)), odisf->nf);  
    
    if(ni < ni_prev) // If it is not the last iteration
    {
      #if IFT_DEBUG //-------------------------------------------------------//
      fprintf(stderr, "iftRunODISF: Ni = %d\n", ni);
      #endif //--------------------------------------------------------------//
      // Select the Ni relevant seeds (and remove the others)
      if(old_seeds != NULL) iftDestroyIntArray(&old_seeds);
      old_seeds = seeds;

      seeds = _iftRemSeeds(odisf, ni, iftdata, old_seeds);
      iter++;
    }
  }
  while(ni < ni_prev);

  // Creates a label image based on the seed id's in the IFT's id map
  label_img = _iftCreateLabelImageFromIdMap(odisf, seeds, iftdata);
  
  _iftDestroyIFTData(&iftdata);
  iftDestroyIntArray(&seeds);

  return label_img;
}