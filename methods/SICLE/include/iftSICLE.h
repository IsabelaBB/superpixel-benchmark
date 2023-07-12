/*****************************************************************************\
* iftSICLE.h
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-07-08
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#ifndef IFT_SICLE_H
#define IFT_SICLE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift.h"

/*****************************************************************************\
*
*                               PUBLIC STRUCTS
*
\*****************************************************************************/
/* Contains all the sampling options avaliable */
typedef enum ift_sicle_sampl_opt
{
  IFT_SICLE_SAMPL_GRID, // Grid sampling (a.k.a. equally distanced seeds)
  IFT_SICLE_SAMPL_RND, // Random sampling
} iftSICLESampl;

/* Contains all the arc-cost options avaliable */
typedef enum ift_sicle_arccost_opt
{
  IFT_SICLE_ARCCOST_ROOT, // Root-based arc-cost estimation
  IFT_SICLE_ARCCOST_DYN, // Dynamic arc-cost estimation
} iftSICLEArcCost;

/*
  Contains all the seed removal criterion avaliable. If a saliency map is
  provided, then they become object-based (except random removal)
*/
typedef enum ift_sicle_rem
{
  IFT_SICLE_REM_MAXCONTR, // Maximum Contrast-based removal
  IFT_SICLE_REM_MINCONTR, // Minimum Contrast-based removal
  IFT_SICLE_REM_SIZE, // Size-based removal
  IFT_SICLE_REM_RND, // Random removal
  IFT_SICLE_REM_MAXSC, // Size- and maximum contrast-based removal
  IFT_SICLE_REM_MINSC, // Size- and minimum contrast-based removal
} iftSICLERem;

/*
  Contains the parameters and auxiliary data structures for facilitating the
  construction of a SICLE variant algorithm.
*/
typedef struct ift_sicle_alg iftSICLE;

/*****************************************************************************\
*
*                               PUBLIC FUNCTIONS
*
\*****************************************************************************/
//===========================================================================//
// CONSTRUCTOR & DESTRUCTOR
//===========================================================================//
/*
  Creates a new SICLE object considering the given image in parameter. The user
  may provide a mask (whose dimensions must be the same as the image's), or
  simply setting it to NULL, indicating that are no spel restrictions.
  Similarly, it is possible to provide an object saliency map or set it to
  NULL.

  The default values are:
    N0 = 3000, Nf = 200, Maximum number of iterations = 5
    Adjacency relation = 8-adjacency
    No scales to output
    Arc-cost Function = Root
    Sampling option = Random
    Seed removal criterion = Size and Minimum Contrast
    No mask and no object saliency map.
*/
iftSICLE *iftCreateSICLE
(const iftImage *img, const iftImage *mask, const iftImage *objsm);

/*
  Deallocates the memory stored for the object and sets it to NULL.
*/
void iftDestroySICLE
(iftSICLE **sicle);

//===========================================================================//
// GETTERS
//===========================================================================//
/*
  Gets the current desired maximum number of iterations for segmentation
*/
int iftSICLEGetMaxIters
(const iftSICLE *sicle);

/*
  Gets the current desired initial quantity of seeds to be sampled
*/
int iftSICLEGetN0
(const iftSICLE *sicle);

/*
  Gets the current desired final number of superpixels in the segmentation
*/
int iftSICLEGetNf
(const iftSICLE *sicle);

/*
  Gets the number of superpixel segmentation scales to be computed
*/
int iftSICLEGetNumScales
(const iftSICLE *sicle);

/*
  Gets the superpixel segmentation scales to be computed
*/
int *iftSICLEGetScales
(const iftSICLE *sicle);

//===========================================================================//
// SETTERS
//===========================================================================//
/*
  Sets the desired maximum number of iterations for segmentation. Note that
  SICLE may require fewer iterations than the established maximum.
*/
void iftSICLESetMaxIters
(iftSICLE **sicle, const int max_iters);

/*
  Sets the desired initial quantity of seeds to be sampled. Depending of the
  algorithm, such number may not be guaranteed. The new quantity must be greater
  than, or equal to, the final number of superpixels.
*/
void iftSICLESetN0
(iftSICLE **sicle, const int n0);

/*
  Sets the desired final number of superpixels in the segmentation. The new
  quantity must be fewer than, or equal to, the initial number of seeds.
*/
void iftSICLESetNf
(iftSICLE **sicle, const int nf);

/*
  Sets if the number of superpixel segmentation scales and their superpixel
  quantity over the default single-result exponential computation.
*/
void iftSICLESetScales
(iftSICLE **sicle, const int num_scales, const int* scales);

/*
  Sets if the diagonal adjacency relation (i.e., 8-adjacency) is to
  be used during the segmentation. If not, then the 4-adjacency
  is used.
*/
void iftSICLEUseDiagAdj
(iftSICLE **sicle, const bool use);

/*
TODO
*/
void iftSICLEEnableBoost
(iftSICLE **sicle, const bool enable);

//===========================================================================//
// SEED SAMPLING
//===========================================================================//
/*
  Gets the current seed sampling algorithm
*/
iftSICLESampl iftSICLEGetSamplOpt
(const iftSICLE *sicle);

/*
  Sets to the desired seed sampling algorithm
*/
void iftSICLESetSamplOpt
(iftSICLE **sicle, const iftSICLESampl sampl_opt);

//===========================================================================//
// ARC COST FUNCTION
//===========================================================================//
/*
  Gets the current arc-cost function
*/
iftSICLEArcCost iftSICLEGetArcCostOpt
(const iftSICLE *sicle);

/*
  Sets to the desired arc-cost function
*/
void iftSICLESetArcCostOpt
(iftSICLE **sicle, const iftSICLEArcCost arc_opt);


//===========================================================================//
// SEED REMOVAL
//===========================================================================//
/*
  Gets the current seed removal criterion
*/
iftSICLERem iftSICLEGetRemOpt
(const iftSICLE *sicle);

/*
  Sets to the desired seed removal criterion
*/
void iftSICLESetRemOpt
(iftSICLE **sicle, const iftSICLERem rem_opt);

//===========================================================================//
// VERIFIERS
//===========================================================================//
/*
  Verifies whether the algorithm is considering the diagonal adjacents (i.e.,
  26- or 8-adjacency).
*/
bool iftSICLEUsingDiagAdj
(const iftSICLE *sicle);

/* TODO */
bool iftSICLEBoostEnabled
(const iftSICLE *sicle);

//===========================================================================//
// RUNNER
//===========================================================================//
/*
  Runs the SICLE algorithm considering the configuration and parameters defined
  in the given object.
*/
iftImage **iftRunSICLE
(const iftSICLE *sicle);

#ifdef __cplusplus
}
#endif

#endif
