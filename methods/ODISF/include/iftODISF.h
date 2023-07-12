/*****************************************************************************\
* iftODISF.h
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-06-15
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#ifndef IFT_ODISF_H
#define IFT_ODISF_H

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
typedef enum ift_odisf_sampl_opt 
{
  IFT_ODISF_SAMPL_GRID, // Grid sampling (a.k.a. equally distanced seeds)
  IFT_ODISF_SAMPL_RND, // Random sampling
} iftODISFSampl;

/*
  Contains the parameters and auxiliary data structures for facilitating the
  construction of an ODISF variant algorithm.
*/
typedef struct ift_odisf_alg iftODISF;

/*****************************************************************************\
*
*                               PUBLIC FUNCTIONS
*
\*****************************************************************************/
//===========================================================================//
// CONSTRUCTOR & DESTRUCTOR
//===========================================================================//
/*
  Description

  The default values are:
    N0 = 8000, Nf = 200
    Adjacency relation = 8-adjacency
    Sampling option = Random
    No mask and no object saliency map.
*/
iftODISF *iftCreateODISF
(const iftImage *img, const iftImage *mask, const iftImage *objsm);

/*
  Description
*/
void iftDestroyODISF
(iftODISF **odisf);

//===========================================================================//
// GETTERS
//===========================================================================//
/*
  Description
*/
int iftODISFGetN0
(const iftODISF *odisf);

/*
  Description
*/
int iftODISFGetNf
(const iftODISF *odisf);

//===========================================================================//
// SETTERS
//===========================================================================//
/*
  Description
*/
void iftODISFSetN0
(iftODISF **odisf, const int n0);

/*
  Description
*/
void iftODISFSetNf
(iftODISF **odisf, const int nf);

/*
  Description
*/
void iftODISFUseDiagAdj
(iftODISF **odisf, const bool use);

//===========================================================================//
// SEED SAMPLING
//===========================================================================//
/*
  Description
*/
iftODISFSampl iftODISFGetSamplOpt
(const iftODISF *odisf);

/*
  Description
*/
void iftODISFSetSamplOpt
(iftODISF **odisf, const iftODISFSampl sampl_opt);

//===========================================================================//
// VERIFIERS
//===========================================================================//
/*
  Description
*/
bool iftODISFUsingDiagAdj
(const iftODISF *odisf);

//===========================================================================//
// RUNNER
//===========================================================================//
/*
  Description
*/
iftImage *iftRunODISF
(const iftODISF *odisf);

#ifdef __cplusplus
}
#endif

#endif