/**
 * @file
 *
 * @brief Extension of the ISF interface to perform superpixel
 *   segmentation over superpixel graphs.
 *
 * @details This file is separated into three parts identified
 *   by their main struct: iftISFOptions, iftForestAnnotationInfo
 *   and iftSuperpixelIGraph.
 *
 * @details The ISFOptions part was originally intended to simplify
 *   the access to different versions of the original ISF algorithm
 *   based on iftIGraph. It currently serves a similar purpose in the
 *   RISF framework by defining how a single segmentation step
 *   (i.e. ISF generalized to superpixel graphs) is to be performed
 *   over the iftSuperpixelIGraph.
 *
 * @details The iftForestAnnotation part is intended to provide
 *   standalone data and operations for IFT (i.e. the framework)
 *   graphs in general. Still needs some improvements to be worth
 *   using in other contexts.
 *
 * @details The iftSuperpixelIGraph defines the image graph used
 *   in single level of the RISF framework, and the relevant operations
 *   to perform  superpixel segmentation with ISF (seen as a single
 *   segmentation step in framework). It is designed to work primarily
 *   over RAGs by accessing a superpixel label map (refImg) to
 *   map the image domain to graph nodes.
 *
 * @details Note that there is no explicit hierarchy struct, each
 *   iftSuperpixelIGraph represents a single level, and adjacent levels
 *   are implicitly linked by their refImg.
 */

#ifndef _IFT_RISF_H_
#define _IFT_RISF_H_

#include "ift.h"

#ifdef __cplusplus
extern "C" {
#endif

// iftISFOptions begin

/**
 * @brief Sampling method for ISF.
 */
typedef enum ift_isf_sampling {
  IFT_ISF_UNDEFINED_SAMPLING = 0,
  IFT_ISF_GRID_SAMPLING = 1,
  IFT_ISF_MIXED_SAMPLING = 2,
  IFT_ISF_RANDOM_SAMPLING = 3,
  IFT_ISF_GEODESIC_SAMPLING = 4,
  IFT_ISF_LATTICE_SAMPLING = 5,
  IFT_ISF_HEXGRID_SAMPLING = 6
} iftISFSampling;

/**
 * @brief Non-trivial path cost calculation method (e.g. max).
 *
 * @details Prefix ADD means it is an additive IFT. Postfix ROOT, MEAN
 *   and ARC means distance computation in parametric space is against
 *   the superpixel root, superpixel mean and node offering the cost,
 *   respectively.
 */
typedef enum ift_isf_pathcost {
  /** Invalid fallback for NULL values. */
  IFT_ISF_UNDEFINED_PATHCOST = 0,
  IFT_ISF_ADD_ROOT_PATHCOST = 1,
  IFT_ISF_ADD_MEAN_PATHCOST = 2,
  IFT_ISF_ADD_ARC_PATHCOST = 3,
  IFT_ISF_MAX_MEAN_PATHCOST = 4
} iftISFPathCost;

/**
 * @brief Seed recomputation method.
 */
typedef enum ift_isf_seedrecomp {
  /** Invalid fallback for NULL values. */
  IFT_ISF_UNDEFINED_SEEDRECOMP = 0,
  /** Medoid according to some pre-specified feature distance. */
  IFT_ISF_MEDOID_SEEDRECOMP = 1,
  /** Geometric medoid. */
  IFT_ISF_CENTROID_SEEDRECOMP = 2
} iftISFSeedRecomp;


typedef struct ISFOptions {
  int nSuperpixels;
  int nIters;
  int nSmoothIters;
  /** Superpixels with size below threshold are removed in the end. */
  float smallSpThreshold; /*< Value relative to mean sp size. */
  iftISFSampling sampling;
  iftISFPathCost pathCost; /*< See iftISFPathCost. */
  iftISFSeedRecomp seedRecomp; /*< See iftISFSeedRecomp. */
} iftISFOptions;

/**
 * @brief Initialize ISFOptions with given parameters.
 *
 * @details The actual ISF methods are initialized with reasonable
 *   default values but should be changed with corresponding set
 *   functions.
 * @details The returned \c iftISFOptions is owned by the caller.
 * @see iftDestroyISFOptions()
 *
 * @param[in] nSuperpixels Target number of superpixels.
 * @param[in] nIters Number of IFT iterations.
 * @param[in] nSmoothIters Number of smoothing iterations.
 * @param[in] sampling Initial sampling method.
 * @param[in] pathCost Non-trivial path cost function.
 * @param[in] seedRecomp Seed recomputation method.
 */
iftISFOptions * iftInitISFOptions(int nSuperpixels, int nIters, int nSmoothIters, float smallSpThreshold, iftISFSampling sampling, iftISFPathCost pathCost, iftISFSeedRecomp seedRecomp);

/**
 * @brief Destroys ISF options.
 *
 * @param[in,out] options Reference to the options being destroyed.
 */
void iftDestroyISFOptions(iftISFOptions **options);

// iftISFOptions end, iftForestAnnotationInfo begin

/**
 * @brief Generic forest annotation info for graphs.
 */
typedef struct ift_forest_ann_info {
  int n;
  int nLabels;
  int *label;
  int *marker;
  int *root;
  int *pred;
  double *pathValue;
} iftForestAnnotationInfo;

/*
 * @brief Create new annotation info for \c n elements.
 *
 * @details The returned \c iftForestAnnotationInfo is owned by the caller.
 * @see iftDestroyForestAnnotationInfo()
 *
 * @param[in] n Number of elements to be anoted (e.g. nodes).
 * @return The new annotation info.
 *
 * @pre \c n is positive.
 */
iftForestAnnotationInfo * iftCreateForestAnnotationInfo(int n);

/**
 * @brief Destroys annotation info.
 *
 * @param[in,out] ann Reference to the annotation to be destroyed.
 */
void iftDestroyForestAnnotationInfo(iftForestAnnotationInfo **ann);

/**
 * @brief Copies annotation info.
 *
 * @details The returned \c iftForestAnnotationInfo is owned by the caller.
 * @see iftDestroyForestAnnotationInfo()
 *
 * @param[in] ann Reference to the annotation being copied.
 * @return Newly allocated copy of \c ann.
 */
iftForestAnnotationInfo * iftCopyForestAnnotationInfo(const iftForestAnnotationInfo *ann);

/**
 * @brief Set all nodes to non-seed initialization values.
 *
 * @param[in,out] ann Annotation to be updated.
 * @param[in] trivialPathCost Initial path cost (e.g. IFT_INFINITY_DLB).
 * @param[in] resetLabelData Flag for including label data.
 *
 * @pre \c ann is not NULL.
 */
void iftInitForestAnnotationInfo(iftForestAnnotationInfo *ann, const double trivialPathCost, const bool resetLabelData);

/**
 * @brief Initialize a single seed.
 *
 * @param[in,out] ann Annotation to be updated.
 * @param[in] seed Single seed to be initialized.
 * @param[in] trivialPathCost Initial seed path cost (e.g. 0).
 * @param[in] createNewLabels Flag to create a new label for the seed.
 *
 * @pre \c ann is not NULL and \c seed is within [0,ann->n].
 */
void iftInitForestAnnotationSingleSeed(iftForestAnnotationInfo *ann, int seed, const double seedPathCost, const bool createNewLabels);

/**
 * @brief Initializes all input seeds.
 *
 * @param[in,out] ann Annotation to be updated.
 * @param[in] seeds A set with the nodes marked as seed.
 * @param[in] trivialPathCost Initial seed path cost (e.g. 0).
 * @param[in] createNewLabels Flag to create new labels for each seed.
 *
 * @pre \c ann is not NULL and \c seeds not empty.
 */
void iftInitForestAnnotationSeeds(iftForestAnnotationInfo *ann, iftSet *seeds, const double seedPathCost, const bool createNewLabels);

// iftForestAnnotationInfo end, iftSuperpixelIGraph begin

/**
 * @brief Variation of IGraph to support superpixels.
 */
typedef struct ift_superpixel_igraph {
  /** Number of nodes. */
  int nNodes;
  /** Node adjacencies for explicit graphs. */
  iftSet **nodeAdj;
  /** Parametric space features (e.g. LAB color). */
  iftMatrix *feats; /*< Node based */
  iftMatrix *seedFeats; /*< Superpixel based */
  /** Function used to calculate distances in parametric space. */
  iftArcWeightFun distFun; /*< float(float, float, float, int) */
  /** Alpha for selecting of weighting feats in distFun. */
  float *distAlpha;
  /** Geometric space features (e.g. centroid). */
  iftMatrix *pos; /*< Node based */
  iftMatrix *seedPos; /*< Superpixel based */
  /** Path cost parameters. */
  float alpha;
  float beta;
  /** Look up table from image to nodes & spatial domain info. */
  iftImage * refImg; /*< Generally a superpixel label map. */
  long long int * spSize; /*< w.r.t. refImg/nNodes. */
  /** General Forest annotation info indexed by nodes. */
  iftForestAnnotationInfo * ann;
  /** Adjacency relation for implicit graphs. */
  iftAdjRel *A;
  /** Flag for COMPLETE, IMPLICIT or EXPLICIT graphs. */
  char type;
} iftSuperpixelIGraph;

/**
 * @brief Creates a new superpixel graph.
 *
 * @details The created graph still lacks edges and features, which
 *  are added through other functions in this header.
 * @details Number of nodes is implicitly defined as the maximum
 *  value in refImg. Additionally, if refImg has missing labels in
 *  between, independent nodes (i.e., no neighbors) might be created
 *  later on. 
 * @details The returned iftSuperpixelIGraph is owned by the caller.
 * @see iftDestroySuperpixelIGraph()
 *
 * @param[in] refImg Label image linking pixels to nodes/superpixels.
 * @return The created superpixel graph.
 *
 * @pre \c refImg is not NULL and contain positive values.
 */
iftSuperpixelIGraph *iftInitSuperpixelIGraph(const iftImage *refImg);

/**
 * @brief Destroys a superpixel graph.
 *
 * @param[in,out] igraph Reference to the graph being destroyed.
 */
void iftDestroySuperpixelIGraph(iftSuperpixelIGraph **igraph);

/**
 * @brief Copies a superpixel graph.
 *
 * @details The returned iftSuperpixelIGraph is owned by the caller.
 * @see iftDestroySuperpixelIGraph()
 *
 * @param[in] igraph Superpixel graph being copied.
 * @return Newly allocated copy of \c igraph.
 */
iftSuperpixelIGraph *iftCopySuperpixelIGraph(const iftSuperpixelIGraph *igraph);

/**
 * @brief Sets the igraph adjacency to be implicit.
 *
 * @details This actually only makes sense for pixel graphs (i.e.
 *  number of nodes equals number of pixels) and the adjacency relation
 *  will use the \c refImg for the spatial information.
 *
 * @param[in,out] igraph The superpixel graph being modified.
 * @param [in] A The implicit adjacency to be used.
 *
 * @pre \c igraph and \c A are not NULL.
 * @pre \c igraph is a pixel graph.
 */
void iftSetSuperpixelIGraphImplicitAdjacency(iftSuperpixelIGraph *igraph, iftAdjRel *A);

/**
 * @brief Sets the igraph adjacency to be neighboring superpixels.
 *
 * @param[in,out] igraph The superpixel graph being modified.
 * @param [in] A Adjacency defining neighboring superpixels.
 *
 * @pre \c igraph and \c A are not NULL.
 */
void iftSetSuperpixelIGraphExplicitRAG(iftSuperpixelIGraph *igraph, const iftAdjRel *A);

/*
 * @brief Assigns features for graph nodes.
 */
void iftSetSuperpixelIGraphFeatures(iftSuperpixelIGraph *igraph, const iftMatrix *superpixelFeatures, iftArcWeightFun distFun, float alpha, float beta);

/**
 * @brief Calculate the RISF with specified graph and options.
 *
 * @param[in,out] igraph The superpixel graph.
 * @param[in] opt The ISF framerwork options.
 * @param[in] img Base image used to build \c igraph.
 * @return The resulting superpixel label image.
 *
 * @pre All parameters are not NULL.
 * @pre igraph has been assigned a feature matrix.
 * @pre \c img and \c igraph->refImg share the same domain. 
 */
iftImage * iftSuperpixelSegmentationByRISF(iftSuperpixelIGraph *igraph, const iftISFOptions *opt, const iftImage *img);

// iftSuperpixelIGraph end
// Additional functions to compute superpixel features which
//   should be moved elsewhere later. See redundancy with
//   iftSupervoxelTo*Dataset functions.

/**
 * @brief Use geometric center as superpixel feature.
 *
 * @todo Move to a proper file dedicated to superpixel features.
 */
iftMatrix * iftComputeSuperpixelFeaturesByGeometricCenter(const iftImage *superpixelLabelMap);

/**
 * @brief Use mean of some color space as superpixel feature.
 *
 * @todo Move to a proper file dedicated to superpixel features.
 */
iftMatrix * iftComputeSuperpixelFeaturesByColorSpaceMean(const iftImage *superpixelLabelMap, const iftImage *img, iftColorSpace colorSpace);  

/**
 * @brief Use mean of mimage values as superpixel feature.
 *
 * @details For each superpixel, takes the mean (average) value of its
 *   pixels' feature vectors (i.e. \c mimg values) as superpixel
 *   feature vector.
 * @details For now we assume that the number of pixels in each
 *   superpixel is small enough so that naive float accumulation
 *   does not lead to significant errors.
 * @details The returned \c iftMatrix is owned by the caller.
 * @see iftDestroyMatrix()
 * @todo Move to a proper file dedicated to superpixel features.
 *
 * @param[in] mimg Base \c iftMImage.
 * @param[in] superpixelLabelMap Some superpixel segmentation of mimg.
 * @return Data matrix where each row corresponds to the corresponding
 * @pre \c mimg and \c superpixelLabelMap are not NULL and share the
 *      same spatial domain.
 * @pre \c superpixelLabelMap does not contain label gaps.
 */
iftMatrix * iftComputeSuperpixelFeaturesByMImageMean(const iftMImage *mimg, const iftImage *superpixelLabelMap);

#ifdef __cplusplus
}
#endif

#endif // _IFT_RISF_H_
