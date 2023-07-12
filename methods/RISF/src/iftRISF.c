#include "iftRISF.h"

// ---------- Private Methods Declaration ----------

// iftISFOptions begin

/**
 * @brief Redirects to appropriate sampling function based on options.
 */
iftImage * iftISFPixelSampling(iftISFSampling sampling, iftMImage *mimg, iftImage *mask, int nSuperpixels);

/**
 * @brief Compute IFT with cost function based on options.
 */
void iftISFRunIFT(const iftISFOptions *options, iftIGraph *igraph, const iftImage *seeds);

// iftISFOptions end, iftForestAnnotationInfo begin

/**
 * @brief Removes a set of marked trees from graph.
 *
 * @details Frontier nodes may include duplicates and removed nodes.
 * @details Dependency to igraph is only for ease of implementation,
 *   can be removed to make the function more general.
 *
 * @param[in, out] igraph Graph whose annotation is updated. 
 * @param[in] markedTrees The set of nodes from trees to be removed.
 * @param[in] trivialPathCost Initial path cost (e.g. IFT_INFINITY_DBL).
 * @param[in, out] removedNodes Concatenates set with removed nodes.
 *                              Ignored if null.
 * @param[in, out] frontierNodes Concatenates set with frontier nodes.
 *                               Ignored if null.
 *
 * @pre \c igraph not NULL and set to IMPLICIT adjacency.
 * @pre \c markedTrees is not empty.
 */
void iftForestTreeRemovalImplicitAdj(iftSuperpixelIGraph *igraph, iftSet *markedTrees, const double trivialPathCost, iftSet **removedNodes, iftSet **frontierNodes);

/**
 * @brief Remove graph sub-tree starting on marked node.
 *
 * @details Frontier nodes may include duplicates and removed nodes.
 * @details Dependency to igraph is only for ease of implementation,
 *   can be removed to make the function more general.
 *
 * @param[in, out] igraph Graph whose annotation is updated. 
 * @param[in] markedNode Root from sub-tree being removed.
 * @param[in] trivialPathCost Cost for non-seed trivial paths.
 * @param[in, out] removedNodes Concatenates set with removed nodes.
 *                              Ignored if null.
 * @param[in, out] frontierNodes Concatenates set with frontier nodes.
 *                               Ignored if null.
 *
 * @pre \c igraph not NULL and set to IMPLICIT adjacency.
 * @pre \c markedNode is valid.
 */
void iftRemoveSubTreeImplicitAdj(iftSuperpixelIGraph *igraph, const int markedNode, const double trivialPathCost, iftSet **removedNodes, iftSet **frontierNodes);

/**
 * @brief Removes a set of marked trees from graph.
 *
 * @details Frontier nodes may include duplicates and removed nodes.
 *
 * @param[in,out] ann Annotation to be updated.
 * @param[in] nodeAdj Array of node indexed adjacencies.
 * @param[in] markedTrees The set of nodes from trees to be removed.
 * @param[in] trivialPathCost Initial path cost (e.g. IFT_INFINITY_DBL).
 * @param[in, out] removedNodes Concatenates set with removed nodes.
 *                              Ignored if null.
 * @param[in, out] frontierNodes Concatenates set with frontier nodes.
 *                               Ignored if null.
 *
 * @pre \c ann and \c nodeAdj are not NULL.
 * @pre \c markedTrees is not empty.
 */
void iftForestTreeRemovalExplicitAdj(iftForestAnnotationInfo *ann, iftSet *markedTrees, iftSet **nodeAdj, const double trivialPathCost, iftSet **removedNodes, iftSet **frontierNodes);

/**
 * @brief Remove graph sub-tree starting on marked node.
 *
 * @details Frontier nodes may include duplicates and removed nodes.
 *
 * @param[in, out] ann Annotation to be updated
 * @param[in] markedNode Root from sub-tree being removed.
 * @param[in] nodeAdj Explicit graph adjacency.
 * @param[in] trivialPathCost Cost for non-seed trivial paths.
 * @param[in, out] removedNodes Concatenates set with removed nodes.
 *                              Ignored if null.
 * @param[in, out] frontierNodes Concatenates set with frontier nodes.
 *                               Ignored if null.
 *
 * @pre \c ann and \c nodeAdj are not NULL.
 * @pre \c markedNode is valid.
 */
void iftRemoveSubTreeExplicitAdj(iftForestAnnotationInfo *ann, const int markedNode, iftSet **nodeAdj, const double trivialPathCost, iftSet **removedNodes, iftSet **frontierNodes);

// iftForestAnnotationInfo end, iftSuperpixelIGraph begin

/*
 * @brief Initializes auxiliary features for mean path cost.
 *
 * @details As there is no superpixel grouping at the start,
 *   initialize using the seed features.
 *
 * @param[in] igraph Working superpixel graph.
 * @param[in] seeds Initial sampling seeds.
 */
void iftInitSuperpixelIGraphMeanAuxFeat(iftSuperpixelIGraph *igraph, iftSet *seeds);

/*
 * @brief Computes auxiliary features for mean path cost.
 *
 * @details Currently uses the mean from previous iteration. Can be
 *   extended to update while running IFT.
 *
 * @param[in] igraph Working superpixel graph.
 */
void iftUpdateSuperpixelIGraphMeanAuxFeat(const iftSuperpixelIGraph *igraph);

/**
 * @brief Returns the medoid node from each superpixel of the igraph.
 *
 * @details The superpixels we mention here are the new ones being
 *  calculated which are contained in igraph's forest annotation.
 *
 * @param[in] igraph Graph from which medoids are being computed.
 * @param[in] recompFeats Selection of features used for seed recomp.
 *                        Can be NULL to use all valid features.
 * @return Set with the medoid of each superpixel. 
 */
iftSet * iftRISFMedoidSeedRecomp(const iftSuperpixelIGraph *igraph, bool isParametric);

/*
 * @brief Filter candidate seeds below adaptative threshold.
 *
 * @details Adaptative threshold is given by the square root of the
 *   mean distance from all nodes to their respective roots. If
 *   candidate seed is below threshold, it is reverted to the previous
 *   iteration seed.
 *
 * @param[in] igraph Working superpixel graph.
 * @param[in,out] seeds Candidate seeds from a seedRecomp function.
 */
void iftSeedUpdateThreshold(const iftSuperpixelIGraph *igraph, iftSet *seeds);

/**
 * @brief Build a mask with only the medoid of each superpixel.
 */
iftImage * iftSelectSuperpixelCenters(const iftMImage *mimg, const iftImage *label);

/**
 * @brief Convert an IGraph generated from superpixels into label image.
 */
iftImage * iftHierarchyIGraphLabel(iftIGraph *igraph, const iftImage *predLabel);

/**
 * @brief Computes the cost of extending the path pi_s to node \c t.
 *
 * @param[in] igraph Graph in which the transaction occurs.
 * @param[in] s Node offering the path.
 * @param[in] t Node receiving the offer.
 * @param[in] optISF ISF options containing path cost information.
 * @return TRUE if the offer is accepted and FALSE otherwise.
 *
 * @pre \c igraph and \c optISF are not NULL.
 * @pre \c s and \c are non-negative (i.e. valid nodes).
 * @post Returns a non-negative distance.
 */
double iftRISFGetPathCost(const iftSuperpixelIGraph *igraph, const int s, const int t, const iftISFOptions *optISF);

/**
 * @brief Compute node \c s offering new path to node \c t in ISF.
 *
 * @param[in,out] igraph Graph in which the transaction occurs.
 * @param[in] s Node offering the path.
 * @param[in] t Node receiving the offer.
 * @param[in] optISF ISF options containing path cost information.
 * @return 1 if offer is better, 0 if equal and -1 if worse.
 *
 * @pre \c igraph and \c optISF are not NULL.
 * @pre \c s and \c are non-negative (i.e. valid nodes).
 */
int iftRISFOfferNewPath(iftSuperpixelIGraph *igraph, const int s, const int t, const iftISFOptions *optISF);

void iftRISFInitializeDIFT(iftSuperpixelIGraph *igraph, const iftISFOptions *opt, iftDHeap *Q, double trivialPathCost);

void iftRISFRemoveInconsistentSubTree(iftSuperpixelIGraph *igraph, iftDHeap *Q, double trivialPathCost, int t);

void iftRISFNeighborProcessingDIFT(iftSuperpixelIGraph *igraph, iftDHeap *Q, double trivialPathCost, int s, int t, const iftISFOptions *opt);

/**
 * @brief Consume removed and frontier nodes to update queue for DIFT.
 */
void iftRISFConsumeRemovedAndFrontierNodes(iftSuperpixelIGraph *igraph, iftDHeap *Q, double trivialPathCost, iftSet **removedNodes, iftSet **frontierNodes);

/**
 * @brief Helper function for DIFT tree removal.
 */
void iftRISFTreeRemoval(iftSuperpixelIGraph *igraph, iftDHeap *Q, iftSet *markedNodes, double trivialPathCost);

/**
 * @brief Remove small trees and their seeds.
 *
 * @param[in] threshold Percent relative to mean superpixel size.
 */
void iftRISFRemoveSmallSuperpixels(iftSuperpixelIGraph *igraph, iftDHeap *Q, float threshold, double trivialPathCost);

/**
 * @brief Fix gaps in label data generated by removed seeds.
 */
void iftRemoveLabelGaps(int *labelArray, int size);

/**
 * @brief Redirect to correct sampling function based on \c opt.
 *
 * @param[in] igraph Graph from which nodes are being sampled.
 * @param[in] opt Options defining how sampling should be done.
 * @param[in] img Original image from which \c igraph was created.
 * @return Set of node indexed seeds.
 *
 * @pre All parameters are not NULL.
 * @pre \c img has the same domain of \c igraph->refImg.
 */
iftSet * iftRISFSampling(const iftSuperpixelIGraph *igraph, const iftISFOptions *opt, const iftImage *img);

/**
 * @brief Samples nodes over a square lattice of the spatial domain.
 *
 * @details We only check in which superpixel each grid pixel fall to
 *  get our samples. If the sampling is too dense, two grid points
 *  might fall into the same superpixel and less samples are returned
 *  as a result.
 *
 * @param[in] igraph Graph from which nodes are being sampled.
 * @param[in] img Original image from which \c igraph was created.
 * @param[in] nSamples Number of nodes being sampled.
 * @return Set of node indexed seeds.
 *
 * @pre \c igraph and \c img parameters are not NULL.
 * @pre \c img has the same domain of \c igraph->refImg.
 * @pre \c nSamples is within (0, igraph->nNodes).
 */
iftSet * iftRISFGridSampling(const iftSuperpixelIGraph *igraph, const iftImage *img, int nSamples);

/**
 * @brief Samples nodes biased by high entropy areas.
 *
 * @param[in] igraph Graph from which nodes are being sampled.
 * @param[in] img Original image from which \c igraph was created.
 * @param[in] nSamples Number of nodes being sampled.
 * @return Set of node indexed seeds.
 *
 * @pre \c igraph and \c img parameters are not NULL.
 * @pre \c img has the same domain of \c igraph->refImg.
 * @pre \c nSamples is within (0, igraph->nNodes).
 */
iftSet * iftRISFMixedSampling(const iftSuperpixelIGraph *igraph, const iftImage *img, int nSamples);

/**
 * @brief Samples random nodes from the igraph.
 *
 * @param[in] igraph Graph from which nodes are being sampled.
 * @param[in] nSamples Number of nodes being sampled.
 * @return Set of node indexed seeds.
 *
 * @pre \c igraph is not NULL.
 * @pre \c nSamples is within (0, igraph->nNodes).
 */
iftSet * iftRISFRandomSampling(const iftSuperpixelIGraph *igraph, int nSamples);

/**
 * @brief Samples all nodes from the igraph.
 *
 * @param[in] igraph Graph from which nodes are being sampled.
 *
 * @pre \c igraph is not NULL.
 */
iftSet * iftRISFAllNodesSampling(const iftSuperpixelIGraph *igraph);

/**
 * @brief Recursively obtains new samples from max geodesic distance.
 *
 * @details This method is defined by induction on the number of
 *   seeds. For the first sample (base case) we pick an arbitrary
 *   node as seed. We obtain a new seed (induction step) by
 *   computing the node with maximum geodesic distance from the
 *   current seed set. This is done with an additive IFT over the
 *   superpixel mean brightness.
 *
 * @param[in] igraph Graph from which nodes are being sampled.
 * @param[in] img Original image for feature extraction.
 * @param[in] nSamples Number of nodes being sampled.
 * @return Set of node indexed seeds.
 *
 * @pre \c igraph is not NULL.
 * @pre \c nSamples is within (0, igraph->nNodes).
 */
iftSet * iftRISFGeodesicSampling(const iftSuperpixelIGraph *igraph, const iftImage *img, int nSamples);

/*
 * @brief Performs a regular lattice sampling on arbitrary masks.
 *
 * @details See iftGridSamplingOnMask. We first estimate a radius
 *   within a 5% error of the requested number of samples and
 *   then obtain the final sampling.
 *
 * @param[in] igraph Graph from which nodes are being sampled.
 * @param[in] nSamples Number of nodes being sampled.
 *
 * @pre \c igraph is not NULL.
 * @pre \c nSamples is within (0, igraph->nNodes).
 */
iftSet * iftRISFGridOnMaskSampling(const iftSuperpixelIGraph *igraph, int nSamples);

/**
 * @brief Samples nodes over a hexagonal lattice of the spatial domain.
 *
 * @details Analogous to square lattice grid sampling.
 *
 * @param[in] igraph Graph from which nodes are being sampled.
 * @param[in] nSamples Number of nodes being sampled.
 * @return Set of node indexed seeds.
 *
 * @pre \c igraph is not NULL.
 * @pre \c nSamples is within (0, igraph->nNodes).
 */
iftSet * iftRISFHexGridSampling(const iftSuperpixelIGraph *igraph, int nSamples);

/*
 * @brief Removes duplicate seeds from selection.
 *
 * @details Temporary method for easier implementation, should make a
 *   more general one that works for any iftSet.
 *
 * @param[in,out] seeds Set of seeds from a seed sampling method.
 * @param[in] nNodes Number of igraph nodes for range of seed values.
 */
void iftRISFRemoveDuplicateSeeds(iftSet **seeds, int nNodes);

/**
 * @brief Graph initialization for Geodesic sampling IFT.
 */
iftSuperpixelIGraph * iftBuildGeodesicSamplingIGraph(const iftSuperpixelIGraph *igraph, const iftImage *img);

/**
 * @brief Seed initialization for Geodesic sampling IFT.
 */
void iftGeodesicSamplingSeedInit(iftSuperpixelIGraph *igraph, iftDHeap *Q, iftDHeap *R, iftSet *seeds, double trivialPathCost, double seedPathCost, int iter);

/**
 * @brief Helper function for Geodesic sampling IFT.
 */
void iftGeodesicSamplingNeighborProcessing(iftSuperpixelIGraph *igraph, iftDHeap *Q, iftDHeap *R, int s, int t);

/**
 * @brief DIFT sub-tree removal procedure for geodesic sampling.
 */
void iftGeodesicSamplingRemoveInconsistentSubTree(iftSuperpixelIGraph *igraph, iftDHeap *Q, iftDHeap *R, int t);

// ---------- Private Methods Implementation ----------

// iftISFOptions begin

iftImage * iftISFPixelSampling(iftISFSampling sampling, iftMImage *mimg, iftImage *mask, int nSuperpixels)
{
  switch(sampling) {
    case IFT_ISF_GRID_SAMPLING:
      return iftGridSampling(mimg, mask, nSuperpixels);
    case IFT_ISF_MIXED_SAMPLING:
      return iftAltMixedSampling(mimg, mask, nSuperpixels);
    default:
      iftError("Sampling strategy not selected.", "iftISFPixelSampling");
  }
  return NULL;
}

void iftISFRunIFT(const iftISFOptions *options, iftIGraph *igraph, const iftImage *seeds)
{
  // Deprecated
  /*
  switch(options->costFunction) {
    case IFT_ISF_ROOT_COSTFUNCTION:
      {
        float alpha = options->costParams.additive.alpha;
        float beta = options->costParams.additive.beta;
        iftIGraphISF_Root(igraph, seeds, alpha, beta, options->nIters);
      }
      break;
    case IFT_ISF_MEAN_COSTFUNCTION:
      {
        float alpha = options->costParams.additive.alpha;
        float beta = options->costParams.additive.beta;
        iftIGraphISF_Mean(igraph, seeds, alpha, beta, options->nIters);
      }
      break;
    case IFT_ISF_REGMIN_COSTFUNCTION:
      iftError("REGMIN not added yet", "iftSuperpixelSegmentationByISF");
    default:
      iftError("Cost function not selected.", "iftISFRunIFT");
  }
  */
}

// iftISFOptions end, iftForestAnnotationInfo begin

void iftForestTreeRemovalImplicitAdj(iftSuperpixelIGraph *igraph, iftSet *markedTrees, const double trivialPathCost, iftSet **removedNodes, iftSet **frontierNodes)
{
  assert(igraph != NULL && igraph->type == IMPLICIT);

  int ct = 0;
  for (iftSet *s = markedTrees; s != NULL; s = s->next) {
    ct++;
    int p = s->elem;

    // Skip already removed trees
    if (igraph->ann->pathValue[p] == trivialPathCost)
      continue;

    int r = igraph->ann->root[p];

    // Sub-tree starting from root is entire tree
    iftRemoveSubTreeImplicitAdj(igraph, r,
        trivialPathCost, removedNodes, frontierNodes);
  }
}

void iftRemoveSubTreeImplicitAdj(iftSuperpixelIGraph *igraph, const int markedNode, const double trivialPathCost, iftSet **removedNodes, iftSet **frontierNodes)
{
  assert(igraph != NULL && igraph->type == IMPLICIT);
  assert(markedNode >= 0 && markedNode < igraph->ann->n);

  iftSet * pendingNodes = NULL;
  iftInsertSet(&pendingNodes, markedNode);

  while (pendingNodes != NULL) {
    int p = iftRemoveSet(&pendingNodes);
    if (removedNodes)
      iftInsertSet(removedNodes, p);

    // Skip already removed nodes
    if (igraph->ann->pathValue[p] == trivialPathCost)
      continue;

    igraph->ann->root[p] = IFT_NIL;
    igraph->ann->pred[p] = IFT_NIL;
    igraph->ann->pathValue[p] = trivialPathCost;
    
    // p is not a proper pixel index when using mask
    iftVoxel u = { 
      iftMatrixRowPointer(igraph->pos, p)[0],
      iftMatrixRowPointer(igraph->pos, p)[1],
      (igraph->pos->ncols > 2) ?
        iftMatrixRowPointer(igraph->pos, p)[2] : 0 };
    for (int i = 1; i < igraph->A->n; ++i) {
      iftVoxel v = iftGetAdjacentVoxel(igraph->A, u, i);
      if (!iftValidVoxel(igraph->refImg, v))
        continue;

      int q = iftGetVoxelIndex(igraph->refImg, v);
      q = igraph->refImg->val[q] - 1;
      if (q < 0)
        continue;

      if (igraph->ann->pred[q] == p)
        iftInsertSet(&pendingNodes, q);
      else if (frontierNodes && igraph->ann->pathValue[q] != trivialPathCost)
        iftInsertSet(frontierNodes, q);
    }
  }
}

void iftForestTreeRemovalExplicitAdj(iftForestAnnotationInfo *ann, iftSet *markedTrees, iftSet **nodeAdj, const double trivialPathCost, iftSet **removedNodes, iftSet **frontierNodes)
{
  assert(ann != NULL && nodeAdj != NULL);

  for (iftSet *s = markedTrees; s != NULL; s = s->next) {
    int p = s->elem;

    // Skip already removed trees
    if (ann->pathValue[p] == trivialPathCost)
      continue;

    int r = ann->root[p];

    // Sub-tree starting from root is entire tree
    iftRemoveSubTreeExplicitAdj(ann, r, nodeAdj,
        trivialPathCost, removedNodes, frontierNodes);
  }
}

void iftRemoveSubTreeExplicitAdj(iftForestAnnotationInfo *ann, const int markedNode, iftSet **nodeAdj, const double trivialPathCost, iftSet **removedNodes, iftSet **frontierNodes)
{
  assert(ann != NULL && nodeAdj != NULL);
  assert(markedNode >= 0 && markedNode < ann->n);

  iftSet * pendingNodes = NULL;
  iftInsertSet(&pendingNodes, markedNode);

  while (pendingNodes != NULL) {
    int p = iftRemoveSet(&pendingNodes);
    if (removedNodes)
      iftInsertSet(removedNodes, p);

    // Skip already removed nodes
    if (ann->pathValue[p] == trivialPathCost)
      continue;

    ann->root[p] = IFT_NIL;
    ann->pred[p] = IFT_NIL;
    ann->pathValue[p] = trivialPathCost;
    
    for (iftSet *adj = nodeAdj[p]; adj != NULL; adj = adj->next) {
      int q = adj->elem;
      if (ann->pred[q] == p)
        iftInsertSet(&pendingNodes, q);
      else if (frontierNodes && ann->pathValue[q] != trivialPathCost)
        iftInsertSet(frontierNodes, q);
    }
  }
}

// iftForestAnnotationInfo end, iftSuperpixelIGraph begin

void iftInitSuperpixelIGraphMeanAuxFeat(iftSuperpixelIGraph *igraph, iftSet *seeds)
{
  assert(igraph != NULL);
  assert(seeds != NULL);

  // Copy seed features directly
  iftMatrix *seedFeats = iftCreateMatrix(igraph->feats->ncols, igraph->ann->nLabels);
  iftMatrix *seedPos = iftCreateMatrix(igraph->pos->ncols, igraph->ann->nLabels);
  int row = 0;
  for (; seeds != NULL; seeds = seeds->next, ++row) {
    for (int col = 0; col < igraph->feats->ncols; ++col)
      iftMatrixElem(seedFeats, col, row) =
        iftMatrixElem(igraph->feats, col, seeds->elem);
    for (int col = 0; col < igraph->pos->ncols; ++col)
      iftMatrixElem(seedPos, col, row) =
        iftMatrixElem(igraph->pos, col, seeds->elem);
  }

  // Assing aux features to graph
  if (igraph->seedFeats != NULL)
    iftDestroyMatrix(&(igraph->seedFeats));
  if (igraph->seedPos != NULL)
    iftDestroyMatrix(&(igraph->seedPos));
  igraph->seedFeats = seedFeats;
  igraph->seedPos = seedPos;
}

void iftUpdateSuperpixelIGraphMeanAuxFeat(const iftSuperpixelIGraph *igraph)
{
  assert(igraph != NULL);
  assert(igraph->seedFeats != NULL);
  assert(igraph->seedPos != NULL);

  // Shorthands and initialization
  iftMatrix *feats = igraph->feats;
  iftMatrix *pos = igraph->pos;
  iftMatrix *meanFeats = igraph->seedFeats;
  iftMatrix *meanPos = igraph->seedPos;
  for (int i = 0; i < meanFeats->n; ++i)
    meanFeats->val[i] = 0.0;

  // TODO Eliminate redundancy with MedoidSeedRecomp & duplicate
  // Find mean for each label
  int *labelCount = iftAllocIntArray(igraph->seedFeats->nrows);
  // Sum all features for each label
  for (int node = 0; node < igraph->nNodes; ++node) {
    int label = igraph->ann->label[node] - 1;
    if (label < 0)
      continue;
    labelCount[label] += 1;
    // Naive accumulation
    for (int col = 0; col < feats->ncols; ++col)
      iftMatrixElem(meanFeats, col, label) +=
        iftMatrixElem(feats, col, node);
    for (int col = 0; col < pos->ncols; ++col)
      iftMatrixElem(meanPos, col, label) +=
        iftMatrixElem(pos, col, node);
  }
  // Divide by count to obtain mean
  for (int row = 0; row < meanFeats->nrows; ++row) {
    for (int col = 0; col < meanFeats->ncols; ++col)
      iftMatrixElem(meanFeats, col, row) /=
        ((float) labelCount[row]);
    for (int col = 0; col < meanPos->ncols; ++col)
      iftMatrixElem(meanPos, col, row) /=
        ((float) labelCount[row]);
  }
}

iftSet * iftRISFMedoidSeedRecomp(const iftSuperpixelIGraph *igraph, bool isParametric)
{
  assert(igraph != NULL);
  assert(igraph->ann->nLabels > 0);

  int nLabels = igraph->ann->nLabels; // shorthand
  double *nodeToMeanDist = iftAllocDoubleArray(igraph->nNodes);

  // Find mean for each label
  iftMatrix *Mx = isParametric ? igraph->feats : igraph->pos;
  iftMatrix *meanMx = iftCreateMatrix(Mx->ncols, nLabels);
  int *labelCount = iftAllocIntArray(nLabels);
  // Sum all features for each label
  for (int node = 0; node < igraph->nNodes; ++node) {
    int label = igraph->ann->label[node] - 1;
    if (label < 0)
      continue;
    labelCount[label] += igraph->spSize[node];
    // Original seed recomputation without size weight
    /*labelCount[label] += 1;*/
    // Naive accumulation
    for (int col = 0; col < Mx->ncols; ++col)
      meanMx->val[iftGetMatrixIndex(meanMx, col, label)] +=
        Mx->val[iftGetMatrixIndex(Mx, col, node)] * igraph->spSize[node];
    // Original seed recomputation without size weight
    /*for (int col = 0; col < Mx->ncols; ++col)
      meanMx->val[iftGetMatrixIndex(meanMx, col, label)] +=
        Mx->val[iftGetMatrixIndex(Mx, col, node)];*/
  }
  // Divide by count to obtain mean
  for (int row = 0; row < meanMx->nrows; ++row)
    for (int col = 0; col < meanMx->ncols; ++col)
      meanMx->val[iftGetMatrixIndex(meanMx, col, row)] /=
        ((float) labelCount[row]);

  // Compute each node's distance to its label mean
  for (int node = 0; node < igraph->nNodes; ++node) {
    int label = igraph->ann->label[node] - 1;
    if (label < 0)
      continue;
    nodeToMeanDist[node] = igraph->distFun(
        iftMatrixRowPointer(Mx, node),
        iftMatrixRowPointer(meanMx, label),
        igraph->distAlpha,
        Mx->ncols);
  }

  // Find node with smallest distance to each label mean (i.e. medoid)
  int *selectedNode = iftAllocIntArray(nLabels);
  double *selectedNodeDistance = iftAllocDoubleArray(nLabels);
  for (int i = 0; i < nLabels; ++i) {
    selectedNode[i] = -1;
    selectedNodeDistance[i] = IFT_INFINITY_DBL;
  }
  for (int node = 0; node < igraph->nNodes; ++node) {
    int label = igraph->ann->label[node] - 1;
    if (label < 0)
      continue;
    if (nodeToMeanDist[node] < selectedNodeDistance[label]) {
      selectedNode[label] = node;
      selectedNodeDistance[label] = nodeToMeanDist[node];
    }
  }

  // Convert to set
  iftSet *samples = NULL;
  for (int i = 0; i < nLabels; ++i)
    if (selectedNode[i] >= 0)
      iftInsertSet(&samples, selectedNode[i]);

  // Clean up
  free(nodeToMeanDist);
  free(selectedNode);
  free(selectedNodeDistance);
  iftDestroyMatrix(&meanMx);
  free(labelCount);

  return samples;
}

void iftSeedUpdateThreshold(const iftSuperpixelIGraph *igraph, iftSet *seeds)
{
  assert(igraph != NULL);
  assert(seeds != NULL);

  // Compute overall mean distance to root
  double pMeanDist = 0.0; // Parametric distance
  double spatialMeanDist = 0.0; // Spatial distance
  // Accumulation
  for (int node = 0; node < igraph->nNodes; ++node) {
    // kek
    pMeanDist += igraph->distFun(
        iftMatrixRowPointer(igraph->feats, node),
        iftMatrixRowPointer(igraph->feats, igraph->ann->root[node]),
        igraph->distAlpha,
        igraph->feats->ncols);
    spatialMeanDist += iftFeatDistance(
      iftMatrixRowPointer(igraph->pos, node),
      iftMatrixRowPointer(igraph->pos, igraph->ann->root[node]),
      igraph->pos->ncols);
  }
  // Actual mean
  pMeanDist /= (float) igraph->nNodes;
  spatialMeanDist /= (float) igraph->nNodes;

  // Obtain actual threshold as square root of mean distances
  double pThresh = sqrtf(pMeanDist);
  double spatialThresh = sqrtf(spatialMeanDist);

  // Revert candidates below threshold
  for (; seeds != NULL; seeds = seeds->next) {
    int node = seeds->elem;
    double pDist = igraph->distFun(
        iftMatrixRowPointer(igraph->feats, node),
        iftMatrixRowPointer(igraph->feats, igraph->ann->root[node]),
        igraph->distAlpha,
        igraph->feats->ncols);
    double spatialDist = iftFeatDistance(
      iftMatrixRowPointer(igraph->pos, node),
      iftMatrixRowPointer(igraph->pos, igraph->ann->root[node]),
      igraph->pos->ncols);

    if (pDist < pThresh && spatialDist < spatialThresh)
      seeds->elem = igraph->ann->root[node];
  }
}

void iftRISFSmoothing(iftSuperpixelIGraph *igraph, const iftImage *img, int nSmoothIters)
{
  if (igraph->nNodes != igraph->refImg->n || nSmoothIters <= 0)
    return;

  // For now we offload to existing igraph implementation
  // -- Create corresponding iftIGraph
  iftMImage *mimg; 
  if (iftIsColorImage(img))
    mimg = iftImageToMImage(img, LABNorm_CSPACE);
  else
    mimg = iftImageToMImage(img, GRAY_CSPACE);  
  iftAdjRel *A = iftIs3DMImage(mimg) ? iftSpheric(1.0) : iftCircular(1.0);
  iftImage *mask = iftSelectImageDomain(mimg->xsize, mimg->ysize, mimg->zsize);
  iftIGraph *pigraph = iftImplicitIGraph(mimg, mask, A);

  // Copy superpixel igraph labels
  for (int node = 0; node < igraph->nNodes; ++node)
    pigraph->label[node] = igraph->ann->label[node];

  // Use existing smoothing functions
  iftIGraphSetWeightForRegionSmoothing(pigraph, img);
  iftIGraphSmoothRegions(pigraph, nSmoothIters);

  // Retrieve smoothed label map
  for (int node = 0; node < igraph->nNodes; ++node)
    igraph->ann->label[node] = pigraph->label[node];

  iftDestroyAdjRel(&A);
  iftDestroyMImage(&mimg);
  iftDestroyImage(&mask);
  iftDestroyIGraph(&pigraph);
}

// TODO merge redundancy with iftIGraphSuperpixelCenters
iftImage * iftSelectSuperpixelCenters(const iftMImage *mimg, const iftImage *label)
{
  int nseeds = iftMaximumValue(label);
  float **feat = (float **) iftAlloc(nseeds, sizeof(*feat));
  float *nelems = iftAllocFloatArray(nseeds); 
  int *center = iftAllocIntArray(nseeds);

  for (int i = 0; i < nseeds; ++i) {
    feat[i] = iftAllocFloatArray(mimg->m);
    center[i] = -1;
  }

  // Find average superpixel feature value
  for (int p = 0; p < mimg->n; ++p) {
    int i = label->val[p] - 1;
    nelems[i]++;
    if (center[i] < 0)
      center[i] = p;
    for (int j = 0; j < mimg->m; ++j)
      feat[i][j] += mimg->val[p][j];
  }
  for (int i = 0; i < nseeds; ++i)
    for (int j = 0; j < mimg->m; ++j)
      feat[i][j] /= nelems[i];

  // Find pixel closest to the center
  for (int p = 0; p < mimg->n; ++p) {
    int i = label->val[p] - 1;
    int q = center[i];
    float dist1 = 0.0;
    float dist2 = 0.0;
    for (int j = 0; j < mimg->m; ++j) {
      float tmp1 = feat[i][j] - mimg->val[p][j];
      dist1 += tmp1 * tmp1;
      float tmp2 = feat[i][j] - mimg->val[q][j];
      dist2 += tmp2 * tmp2;
    }
    dist1 = sqrt(dist1);
    dist2 = sqrt(dist2);
    if (dist1 < dist2)
      center[i] = p;
  }

  // Build corresponding mask
  iftImage *mask = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
  for (int i = 0; i < nseeds; ++i)
    mask->val[center[i]] = 1;

  for (int i = 0; i < nseeds; ++i)
    iftFree(feat[i]);
  iftFree(feat);
  iftFree(nelems);
  iftFree(center);

  return mask;
}

iftImage * iftHierarchyIGraphLabel(iftIGraph *igraph, const iftImage *predLabel)
{
  int npred = iftMaximumValue(predLabel);
  int *conversionTb = iftAllocIntArray(npred);
  for (int s = 0; s < igraph->nnodes; ++s) {
    int p = igraph->node[s].voxel;
    int r = igraph->root[p]; 
    conversionTb[predLabel->val[p] - 1] = igraph->label[r];
  }

  iftImage *label = iftCopyImage(predLabel);
  for (int i = 0; i < label->n; ++i)
    label->val[i] = conversionTb[label->val[i]-1];

  iftFree(conversionTb);
  return label;
}

double iftRISFGetPathCost(const iftSuperpixelIGraph *igraph, const int s, const int t, const iftISFOptions *optISF)
{
  const char *funName = "iftRISFGetPathCost";
  assert(igraph != NULL && optISF != NULL);
  assert(s >= 0 && t >= 0);

  // Cost of path up to s
  double predCost = igraph->ann->pathValue[s];

  // Parametric space distance calculation
  double pCost = 0.0;
  switch(optISF->pathCost) {
    case IFT_ISF_ADD_ROOT_PATHCOST:
      pCost = igraph->distFun(
          iftMatrixRowPointer(igraph->feats, igraph->ann->root[s]),
          iftMatrixRowPointer(igraph->feats, t),
          igraph->distAlpha,
          igraph->feats->ncols);
      break;
    case IFT_ISF_ADD_MEAN_PATHCOST:
        pCost = igraph->distFun(
            iftMatrixRowPointer(igraph->seedFeats, igraph->ann->label[s]-1),
            iftMatrixRowPointer(igraph->feats, t),
            igraph->distAlpha,
            igraph->feats->ncols);
        break;
    case IFT_ISF_ADD_ARC_PATHCOST:
      pCost = igraph->distFun(
          iftMatrixRowPointer(igraph->feats, s),
          iftMatrixRowPointer(igraph->feats, t),
          igraph->distAlpha,
          igraph->feats->ncols);
      break;
    default:
      iftError("Invalid pathCostNodes.", funName);
  }

  double distCost = iftFeatDistance(
      iftMatrixRowPointer(igraph->pos, s),
      iftMatrixRowPointer(igraph->pos, t),
      igraph->pos->ncols);

  double edgeCost = pow(igraph->alpha * pCost, igraph->beta) + distCost;
  // Experimental, weight superpixel size
  //edgeCost *= igraph->spSize[s];
  //edgeCost *= iftMin(igraph->spSize[s], igraph->spSize[t]);

  // Compute final cost
  double cost = 0.0;
  switch(optISF->pathCost) {
    case IFT_ISF_ADD_ROOT_PATHCOST:
    case IFT_ISF_ADD_MEAN_PATHCOST:
    case IFT_ISF_ADD_ARC_PATHCOST:
    default:
      cost = predCost + edgeCost;
  }

  assert(cost >= 0.0);
  return cost;
}

int iftRISFOfferNewPath(iftSuperpixelIGraph *igraph, const int s, const int t, const iftISFOptions *optISF)
{
  assert(igraph != NULL && optISF != NULL);
  assert(s >= 0 && t >= 0);

  double offeredCost = iftRISFGetPathCost(igraph, s, t, optISF);

  if (offeredCost > igraph->ann->pathValue[t])
    return -1;
  else if (offeredCost == igraph->ann->pathValue[t])
    return 0;

  igraph->ann->pathValue[t] = offeredCost;
  igraph->ann->root[t] = igraph->ann->root[s];
  igraph->ann->label[t] = igraph->ann->label[s];
  igraph->ann->pred[t] = s;

  return 1;
}

void iftRISFInitializeDIFT(iftSuperpixelIGraph *igraph, const iftISFOptions *opt, iftDHeap *Q, double trivialPathCost)
{
  const char *funName = "iftRISFInitializeDIFT";

  // Seed recomputation step
  iftSet *seeds = NULL;
  switch(opt->seedRecomp) {
    case IFT_ISF_MEDOID_SEEDRECOMP:
      seeds = iftRISFMedoidSeedRecomp(igraph, true);
      break;
    case IFT_ISF_CENTROID_SEEDRECOMP:
      // Centroid is equivalent to spatial medoid
      seeds = iftRISFMedoidSeedRecomp(igraph, false);
      break;
    default:
      iftError("Invalid seed recomputation method.", funName);
  }

  // Revert seed updates under default threshold
  iftSeedUpdateThreshold(igraph, seeds);
  
  // DIFT - Tree Removal & queue update
  // Isolate updated seeds
  iftSet *updatedSeeds = NULL;
  for (iftSet *s = seeds; s != NULL; s = s->next)
    if (s->elem != igraph->ann->root[s->elem])
      iftInsertSet(&updatedSeeds, s->elem);

  // Remove trees with update seeds and re-initialize
  //   queue accordingly.
  iftRISFTreeRemoval(igraph, Q, updatedSeeds, trivialPathCost);

  // Add updated seeds
  // TODO add seedPathCost as a variable
  if (updatedSeeds != NULL)
    iftInitForestAnnotationSeeds(igraph->ann, updatedSeeds, 0.0, false);
  while (updatedSeeds != NULL)
    iftInsertDHeap(Q, iftRemoveSet(&updatedSeeds));

  iftDestroySet(&seeds);
}

void iftRISFRemoveInconsistentSubTree(iftSuperpixelIGraph *igraph, iftDHeap *Q, double trivialPathCost, int t)
{
  const char *funName = "iftRISFRemoveInconsistentSubTree";
  iftSet *removedNodes = NULL;
  iftSet *frontierNodes = NULL;
  switch (igraph->type) {
    case IMPLICIT:
      iftRemoveSubTreeImplicitAdj(igraph, t,
          trivialPathCost, &removedNodes, &frontierNodes);
      break;
    case EXPLICIT:
      iftRemoveSubTreeExplicitAdj(igraph->ann, t,
          igraph->nodeAdj, trivialPathCost,
          &removedNodes, &frontierNodes);
      break;
    case COMPLETE:
      // TODO
    default:
      iftError("Invalid graph type.", funName);
  }
  iftRISFConsumeRemovedAndFrontierNodes(igraph, Q, 
      trivialPathCost, &removedNodes, &frontierNodes);
}

void iftRISFNeighborProcessingDIFT(iftSuperpixelIGraph *igraph, iftDHeap *Q, double trivialPathCost, int s, int t, const iftISFOptions *opt)
{
  // TEMPORARY TO REPRODUCE OLD RESULT
  //if (Q->color[t] == IFT_BLACK)
  //  return;

  int offerStatus = iftRISFOfferNewPath(igraph, s, t, opt);
  if (offerStatus > 0) {
    // Better cost, standard IFT behavior
    if (Q->color[t] == IFT_GRAY)
      iftGoUpDHeap(Q, Q->pos[t]);
    else
      iftInsertDHeap(Q, t);
  } else {
    if (igraph->ann->pred[t] == s) {
      // DIFTmod with simplified tie processing
      int sLabel = igraph->ann->label[s];
      int tLabel = igraph->ann->label[t];
      if (offerStatus < 0 || sLabel != tLabel)
        iftRISFRemoveInconsistentSubTree(igraph, Q, trivialPathCost, t);
    }
  }
}

void iftRISFConsumeRemovedAndFrontierNodes(iftSuperpixelIGraph *igraph, iftDHeap *Q, double trivialPathCost, iftSet **removedNodes, iftSet **frontierNodes)
{
  while(*removedNodes != NULL) {
    int p = iftRemoveSet(removedNodes);
    if (Q->color[p] == IFT_GRAY)
      iftRemoveDHeapElem(Q, p);
    else
      Q->color[p] = IFT_WHITE;
  }

  while(*frontierNodes != NULL) {
    int p = iftRemoveSet(frontierNodes);
    
    // Frontier set might include removed nodes, skip those
    if (igraph->ann->pathValue[p] == trivialPathCost)
      continue; 

    // Skip nodes already on the queue and frontier set duplicates
    if (Q->color[p] == IFT_GRAY)
      continue; 

    iftInsertDHeap(Q, p);
  }
}

void iftRISFTreeRemoval(iftSuperpixelIGraph *igraph, iftDHeap *Q, iftSet *markedNodes, double trivialPathCost)
{
  const char *funName = "iftRISFTreeRemoval";
  iftSet *removedNodes = NULL;
  iftSet *frontierNodes = NULL;
  switch (igraph->type) {
    case IMPLICIT:
      iftForestTreeRemovalImplicitAdj(igraph, markedNodes,
          trivialPathCost, &removedNodes, &frontierNodes);
      break;
    case EXPLICIT:
      iftForestTreeRemovalExplicitAdj(igraph->ann, markedNodes,
          igraph->nodeAdj, trivialPathCost,
          &removedNodes, &frontierNodes);
      break;
    case COMPLETE:
      // TODO
    default:
      iftError("Invalid graph type.", funName);
  }
  iftRISFConsumeRemovedAndFrontierNodes(igraph, Q, 
      trivialPathCost, &removedNodes, &frontierNodes);
}

void iftRISFRemoveSmallSuperpixels(iftSuperpixelIGraph *igraph, iftDHeap *Q, float threshold, double trivialPathCost)
{
  // Compute superpixel sizes and store their seeds
  int nSuperpixels = igraph->ann->nLabels;
  int *spSize = iftAllocIntArray(nSuperpixels);
  int *labelSeed = iftAllocIntArray(nSuperpixels);
  for (int p = 0; p < igraph->refImg->n; ++p) {
    int node = igraph->refImg->val[p] - 1;
    if (node < 0)
      continue;
    int label = igraph->ann->label[node] - 1;
    if (label < 0)
      continue;
    spSize[label] += 1;
    labelSeed[label] = igraph->ann->root[node];
  }

  // Create set with superpixels below threshold
  float meanSpSize = (float)igraph->refImg->n / (float)nSuperpixels;
  iftSet *markedNodes = NULL;
  for (int s = 0; s < nSuperpixels; ++s)
    if (spSize[s] < threshold * meanSpSize)
      iftInsertSet(&markedNodes, labelSeed[s]);

  // Remove trees with update seeds and re-initialize
  //   queue accordingly.
  iftRISFTreeRemoval(igraph, Q, markedNodes, trivialPathCost);

  free(spSize);
  free(labelSeed);
}

void iftRemoveLabelGaps(int *labelArray, int size)
{
  assert(labelArray != NULL);
  assert(size > 0);

  // Find maximum label value
  int maxLabel = 0;
  for (int i = 0; i < size; ++i)
    if (labelArray[i] > maxLabel)
      maxLabel = labelArray[i];
  assert(maxLabel > 0);

  // Find number of nodes with a given label
  int *nodeCount = iftAllocIntArray(maxLabel);
  for (int i = 0; i < size; ++i)
    if (labelArray[i] > 0)
      nodeCount[labelArray[i] - 1] += 1;

  // Create a conversion table for non-empty label groups
  int *newLabelTb = iftAllocIntArray(maxLabel);
  int newLabel = 1;
  for (int i = 0; i < maxLabel; ++i)
    if (nodeCount[i] > 0)
      newLabelTb[i] = newLabel++;

  // Convert old labels
  for (int i = 0; i < size; ++i)
    labelArray[i] = newLabelTb[labelArray[i] - 1];

  free(nodeCount);
  free(newLabelTb);
}

iftSet * iftRISFSampling(const iftSuperpixelIGraph *igraph, const iftISFOptions *opt, const iftImage *img)
{
  const char *functionName = "iftRISFSampling";
  assert(igraph != NULL);
  assert(opt != NULL);
  assert(iftIsDomainEqual(img, igraph->refImg));

  // Check if specific selection method is even needed
  if (opt->nSuperpixels >= igraph->nNodes)
    return iftRISFAllNodesSampling(igraph);

  // Retrieve seed set from appropriate function
  iftSet *seeds = NULL;
  switch(opt->sampling) {
    case IFT_ISF_GRID_SAMPLING:
      seeds = iftRISFGridSampling(igraph, img, opt->nSuperpixels);
      break;
    case IFT_ISF_MIXED_SAMPLING:
      seeds = iftRISFMixedSampling(igraph, img, opt->nSuperpixels);
      break;
    case IFT_ISF_RANDOM_SAMPLING:
      seeds = iftRISFRandomSampling(igraph, opt->nSuperpixels);
      break;
    case IFT_ISF_GEODESIC_SAMPLING:
      seeds = iftRISFGeodesicSampling(igraph, img, opt->nSuperpixels);
      break;
    case IFT_ISF_LATTICE_SAMPLING:
      seeds = iftRISFGridOnMaskSampling(igraph, opt->nSuperpixels);
      break;
    case IFT_ISF_HEXGRID_SAMPLING:
      seeds = iftRISFHexGridSampling(igraph, opt->nSuperpixels);
      break;
    default:
      iftError("Sampling strategy not selected.", functionName);
  };

  assert(seeds != NULL);
  return seeds;
}

iftSet * iftRISFGridSampling(const iftSuperpixelIGraph *igraph, const iftImage *img, int nSamples)
{
  assert(igraph != NULL);
  assert(img != NULL);
  assert(iftIsDomainEqual(img, igraph->refImg));
  assert(nSamples > 0 && nSamples < igraph->nNodes);
  
  iftSet *samples = NULL;
  bool is3D = iftIs3DImage(igraph->refImg);
  //iftAdjRel *A = is3D ? iftSpheric(sqrtf(3.0)) : iftCircular(sqrtf(2.0));

  if (igraph->type == EXPLICIT) {
    // Actual grid sampling
    iftBoundingBox minBounds = iftMinBoundingBox(igraph->refImg, NULL);
    int xSize = minBounds.end.x - minBounds.begin.x + 1;
    int ySize = minBounds.end.y - minBounds.begin.y + 1;
    int zSize = minBounds.end.z - minBounds.begin.z + 1;
    float superpixelSize = 0.5 + (float) (xSize * ySize * zSize) /
      (float) nSamples;

    int step = is3D ? (int) ceil(powf(superpixelSize, 1.0/3.0)) :
      (int) ceil(sqrtf(superpixelSize));

    if (step < 1)
      iftError("Number of samples is too high.", "iftRISFGridSampling");

    iftVoxel uBegin, uEnd; 
    uBegin.x = minBounds.begin.x + (step/2);
    uBegin.y = minBounds.begin.y + (step/2);
    uBegin.z = is3D ?  minBounds.begin.z + (step/2) : 0;
    uEnd.x = minBounds.end.x - (step/2) + 1;
    uEnd.y = minBounds.end.y - (step/2) + 1;
    uEnd.z = is3D ? minBounds.end.z - (step/2) + 1 : 1;
    iftVoxel u;
    for (u.z = uBegin.z; u.z < uEnd.z; u.z += step) {
      for (u.y = uBegin.y; u.y < uEnd.y; u.y += step) {
        for (u.x = uBegin.x; u.x < uEnd.x; u.x += step) {
          int p = iftGetVoxelIndex(igraph->refImg, u);
          if (igraph->refImg->val[p] == 0)
            continue;

          // TODO Add local low grad movement
          /*
           * for (int i = 0; i < A->n; ++i) {
           *   iftVoxel v = iftGetAdjacentVoxel(A, u, i);
           *   // Note: Always a valid voxel
           *   int t = iftGetVoxelIndex(
           * }
           */

          iftInsertSet(&samples, igraph->refImg->val[p] - 1);
        }
      }
    }
  } else {
    // For now we redirect to the old implementation for pixels
    // TODO Delete this part once local gradient computation is done
    iftMImage *mimg; 
    if (iftIsColorImage(img))
      mimg = iftImageToMImage(img, LABNorm_CSPACE);
    else
      mimg = iftImageToMImage(img, GRAY_CSPACE);  

    // Image mask for iftGridSampling
    iftImage *tmp = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
    for (int i = 0; i < tmp->n; ++i)
      tmp->val[i] = (igraph->refImg->val[i] <= 0) ? 0 : 255;

    iftImage *grid = iftGridSampling(mimg, tmp, nSamples);

    // Retrieve nodes from the sampling label image
    for (int p = 0; p < grid->n; ++p)
      if (grid->val[p] != 0 && igraph->refImg->val[p] != 0)
        iftInsertSet(&samples, igraph->refImg->val[p] - 1);

    iftDestroyMImage(&mimg);
    iftDestroyImage(&grid);
    iftDestroyImage(&tmp);
  }

  assert(samples != NULL);
  return samples;
}

iftSet * iftRISFMixedSampling(const iftSuperpixelIGraph *igraph, const iftImage *img, int nSamples)
{
  assert(igraph != NULL);
  assert(img != NULL);
  assert(iftIsDomainEqual(img, igraph->refImg));
  assert(nSamples > 0 && nSamples < igraph->nNodes);

  iftMImage *mimg; 
  if (iftIsColorImage(img))
    mimg = iftImageToMImage(img, LABNorm_CSPACE);
  else
    mimg = iftImageToMImage(img, GRAY_CSPACE);  

  // Full image mask for iftAltMixedSampling
  iftImage *tmp = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);
  for (int i = 0; i < tmp->n; ++i)
    tmp->val[i] = 1;
  /* iftAltMixedSampling only considers pixel features, a better
   *  version adapted for superpixels might be implemented later. */
  iftImage *grid = iftAltMixedSampling(mimg, tmp, nSamples);

  // Retrieve nodes from the sampling label image
  // TODO create a function for this step and merge with grid sampling
  iftSet *samples = NULL;
  for (int p = 0; p < grid->n; ++p) {
    if (grid->val[p] == 0)
      continue;
    // Check if it is a pixel graph (which does not use refImg)
    if (igraph->nNodes == igraph->refImg->n)
      iftInsertSet(&samples, p);
    else
      iftInsertSet(&samples, igraph->refImg->val[p] - 1);
  }

  iftDestroyMImage(&mimg);
  iftDestroyImage(&grid);
  iftDestroyImage(&tmp);
  assert(samples != NULL);
  return samples;
}

iftSet * iftRISFRandomSampling(const iftSuperpixelIGraph *igraph, int nSamples)
{
  assert(igraph != NULL);
  assert(nSamples > 0 && nSamples < igraph->nNodes);

  // Build random selection array with all nodes
  int *idxArray = iftAllocIntArray(igraph->nNodes);
  for (int i = 0; i < igraph->nNodes; ++i)
    idxArray[i] = i;

  // Randomly pick the desired number of samples
  iftSet *samples = NULL;
  for (int i = 0; i < nSamples; ++i) {
    int high = nSamples - i - 1;
    int num = iftRandomInteger(0, high);
    iftInsertSet(&samples, idxArray[num]);
    iftSwap(idxArray[num], idxArray[high]); 
  }

  free(idxArray);

  return samples;
}

iftSet * iftRISFAllNodesSampling(const iftSuperpixelIGraph *igraph)
{
  assert(igraph != NULL);

  iftSet *samples = NULL;
  for (int node = 0; node < igraph->nNodes; ++node)
    iftInsertSet(&samples, node);

  return samples;
}

iftSet * iftRISFGeodesicSampling(const iftSuperpixelIGraph *igraph, const iftImage *img, int nSamples)
{
  const char * funName = "iftRISFGeodesicSampling";
  const double trivialPathCost = IFT_INFINITY_DBL;
  const double seedPathCost = 0.0;
  assert(igraph != NULL);
  assert(nSamples > 0 && nSamples < igraph->nNodes);

  // Builds a custom superpixel graph off the original one
  iftSuperpixelIGraph *geograph = iftBuildGeodesicSamplingIGraph(igraph, img);

  iftSet *seeds = NULL;

  // Priority queue init
  // Due to DIFT processing, we need to efficiently keep track of
  //  maximum cost node as well. In a standard IFT it would be the
  //  last node removed from the queue.
  iftDHeap *Q = iftCreateDHeap(geograph->nNodes, geograph->ann->pathValue);
  iftDHeap *R = iftCreateDHeap(geograph->nNodes, geograph->ann->pathValue);
  R->removal_policy = MAXVALUE; 

  // Run an IFT for each new seed
  for (int iter = 0; iter < nSamples; ++iter) {
    // Init seeds
    iftGeodesicSamplingSeedInit(geograph, Q, R, seeds,
        trivialPathCost, seedPathCost, iter);

    // Run IFT
    switch(geograph->type) {
      case IMPLICIT:
        while (!iftEmptyDHeap(Q)) {
          int s = iftRemoveDHeap(Q);
          iftVoxel u = { 
            iftMatrixRowPointer(geograph->pos, s)[0],
            iftMatrixRowPointer(geograph->pos, s)[1],
            (geograph->pos->ncols > 2) ?
              iftMatrixRowPointer(geograph->pos, s)[2] : 0 };

          for (int i = 1; i < geograph->A->n; ++i) {
            iftVoxel v = iftGetAdjacentVoxel(geograph->A, u, i);
            if (!iftValidVoxel(geograph->refImg, v))
              continue;

            int t = iftGetVoxelIndex(geograph->refImg, v); 
            t = geograph->refImg->val[t] - 1;
            if (t < 0)
              continue;

            iftGeodesicSamplingNeighborProcessing(geograph, Q, R, s, t);
          }
        }
        break;
      case EXPLICIT:
        while (!iftEmptyDHeap(Q)) {
          int s = iftRemoveDHeap(Q);
          for (iftSet *adj = geograph->nodeAdj[s]; adj != NULL; adj = adj->next) {
            int t = adj->elem;
            iftGeodesicSamplingNeighborProcessing(geograph, Q, R, s, t);
          }
        }
        break;
      case COMPLETE:
        while (!iftEmptyDHeap(Q)) {
          int s = iftRemoveDHeap(Q);
          for (int t = 0; t < geograph->nNodes; ++t)
            if (s != t)
              iftGeodesicSamplingNeighborProcessing(geograph, 
                  Q, R, s, t);
        }
        break;
      default:
        iftError("Invalid graph type.", funName);
    }
    iftInsertSet(&seeds, iftRemoveDHeap(R));
    iftResetDHeap(Q);
  }
  
  iftDestroyDHeap(&Q);
  iftDestroySuperpixelIGraph(&geograph);

  // Re-order to chronological order
  iftInvertSet(&seeds);

  // TMP
  /*iftImage *vis = iftCreateImage(igraph->refImg->xsize, igraph->refImg->ysize, igraph->refImg->zsize);
  for (iftSet *s = seeds; s != NULL; s = s->next)
    for (int i = 0; i < vis->n; ++i)
      if (igraph->refImg->val[i] - 1 == s->elem)
        vis->val[i] = 255;
  iftWriteImageByExt(vis, "lul.png");
  iftDestroyImage(&vis);*/

  return seeds;
}

iftSet * iftRISFGridOnMaskSampling(const iftSuperpixelIGraph *igraph, int nSamples)
{
  assert(igraph != NULL);
  assert(nSamples > 0 && nSamples < igraph->nNodes);

  float rad = iftEstimateGridOnMaskSamplingRadius(igraph->refImg, -1, nSamples);
  iftIntArray *seeds = iftGridSamplingOnMask(igraph->refImg, rad, -1, 0);
  iftSet *samples = NULL;
  for (int i = 0; i < seeds->n; ++i)
    iftInsertSet(&samples, igraph->refImg->val[seeds->val[i]] - 1);

  iftDestroyIntArray(&seeds);

  assert(samples != NULL);
  return samples;
}

iftSet * iftRISFHexGridSampling(const iftSuperpixelIGraph *igraph, int nSamples)
{
  assert(igraph != NULL);
  assert(nSamples > 0 && nSamples < igraph->nNodes);
  
  iftSet *samples = NULL;
  bool is3D = iftIs3DImage(igraph->refImg);
  if (is3D)
    iftError("iftRISFHexGridSampling", "Only implemented for 2D");
  //iftAdjRel *A = is3D ? iftSpheric(sqrtf(3.0)) : iftCircular(sqrtf(2.0));

  // Actual grid sampling
  iftBoundingBox minBounds = iftMinBoundingBox(igraph->refImg, NULL);
  int xSize = minBounds.end.x - minBounds.begin.x + 1;
  int ySize = minBounds.end.y - minBounds.begin.y + 1;
  float superpixelSize = 0.5 + (float) (xSize * ySize) / (float) nSamples;

  int xstep = (int) ceil(sqrtf(superpixelSize));
  int ystep = (int) ceil(sqrtf(3.0f) * (float)xstep);

  if (ystep < 1)
    iftError("Number of samples is too high.", "iftRISFHexGridSampling");

  iftVoxel uBegin, uEnd; 
  // First alignment
  uBegin.x = minBounds.begin.x + (xstep/4);
  uBegin.y = minBounds.begin.y + (ystep/4);
  uBegin.z = 0;
  uEnd.x = minBounds.end.x - (xstep/4) + 1;
  uEnd.y = minBounds.end.y - (ystep/4) + 1;
  uEnd.z = 1;
  iftVoxel u;
  u.z = 0;
  for (u.y = uBegin.y; u.y < uEnd.y; u.y += ystep) {
    for (u.x = uBegin.x; u.x < uEnd.x; u.x += xstep) {
      int p = iftGetVoxelIndex(igraph->refImg, u);
      if (igraph->refImg->val[p] == 0)
        continue;
      iftInsertSet(&samples, igraph->refImg->val[p] - 1);
    }
  }
  // Second alignment
  uBegin.x = minBounds.begin.x + ((xstep/4)*3);
  uBegin.y = minBounds.begin.y + ((ystep/4)*3);
  uEnd.x = minBounds.end.x - ((xstep/4)*3) + 1;
  uEnd.y = minBounds.end.y - ((ystep/4)*3) + 1;
  for (u.y = uBegin.y; u.y < uEnd.y; u.y += ystep) {
    for (u.x = uBegin.x; u.x < uEnd.x; u.x += xstep) {
      int p = iftGetVoxelIndex(igraph->refImg, u);
      if (igraph->refImg->val[p] == 0)
        continue;
      iftInsertSet(&samples, igraph->refImg->val[p] - 1);
    }
  }

  assert(samples != NULL);
  return samples;
}

void iftRISFRemoveDuplicateSeeds(iftSet **seeds, int nNodes)
{
  // TODO Avoid the need for nNodes and make the method more general
  int *flag = iftAllocIntArray(nNodes);
  iftSet *tmpSeeds = NULL;
  for (iftSet *s = *seeds; s != NULL; s = s->next) {
    if (flag[s->elem] == 0) {
      flag[s->elem] = 1;
      iftInsertSet(&tmpSeeds, s->elem);
    }
  }

  free(flag);
  iftDestroySet(seeds);
  *seeds = tmpSeeds;
}

iftSuperpixelIGraph * iftBuildGeodesicSamplingIGraph(const iftSuperpixelIGraph *igraph, const iftImage *img)
{
  iftSuperpixelIGraph *geograph = iftCopySuperpixelIGraph(igraph);

  // We compute color mean feats again even if igraph already has
  //   them because there is no convenient way to check for
  //   alternative superpixel features (e.g. color histogram).
  iftMatrix *pFeats = NULL;
  if (iftIsColorImage(img))
    pFeats = iftComputeSuperpixelFeaturesByColorSpaceMean(geograph->refImg, img, LABNorm_CSPACE);  
  else
    pFeats = iftComputeSuperpixelFeaturesByColorSpaceMean(geograph->refImg, img, GRAY_CSPACE);  
  iftSetSuperpixelIGraphFeatures(geograph, pFeats, iftDistance1, 0.0, 0.0);

  iftDestroyMatrix(&pFeats);
  return geograph;
}

void iftGeodesicSamplingSeedInit(iftSuperpixelIGraph *igraph, iftDHeap *Q, iftDHeap *R, iftSet *seeds, double trivialPathCost, double seedPathCost, int iter)
{
  if (iter == 0 || iter == 1) {
    // Standard IFT initialization with a single seed
    iftInitForestAnnotationInfo(igraph->ann, trivialPathCost, false);

    int firstSeed;
    if (iter == 0) {
      // First one is arbitrary and we ignore it on later iterations
      firstSeed = igraph->nNodes / 2;
    } else {
      // We start again after removing arbitrary seed "zero"
      iftResetDHeap(R);
      firstSeed = seeds->elem;
    }

    iftInitForestAnnotationSingleSeed(igraph->ann, firstSeed, seedPathCost, false);
    iftInsertDHeap(Q, firstSeed);
  /*if (iter == 0) {
    iftSet *borderSeeds = NULL;
    iftInitForestAnnotationInfo(igraph->ann, trivialPathCost, false);
    iftVoxel u = {0,0,0};
    int p, node;
    for (u.x = 0; u.x < igraph->refImg->xsize; ++(u.x)) {
      u.y = 0;
      p = iftGetVoxelIndex(igraph->refImg, u);
      node = igraph->refImg->val[p] - 1;
      iftUnionSetElem(&borderSeeds, node);
      u.y = igraph->refImg->ysize - 1;
      p = iftGetVoxelIndex(igraph->refImg, u);
      node = igraph->refImg->val[p] - 1;
      iftUnionSetElem(&borderSeeds, node);
    }
    for (u.y = 0; u.y < igraph->refImg->ysize; ++(u.y)) {
      u.x = 0;
      p = iftGetVoxelIndex(igraph->refImg, u);
      node = igraph->refImg->val[p] - 1;
      iftUnionSetElem(&borderSeeds, node);
      u.x = igraph->refImg->xsize - 1;
      p = iftGetVoxelIndex(igraph->refImg, u);
      node = igraph->refImg->val[p] - 1;
      iftUnionSetElem(&borderSeeds, node);
    }
    while (borderSeeds != NULL) {
      int node = iftRemoveSet(&borderSeeds);
      if (node >= 0 && node < igraph->nNodes) {
        iftInitForestAnnotationSingleSeed(igraph->ann, node, seedPathCost, false);
        iftInsertDHeap(Q, node);
      }
    }*/
  } else {
    // Add last chosen seed on remaining iterations
    int lastSeed = seeds->elem;
    iftInitForestAnnotationSingleSeed(igraph->ann, lastSeed, seedPathCost, false);
    iftInsertDHeap(Q, lastSeed);
  }
}

void iftGeodesicSamplingNeighborProcessing(iftSuperpixelIGraph *igraph, iftDHeap *Q, iftDHeap *R, int s, int t)
{
  if (Q->color[t] == IFT_BLACK)
    return;

  float pCost = igraph->distFun(
      iftMatrixRowPointer(igraph->feats, igraph->ann->root[s]),
      iftMatrixRowPointer(igraph->feats, t),
      igraph->distAlpha,
      igraph->feats->ncols); 
  /*double distCost = iftFeatDistance(
    iftMatrixRowPointer(igraph->pos, igraph->ann->root[s]),
    iftMatrixRowPointer(igraph->pos, t),
    igraph->pos->ncols);*/

  //double cost = iftMax(igraph->ann->pathValue[s], pCost);
  // Segmentation function
  //double cost = igraph->ann->pathValue[s] +  0.5f * pCost + distCost;
  // Original geodesic
  double cost = igraph->ann->pathValue[s] +  pCost;
  // Regular distance
  //double cost = igraph->ann->pathValue[s] + distCost;

  double prevCost = igraph->ann->pathValue[t];
  if (cost < prevCost) {
    // Better cost, standard IFT behavior
    igraph->ann->pathValue[t] = cost;
    igraph->ann->root[t] = igraph->ann->root[s];
    igraph->ann->pred[t] = s;
    if (Q->color[t] == IFT_GRAY)
      iftGoUpDHeap(Q, Q->pos[t]);
    else
      iftInsertDHeap(Q, t);
    if (R->color[t] == IFT_GRAY)
      iftGoDownDHeap(R, R->pos[t]);
    else if (R->color[t] == IFT_WHITE)
      iftInsertDHeap(R, t);
  } else {
    if (igraph->ann->pred[t] == s) {
      // DIFTmod with simplified tie processing
      int sRoot = igraph->ann->root[s];
      int tRoot = igraph->ann->root[t];
      if (cost > prevCost || sRoot != tRoot)
        iftGeodesicSamplingRemoveInconsistentSubTree(igraph, Q, R, t);
    }
  }
}

void iftGeodesicSamplingRemoveInconsistentSubTree(iftSuperpixelIGraph *igraph, iftDHeap *Q, iftDHeap *R, int t)
{
  const char *funName = "iftGeodesicSamplingRemoveInconsistentSubTree";
  // We divide by 2 to avoid conflict with iftRemoveDHeapElem on R
  const double trivialPathCost = IFT_INFINITY_DBL / 2.0;

  // TODO Remove redundancy with iftRISFRemoveInconsistentSubTree
  iftSet *removedNodes = NULL;
  iftSet *frontierNodes = NULL;
  switch (igraph->type) {
    case IMPLICIT:
      iftRemoveSubTreeImplicitAdj(igraph, t,
          trivialPathCost, &removedNodes, &frontierNodes);
      break;
    case EXPLICIT:
      iftRemoveSubTreeExplicitAdj(igraph->ann, t,
          igraph->nodeAdj, trivialPathCost,
          &removedNodes, &frontierNodes);
      break;
    case COMPLETE:
      // TODO
    default:
      iftError("Invalid graph type.", funName);
  }
  for (iftSet *s = removedNodes; s != NULL; s = s->next)
    if (R->color[s->elem] == IFT_GRAY)
      iftRemoveDHeapElem(R, s->elem);
  iftRISFConsumeRemovedAndFrontierNodes(igraph, Q, 
      trivialPathCost, &removedNodes, &frontierNodes);
}

// ---------- Public Methods Implementation ----------

// iftISFOptions begin

iftISFOptions * iftInitISFOptions(int nSuperpixels, int nIters, int nSmoothIters, float smallSpThreshold, iftISFSampling sampling, iftISFPathCost pathCost, iftISFSeedRecomp seedRecomp)
{
  iftISFOptions *opt = (iftISFOptions *) iftAlloc(1, sizeof(*opt));
  opt->nSuperpixels = nSuperpixels;
  opt->nIters = nIters;
  opt->nSmoothIters = nSmoothIters;
  opt->smallSpThreshold = smallSpThreshold;
  opt->sampling = sampling; 
  opt->pathCost = pathCost; 
  opt->seedRecomp = seedRecomp; 
  
  return opt;
}

void iftDestroyISFOptions(iftISFOptions **options)
{
  if (options == NULL)
    return;

  iftISFOptions *tmp = *options;
  if (tmp == NULL)
    return;

  free(tmp);
  *options = NULL;
}

// iftISFOptions end, iftForestAnnotationInfo begin

iftForestAnnotationInfo * iftCreateForestAnnotationInfo(int n)
{
  assert(n > 0);

  iftForestAnnotationInfo *ann = iftAlloc(1,sizeof(*ann));
  ann->n = n;
  ann->nLabels = 0;
  ann->label = iftAllocIntArray(n);
  ann->marker = iftAllocIntArray(n);
  ann->root = iftAllocIntArray(n);
  ann->pred = iftAllocIntArray(n);
  ann->pathValue = iftAllocDoubleArray(n);

  return ann;
}

void iftDestroyForestAnnotationInfo(iftForestAnnotationInfo **ann)
{
  if (ann == NULL || *ann == NULL)
    return;

  iftForestAnnotationInfo *tmp = *ann;
  free(tmp->label);
  free(tmp->marker);
  free(tmp->root);
  free(tmp->pred);
  free(tmp->pathValue);
  free(tmp);

  *ann = NULL;
}

iftForestAnnotationInfo * iftCopyForestAnnotationInfo(const iftForestAnnotationInfo *ann)
{
  assert(ann != NULL);

  iftForestAnnotationInfo *res = iftAlloc(1,sizeof(*ann));
  res->n = ann->n;
  res->nLabels = ann->nLabels;
  res->label = iftAllocIntArray(res->n);
  iftCopyIntArray(res->label, ann->label, res->n);
  res->marker = iftAllocIntArray(res->n);
  iftCopyIntArray(res->marker, ann->marker, res->n);
  res->root = iftAllocIntArray(res->n);
  iftCopyIntArray(res->root, ann->root, res->n);
  res->pred = iftAllocIntArray(res->n);
  iftCopyIntArray(res->root, ann->root, res->n);
  res->pathValue = iftAllocDoubleArray(res->n);
  iftCopyDoubleArray(res->pathValue, ann->pathValue, res->n);

  return res;
}

void iftInitForestAnnotationInfo(iftForestAnnotationInfo *ann, const double trivialPathCost, bool resetLabelData)
{
  assert(ann != NULL);
  if (resetLabelData)
    ann->nLabels = 0;

  for (int i = 0; i < ann->n; ++i) {
    if (resetLabelData)
      ann->label[i] = 0;
    ann->marker[i] = 0;
    ann->root[i] = IFT_NIL;
    ann->pred[i] = IFT_NIL;
    ann->pathValue[i] = trivialPathCost;
  }
}

void iftInitForestAnnotationSingleSeed(iftForestAnnotationInfo *ann, int seed, const double seedPathCost, const bool createNewLabels)
{
  assert(ann != NULL);
  assert(seed >= 0 && seed < ann->n);

  if (createNewLabels)
    ann->label[seed] = ++(ann->nLabels);
  ann->root[seed] = seed;
  ann->pred[seed] = IFT_NIL;
  ann->pathValue[seed] = seedPathCost;
}

void iftInitForestAnnotationSeeds(iftForestAnnotationInfo *ann, iftSet *seeds, const double seedPathCost, const bool createNewLabels)
{
  assert(ann != NULL && seeds != NULL);

  iftSet *S = iftSetCopy(seeds);
  while(S != NULL) {
    iftInitForestAnnotationSingleSeed(
        ann, iftRemoveSet(&S), seedPathCost, createNewLabels);
  }
}

// iftForestAnnotationInfo end, iftSuperpixelIGraph begin

iftImage * iftSuperpixelSegmentationByISF(iftImage *img, iftMImage *mimg, iftISFOptions *options)
{
  assert(mimg != NULL && options != NULL);

  iftAdjRel *A = iftIs3DMImage(mimg) ? iftSpheric(1.0) : iftCircular(1.0);

  // Select entire image for both sampling and igraph
  iftImage *mask = iftSelectImageDomain(mimg->xsize, mimg->ysize, mimg->zsize);

  // Seed sampling step
  iftImage *seeds = iftISFPixelSampling(options->sampling, mimg, mask, options->nSuperpixels);

  // IFT step
  iftIGraph *igraph = iftImplicitIGraph(mimg, mask, A);
  iftISFRunIFT(options, igraph, seeds);

  // Smooth regions in the label map of igraph
  if (options->nSmoothIters > 0){
    iftIGraphSetWeightForRegionSmoothing(igraph, img);
    iftIGraphSmoothRegions(igraph, options->nSmoothIters);
  }

  // Convert into final segmentation
  iftImage *label = iftIGraphLabel(igraph);

  // Clean up
  iftDestroyImage(&mask);
  iftDestroyImage(&seeds);
  iftDestroyIGraph(&igraph);
  iftDestroyAdjRel(&A);

  return label;
}

iftSuperpixelIGraph *iftInitSuperpixelIGraph(const iftImage *refImg)
{
  assert(refImg != NULL);
  int nNodes = iftMaximumValue(refImg);
  assert(nNodes > 0);

  iftSuperpixelIGraph *igraph = iftAlloc(1, sizeof(*igraph));
  igraph->nNodes = nNodes;
  igraph->nodeAdj = iftAlloc(nNodes, sizeof(*(igraph->nodeAdj)));
  igraph->feats = NULL;
  igraph->seedFeats = NULL;
  igraph->distFun = iftDistance1;
  igraph->distAlpha = NULL;
  igraph->pos = iftComputeSuperpixelFeaturesByGeometricCenter(refImg);
  igraph->seedPos = NULL;
  igraph->alpha = 0.0;
  igraph->beta = 1.0;
  igraph->refImg = iftCopyImage(refImg);
  igraph->spSize = iftAllocLongLongIntArray(nNodes);
  for (int p = 0; p < refImg->n; ++p) {
    int label = refImg->val[p] - 1;
    if (label >= 0)
      igraph->spSize[label] += 1;
  }
  igraph->ann = iftCreateForestAnnotationInfo(nNodes);
  igraph->A = NULL;
  igraph->type = COMPLETE;
  
  return igraph;
}

void iftDestroySuperpixelIGraph(iftSuperpixelIGraph **igraph)
{
  if (igraph == NULL || *igraph == NULL)
    return;

  iftSuperpixelIGraph *tmp = *igraph;
  for (int i = 0; i < tmp->nNodes; ++i)
    iftDestroySet(&(tmp->nodeAdj[i]));
  free(tmp->nodeAdj);
  iftDestroyMatrix(&(tmp->feats));
  iftDestroyMatrix(&(tmp->seedFeats));
  free(tmp->distAlpha);
  iftDestroyMatrix(&(tmp->pos));
  iftDestroyMatrix(&(tmp->seedPos));
  iftDestroyImage(&(tmp->refImg));
  free(tmp->spSize);
  iftDestroyForestAnnotationInfo(&(tmp->ann));
  iftDestroyAdjRel(&(tmp->A));
  
  free(tmp);
  *igraph = NULL;
}

iftSuperpixelIGraph *iftCopySuperpixelIGraph(const iftSuperpixelIGraph *igraph)
{
  assert(igraph != NULL);
  assert(igraph->nodeAdj != NULL);
  assert(igraph->refImg != NULL);
  assert(igraph->ann != NULL);

  iftSuperpixelIGraph *res = iftAlloc(1, sizeof(*res));
  res->nNodes = igraph->nNodes;
  res->nodeAdj = iftAlloc(res->nNodes, sizeof(*(res->nodeAdj)));
  for (int s = 0; s < res->nNodes; ++s)
    res->nodeAdj[s] = iftSetCopy(igraph->nodeAdj[s]);

  if (igraph->feats != NULL) {
    res->feats = iftCopyMatrix(igraph->feats);
    res->distAlpha = iftAllocFloatArray(igraph->feats->ncols);
    iftCopyFloatArray(res->distAlpha, igraph->distAlpha,
        igraph->feats->ncols);
  } else {
    res->feats = NULL;
    res->distAlpha = NULL;
  }

  if (igraph->seedFeats != NULL)
    res->seedFeats = iftCopyMatrix(igraph->seedFeats);
  else
    res->seedFeats = NULL;

  res->distFun = igraph->distFun;
  res->pos = iftCopyMatrix(igraph->pos);

  if (igraph->seedPos != NULL)
    res->seedPos = iftCopyMatrix(igraph->seedPos);
  else
    res->seedPos = NULL;

  res->alpha = igraph->alpha;
  res->beta = igraph->beta;
  res->refImg = iftCopyImage(igraph->refImg);
  res->spSize = iftAllocLongLongIntArray(igraph->nNodes);
  iftCopyLongLongIntArray(res->spSize, igraph->spSize, igraph->nNodes);
  res->ann = iftCopyForestAnnotationInfo(igraph->ann);
  
  if (igraph->A != NULL)
    res->A = iftCopyAdjacency(igraph->A);
  else
    res->A = NULL;

  res->type = igraph->type;
  
  return res;
}

void iftSetSuperpixelIGraphImplicitAdjacency(iftSuperpixelIGraph *igraph, iftAdjRel *A)
{
  assert(igraph != NULL);
  assert(A != NULL);
  //assert(igraph->nNodes == igraph->refImg->n);

  igraph->A = iftCopyAdjacency(A);
  igraph->type = IMPLICIT;
}

void iftSetSuperpixelIGraphExplicitRAG(iftSuperpixelIGraph *igraph, const iftAdjRel *A)
{
  assert(igraph != NULL);
  assert(A != NULL);

  igraph->type = EXPLICIT;

  // Find neighboring superpixels and update graph edges accordingly
  iftImage *img = igraph->refImg; // shorthand
  for (int p = 0; p < img->n; ++p) {
    if (img->val[p] == 0)
      continue;
    iftVoxel u = iftGetVoxelCoord(img, p);
    for (int i = 1; i < A->n; ++i) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, i);
      if (!iftValidVoxel(img, v))
        continue;
      int q = iftGetVoxelIndex(img, v);
      if (img->val[q] == 0 || img->val[p] == img->val[q])
        continue;
      /* We assume that the number of neighbors << number of nodes
       *  so that traversing the entire node adjacency each time
       *  is not a problem. */
      iftUnionSetElem(&(igraph->nodeAdj[img->val[p]-1]),
          img->val[q]-1);
    }
  }
}

void iftSetSuperpixelIGraphFeatures(iftSuperpixelIGraph *igraph, const iftMatrix *superpixelFeatures, iftArcWeightFun distFun, float alpha, float beta)
{
  // Clear any old data
  iftDestroyMatrix(&(igraph->feats));
  free(igraph->distAlpha);

  // Assing new features
  igraph->feats = iftCopyMatrix(superpixelFeatures);
  igraph->distFun = distFun;
  igraph->distAlpha = iftAllocFloatArray(igraph->feats->ncols);
  for (int i = 0; i < igraph->feats->ncols; ++i)
    igraph->distAlpha[i] = 1.0;
  igraph->alpha = alpha;
  igraph->beta = beta;
}

iftImage * iftSuperpixelSegmentationByRISF(iftSuperpixelIGraph *igraph, const iftISFOptions *opt, const iftImage *img)
{
  const char *funName = "iftSuperpixelSegmentationByRISF";
  const double trivialPathCost = IFT_INFINITY_DBL;
  const double seedPathCost = 0.0;
  assert(igraph != NULL);
  assert(igraph->feats != NULL);
  assert(opt != NULL);
  assert(img != NULL);
  assert(iftIsDomainEqual(img, igraph->refImg));

  // Sampling step
  iftSet *seeds = iftRISFSampling(igraph, opt, img);
  iftRISFRemoveDuplicateSeeds(&seeds, igraph->nNodes);

  iftDHeap *Q = iftCreateDHeap(igraph->nNodes, igraph->ann->pathValue);

  // Initialization including label data
  iftInitForestAnnotationInfo(igraph->ann, trivialPathCost, true);
  iftInitForestAnnotationSeeds(igraph->ann, seeds, seedPathCost, true);

  // Main ISF loop
  int nIters = opt->smallSpThreshold > 0.0 ? opt->nIters + 1 : opt->nIters;
  for (int iter = 0; iter < nIters; ++iter) {
    // Init or update auxiliary features for path computation
    // TODO Make this step more general, for now only MEAN uses it
    if (opt->pathCost == IFT_ISF_ADD_MEAN_PATHCOST) {
      if (iter == 0)
        iftInitSuperpixelIGraphMeanAuxFeat(igraph, seeds);
      else
        iftUpdateSuperpixelIGraphMeanAuxFeat(igraph);
    }

    if (iter == 0) {
      // Prepare for standard IFT
      while (seeds != NULL)
        iftInsertDHeap(Q, iftRemoveSet(&seeds));
    } else if (iter == opt->nIters) {
      // Extra iteration to remove small superpixels
      iftRISFRemoveSmallSuperpixels(igraph, Q,
          opt->smallSpThreshold, trivialPathCost);
    } else {
      // Prepare for DIFT
      // Seed recomputation, tree removal and priority queue update
      iftRISFInitializeDIFT(igraph, opt, Q, trivialPathCost);
    }

    // IFT algorithm
    switch(igraph->type) {
      case IMPLICIT: 
        while (!iftEmptyDHeap(Q)) {
          int s = iftRemoveDHeap(Q);
          iftVoxel u = { 
            iftMatrixRowPointer(igraph->pos, s)[0],
            iftMatrixRowPointer(igraph->pos, s)[1],
            (igraph->pos->ncols > 2) ?
              iftMatrixRowPointer(igraph->pos, s)[2] : 0 };

          for (int i = 1; i < igraph->A->n; ++i) {
            iftVoxel v = iftGetAdjacentVoxel(igraph->A, u, i);
            if (!iftValidVoxel(igraph->refImg, v))
              continue;

            int t = iftGetVoxelIndex(igraph->refImg, v); 
            t = igraph->refImg->val[t] - 1;
            if (t < 0)
              continue;

            iftRISFNeighborProcessingDIFT(igraph, Q, 
                trivialPathCost, s, t, opt);
          }
        }
        break;
      case EXPLICIT:
        while (!iftEmptyDHeap(Q)) {
          int s = iftRemoveDHeap(Q);
          iftSet *adj = igraph->nodeAdj[s];
          for (; adj != NULL; adj = adj->next) {
            int t = adj->elem;
            iftRISFNeighborProcessingDIFT(igraph, Q, 
                trivialPathCost, s, t, opt);
          }
        }
        break;
      case COMPLETE:
        while (!iftEmptyDHeap(Q)) {
          int s = iftRemoveDHeap(Q);
          for (int t = 0; t < igraph->nNodes; ++t) {
            if (s != t)
              iftRISFNeighborProcessingDIFT(igraph, Q, 
                  trivialPathCost, s, t, opt);
          }
        }
        break;
      default:
        iftError("Invalid graph type.", funName);
    }
    // TEMPORARILY COMMENTED TO REPRODUCE OLD RESULTS
    //iftResetDHeap(Q);
  }

  // Extract resulting segmentation label image
  iftRemoveLabelGaps(igraph->ann->label, igraph->nNodes);
  iftImage *res = iftCopyImage(igraph->refImg);
  for (int p = 0; p < res->n; ++p) {
    int node = igraph->refImg->val[p] - 1;
    res->val[p] = igraph->ann->label[node];
  }

  // Clean up
  iftDestroySet(&seeds);
  iftDestroyDHeap(&Q);

  return res;
}

iftImage * iftSuperpixelSegmentationByRISFwithoutDIFT(iftSuperpixelIGraph *igraph, const iftISFOptions *opt, const iftImage *img)
{
  const char *funName = "iftSuperpixelSegmentationByRISF";
  assert(igraph != NULL);
  assert(opt != NULL);
  assert(img != NULL);
  assert(iftIsDomainEqual(img, igraph->refImg));

  // Sampling step
  iftSet *origSeeds = iftRISFSampling(igraph, opt, img);
  iftSet *seeds = iftSetCopy(origSeeds); // Modified through iterations

  iftDHeap *Q = iftCreateDHeap(igraph->nNodes, igraph->ann->pathValue);

  // Initialization including label data
  iftInitForestAnnotationInfo(igraph->ann, IFT_INFINITY_DBL, true);
  iftInitForestAnnotationSeeds(igraph->ann, seeds, 0.0, true);

  // Main ISF loop
  for (int iter = 0; iter < opt->nIters; ++iter) {
    //timer *t1 = iftTic();
    // Init or update auxiliary features for path computation
    // TODO Make this step more general, for now only MEAN uses it
    if (opt->pathCost == IFT_ISF_ADD_MEAN_PATHCOST) {
      if (iter == 0)
        iftInitSuperpixelIGraphMeanAuxFeat(igraph, seeds);
      else
        iftUpdateSuperpixelIGraphMeanAuxFeat(igraph);
    }

    // For subsequent computations, we need to keep label data and
    //  reset IFT information.
    if (iter > 0) {
      // Seed recomputation step
      switch(opt->seedRecomp) {
        case IFT_ISF_MEDOID_SEEDRECOMP:
          seeds = iftRISFMedoidSeedRecomp(igraph, true);
          break;
        case IFT_ISF_CENTROID_SEEDRECOMP:
          // Centroid is equivalent to spatial medoid
          seeds = iftRISFMedoidSeedRecomp(igraph, false);
          break;
        default:
          iftError("Invalid seed recomputation method.", funName);
      }

      // Threshold features
      iftSeedUpdateThreshold(igraph, seeds);

      // IMPORTANT: we currently assume the new seeds set follows the
      //   same ordering from previous iteration. Otherwise anything
      //   indexed from previous iteration seeds (e.g. seedFeats)
      //   will be mixed up.
      iftInitForestAnnotationInfo(igraph->ann, IFT_INFINITY_DBL, true);
      iftInitForestAnnotationSeeds(igraph->ann, seeds, 0.0, true);
      iftResetDHeap(Q);
    }

    // Add seeds to heap
    while (seeds != NULL)
      iftInsertDHeap(Q, iftRemoveSet(&seeds));

    // IFT algorithm
    switch(igraph->type) {
      case IMPLICIT: 
        assert(igraph->nNodes == igraph->refImg->n); // Pixel graph
        while (!iftEmptyDHeap(Q)) {
          int s = iftRemoveDHeap(Q);
          iftVoxel u = iftGetVoxelCoord(igraph->refImg, s);
          for (int i = 1; i < igraph->A->n; ++i) {
            iftVoxel v = iftGetAdjacentVoxel(igraph->A, u, i);
            if (!iftValidVoxel(igraph->refImg, v))
              continue;
            int t = iftGetVoxelIndex(igraph->refImg, v); 
            if (Q->color[t] == IFT_BLACK)
              continue;
            if (iftRISFOfferNewPath(igraph, s, t, opt) > 0) {
              if (Q->color[t] == IFT_GRAY)
                iftGoUpDHeap(Q, Q->pos[t]);
              else
                iftInsertDHeap(Q, t);
            }
          }
        }
        break;
      case EXPLICIT:
        while (!iftEmptyDHeap(Q)) {
          int s = iftRemoveDHeap(Q);
          iftSet *adj = igraph->nodeAdj[s];
          for (; adj != NULL; adj = adj->next) {
            int t = adj->elem;
            if (Q->color[t] == IFT_BLACK)
              continue;
            if (iftRISFOfferNewPath(igraph, s, t, opt) > 0) {
              if (Q->color[t] == IFT_GRAY)
                iftGoUpDHeap(Q, Q->pos[t]);
              else
                iftInsertDHeap(Q, t);
            }
          }
        }
        break;
      case COMPLETE:
        while (!iftEmptyDHeap(Q)) {
          int s = iftRemoveDHeap(Q);
          for (int t = 0; t < igraph->nNodes; ++t) {
            if (s == t || Q->color[t] == IFT_BLACK)
              continue;
            if (iftRISFOfferNewPath(igraph, s, t, opt) > 0) {
              if (Q->color[t] == IFT_GRAY)
                iftGoUpDHeap(Q, Q->pos[t]);
              else
                iftInsertDHeap(Q, t);
            }
          }
        }
        break;
      default:
        iftError("Invalid graph type.", funName);
    }
    //timer* t2 = iftToc();
    //printf("Iter %d time = %fs\n", iter, 0.001*iftCompTime(t1, t2));
  }

  // Smoothing step
  iftRISFSmoothing(igraph, img, opt->nSmoothIters);

  // Extract resulting segmentation label image
  iftImage *res = iftCopyImage(igraph->refImg);
  for (int p = 0; p < res->n; ++p) {
    int node = igraph->refImg->val[p] - 1;
    res->val[p] = igraph->ann->label[node];
  }

  // Clean up
  iftDestroySet(&origSeeds);
  iftDestroySet(&seeds);
  iftDestroyDHeap(&Q);

  return res;
}

iftMatrix * iftComputeSuperpixelFeaturesByGeometricCenter(const iftImage *superpixelLabelMap)
{
  assert(superpixelLabelMap != NULL);

  // Build an MImage containing voxel coordinates as features
  const iftImage *map = superpixelLabelMap; // shorthand
  int dim = iftIs3DImage(map) ? 3 : 2;
  iftMImage *mimg = iftCreateMImage(map->xsize, map->ysize, map->zsize, dim);
  for (int p = 0; p < map->n; ++p) {
    iftVoxel u = iftGetVoxelCoord(map, p); 
    mimg->val[p][0] = (float) u.x;
    mimg->val[p][1] = (float) u.y;
    if (dim == 3)
      mimg->val[p][2] = (float) u.z;
  }

  // Use existing MImage Mean function to obtain superpixel mean
  iftMatrix *mx = iftComputeSuperpixelFeaturesByMImageMean(mimg, map);

  iftDestroyMImage(&mimg);
  return mx;
}

iftMatrix * iftComputeSuperpixelFeaturesByColorSpaceMean(const iftImage *superpixelLabelMap, const iftImage *img, iftColorSpace colorSpace)
{
  iftMImage *mimg = iftImageToMImage(img, colorSpace);
  iftMatrix *res = iftComputeSuperpixelFeaturesByMImageMean(mimg, superpixelLabelMap);
  iftDestroyMImage(&mimg);
  return res;
}

