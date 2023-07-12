#ifndef _IFT_EXPERIMENT_UTILITY_H_
#define _IFT_EXPERIMENT_UTILITY_H_

#include "ift.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Overlays an image with its segmentation borders.
 *
 * @param[in] img The base image.
 * @param[in] segLabelMap Some segmentation of \c img.
 * @param[in] YCbCr Resulting border color.
 *
 * @pre Both images are non-NULL and share the same domain.
 */
iftImage *iftOverlaySegmentationBorders(iftImage *img, iftImage *segLabelMap, iftColor YCbCr);

/**
 * @brief Assigns a segmentation label to each superpixel based on
 *   the groundtruth.
 *
 * @param[in] segLabelMap Some superpixel segmentation.
 * @param[in] gtLabelMap Groundtruth segmentation of objects.
 *
 * @pre Both maps are non-NULL and share the same domain.
 */
iftImage *iftSuperpixelToMajoritySegmentation(const iftImage *segLabelMap, const iftImage *gtLabelMap);

/**
 * @brief Computes ASA (achievable segmentation accuracy).
 *
 * @param[in] segLabelMap Some superpixel segmentation.
 * @param[in] gtLabelMap Groundtruth segmentation of objects.
 *
 * @pre Both are non-NULL and share the same domain.
 */
float iftAchievableSegmentationAccuracy(const iftImage *segLabelMap, const iftImage *gtLabelMap);

/**
 * @brief Computes Boundary Precision.
 *
 * @param[in] segBorders Some superpixel segmentation borders.
 * @param[in] gtBorders Groundtruth segmentation borders.
 * @param[in] toleranceDist Max distance a \c segBorders pixel can be from 
 *                          a \c gtBorders pixel and be considered correct.
 *
 * @pre Both borders are non-NULL and share the same domain.
 * @pre \c toleranceDist is non-negative.
 * @pre \c segBorders is non-empty.
 */
float iftBoundaryPrecision(const iftImage *segBorders, const iftImage *gtBorders, float toleranceDist); 

/**
 * @brief Computes Boundary Fscore.
 *
 * @param[in] segBorders Some superpixel segmentation borders.
 * @param[in] gtBorders Groundtruth segmentation borders.
 * @param[in] toleranceDist Max distance a \c segBorders pixel can be from 
 *                          a \c gtBorders pixel and be considered correct.
 *
 * @pre Both borders are non-NULL and share the same domain.
 * @pre \c toleranceDist is non-negative.
 */
float iftBoundaryFScore(const iftImage *segBorders, const iftImage *gtBorders, float toleranceDist); 

int iftRandomSelectionWeightedByOrder(int n);

/*
 * @brief SLIC-like post-processing.
 *
 *   Turns each connected component into an unique label and relabel components below
 *     \c minSize (in pixels) to a neighbor's label.
 */
iftImage *iftForceLabelMapConnectivity(iftImage *labelMap, int minSize);

/*
 * @brief Reads CSV including header, if it exists.
 *
 *   Direct edit of original iftReadCSV code.
 */
iftCSV *iftReadCSVWithHeader(const char *csv_pathname, const char separator, bool *has_header);

#ifdef __cplusplus
}
#endif

#endif // _IFT_EXPERIMENT_UTILITY_H_ 
