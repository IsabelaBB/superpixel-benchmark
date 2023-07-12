/**
 * @file
 *
 * @details General utility to extract features from superpixel
 *   segmentations. All methods work for both 2D and 3D unless
 *   explicitly noticed. To avoid using the non-standard terminology
 *   'supervoxel', we only refer to 'superpixels' and 'pixels' even
 *   in methods that work in 3D.
 * @details We assume superpixel segmentations are represented by
 *  a \c iftImage containing the superpixel label map.
 * @details For now results are always stored into an \c iftMatrix
 *   where the row r corresponds to the superpixel of label r.
 */

#ifndef _IFT_SUPERPIXEL_FEATURES_H_
#define _IFT_SUPERPIXEL_FEATURES_H_

#include "ift.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Use medoid of some color space as superpixel feature.
 */
iftMatrix * iftComputeSuperpixelFeaturesByColorSpaceMedoid(const iftImage *superpixelLabelMap, const iftImage *img, iftColorSpace colorSpace);  

/**
 * @brief Use medoid of pixel values as superpixel feature.
 *
 * @details For each superpixel, calculate the mean (average) value of
 *   its pixels' feature vectors (i.e. \c mimg values) and then use
 *   the pixel feature vector closest to the mean as the superpixel
 *   feature vector.
 * @details The returned \c iftDataSet is owned by the caller.
 * @see iftDestroyDataSet()
 *
 * @param[in] mimg Base \c iftMImage.
 * @param[in] superpixelLabelMap Some superpixel segmentation of mimg.
 * @return Data matrix where each row corresponds to the corresponding
 * @pre \c mimg and \c superpixelLabelMap are not NULL and share the
 *      same spatial domain.
 * @pre \c superpixelLabelMap does not contain label gaps.
 */
iftMatrix * iftComputeSuperpixelFeaturesByMImageMedoid(const iftMImage *mimg, const iftImage *superpixelLabelMap);

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
 *
 * @param[in] mimg Base \c iftMImage.
 * @param[in] superpixelLabelMap Some superpixel segmentation of mimg.
 * @return Data matrix where each row corresponds to the corresponding
 * @pre \c mimg and \c superpixelLabelMap are not NULL and share the
 *      same spatial domain.
 * @pre \c superpixelLabelMap does not contain label gaps.
 */
iftMatrix * iftComputeSuperpixelFeaturesByMImageMean(const iftMImage *mimg, const iftImage *superpixelLabelMap);

/**
 * @brief Computes color histogram for superpixels.
 * @author Felipe Lemes Galvao
 * @date January 5, 2017
 *
 * @details There is a very similar function that creates a iftDataSet
 *   of superpixel color histogram features. This one exists more as a
 *   supplement for the superpixel contextual features.
 *
 * @return Each row is a superpixel histogram (NOT normalized).
 */
iftMatrix * iftComputeSuperpixelFeaturesByColorHistogram(const iftImage *img, const iftImage *labelMap, int binsPerColor);

#ifdef __cplusplus
}
#endif

#endif // _IFT_SUPERPIXEL_FEATURES_H_
