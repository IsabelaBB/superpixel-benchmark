/*****************************************************************************\
* iftMetrics.h
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-03-08
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#ifndef IFT_METRICS_H
#define IFT_METRICS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift.h"

/*****************************************************************************\
*
*                               PUBLIC FUNCTIONS
*
\*****************************************************************************/
/*
	Calculates the boundary recall as presented by [1]. This metric is applicable
	 to any image (2D, 3D or video)

	[1] - D. Stutz, A. Hermans and B. Leibe. "Superpixels: An evaluation of the 
	state-of-the-art" Computer Vision and Image Understanding, 2018, pp. 1-27 
*/
float iftBoundRecall
(const iftImage *label_img, const iftImage *gt_img);

/*
	Calculates the under-segmentation error as proposed by [1]. This metric is 
	applicable to any image (2D, 3D or video)

	[1] - P. Neubert, P. Protzel, "Superpixel benchmark and comparison", Forum 
	Bildverarbeitung, 2012
*/
float iftUnderSegmError
(const iftImage *label_img, const iftImage *gt_img);

#ifdef __cplusplus
}
#endif

#endif