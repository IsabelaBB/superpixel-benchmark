#ifndef EVAL_H
#define EVAL_H

#ifdef __cplusplus
extern "C" {
#endif


#include "PrioQueue.h"
#include "Image.h"
#include "Utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "ift.h"

void RGBtoYCbCr(int* cin, int num_channels, int normalization_value, int *cout);

int getNumSuperpixels(int *L, int num_pixels, int maxLabel);
double *getImageVariance_channels(Image *image);
int relabelSuperpixels(int *labels, int num_rows, int num_cols, int connectivity);
int enforceNumSuperpixel(int *labels, int num_rows, int num_cols, Image *image, int numDesiredSpx);

void RBD(Image *image, int *labels, int label, int nbuckets, int *alpha, double **Descriptor);
double *SIRS(int *labels, Image *image, int alpha, int nbuckets, char *reconFile, double gauss_variance, double *score);
double *computeExplainedVariation(int *labels, Image *image, char *reconFile, double *score);

bool is4ConnectedBoundaryPixel(int *labels, int num_rows, int num_cols, int i, int j);
bool is4ConnectedBoundaryPixel2(iftImage *gt, int i, int j);

void computeIntersectionMatrix(int *labels, iftImage *gt,
                               int **intersection_matrix, int *superpixel_sizes, int *gt_sizes, int cols, int rows);
double computeBoundaryRecall(int *labels, iftImage *gt, float d);
double computeUndersegmentationError(int *labels, iftImage *gt);
double computeCompactness(int *labels, int num_rows, int num_cols);

#ifdef __cplusplus
}
#endif

#endif // EVAL_H