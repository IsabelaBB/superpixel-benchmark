#ifndef IMAGE_H
#define IMAGE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "Utils.h"
#include "Color.h"
#include <ctype.h>

//=============================================================================
// Macros
//=============================================================================
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define STR_DEFAULT_SIZE 4096
#define ImgVal(img, x, y, z) img->val[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]

//=============================================================================
// Structures
//=============================================================================
typedef struct
{
    int x, y;
} NodeCoords;

typedef struct
{
    int size;
    int *dx, *dy; // Coordinate shifts in each axis
} NodeAdj;

typedef struct
{
    int num_cols, num_rows, num_channels, num_pixels;
    int **val; // Access by val[i < num_pixels][f < num_channels]
} Image;

typedef struct _PGMData {
    int row;
    int col;
    int max_gray;
    int **matrix;
} PGMData;

//=============================================================================
// Constructors & Deconstructors
//=============================================================================
Image *createImage(int num_rows, int num_cols, int num_channels); // Zero-filled
void freeImage(Image **img);
NodeAdj *create4NeighAdj(); // 4-neighborhood
NodeAdj *create8NeighAdj(); // 8-neighborhood
void freeNodeAdj(NodeAdj **adj_rel);
//=============================================================================


int getMaximumValue(Image *img, int channel); // For all channels, set channel = -1
int getMinimumValue(Image *img, int channel); //
int getNormValue(Image *img); // For 8- and 16-bit, norm is 255 and 65535

bool areValidNodeCoordsImage(int num_rows, int num_cols, NodeCoords coords);
int getNodeIndexImage(int num_cols, NodeCoords coords);
int getIndexImage(int num_cols, int row_index, int col_index);
double euclDistance(float *feat1, float *feat2, int num_feats); // L2-norm
NodeCoords getAdjacentNodeCoords(NodeAdj *adj_rel, NodeCoords coords, int id);
NodeCoords getNodeCoordsImage(int num_cols, int index);

double euclDistanceCoords(NodeCoords feat1, NodeCoords feat2);
int* readPGM(const char *file_name);

#ifdef __cplusplus
}
#endif

#endif // IMAGE_H