#include "Image.h"

//=============================================================================
// Constructors & Deconstructors
//=============================================================================
Image *createImage(int num_rows, int num_cols, int num_channels)
{
    Image *new_img;

    new_img = (Image *)calloc(1, sizeof(Image));

    new_img->num_rows = num_rows;
    new_img->num_cols = num_cols;
    new_img->num_pixels = num_rows * num_cols;
    new_img->num_channels = num_channels;

    new_img->val = (int **)calloc(new_img->num_pixels, sizeof(int *));
#pragma omp parallel for
    for (int i = 0; i < new_img->num_pixels; i++)
        new_img->val[i] = (int *)calloc(num_channels, sizeof(int));

    return new_img;
}

void freeImage(Image **img)
{
    if (*img != NULL)
    {
        Image *tmp;

        tmp = *img;

        for (int i = 0; i < tmp->num_pixels; i++)
            free(tmp->val[i]);
        free(tmp->val);

        free(tmp);

        *img = NULL;
    }
}

NodeAdj *create4NeighAdj()
{
    NodeAdj *adj_rel;

    adj_rel = (NodeAdj *)calloc(1, sizeof(NodeAdj));

    adj_rel->size = 4;
    adj_rel->dx = (int *)calloc(4, sizeof(int));
    adj_rel->dy = (int *)calloc(4, sizeof(int));

    adj_rel->dx[0] = -1;
    adj_rel->dy[0] = 0; // Left
    adj_rel->dx[1] = 1;
    adj_rel->dy[1] = 0; // Right

    adj_rel->dx[2] = 0;
    adj_rel->dy[2] = -1; // Top
    adj_rel->dx[3] = 0;
    adj_rel->dy[3] = 1; // Bottom

    return adj_rel;
}

NodeAdj *create8NeighAdj()
{
    NodeAdj *adj_rel;

    adj_rel = (NodeAdj *)calloc(1, sizeof(NodeAdj));
    adj_rel->size = 8;
    adj_rel->dx = (int *)calloc(8, sizeof(int));
    adj_rel->dy = (int *)calloc(8, sizeof(int));

    adj_rel->dx[0] = -1;
    adj_rel->dy[0] = 0; // Center-Left
    adj_rel->dx[1] = 1;
    adj_rel->dy[1] = 0; // Center-Right

    adj_rel->dx[2] = 0;
    adj_rel->dy[2] = -1; // Top-Center
    adj_rel->dx[3] = 0;
    adj_rel->dy[3] = 1; // Bottom-Center

    adj_rel->dx[4] = -1;
    adj_rel->dy[4] = 1; // Bottom-Left
    adj_rel->dx[5] = 1;
    adj_rel->dy[5] = -1; // Top-Right

    adj_rel->dx[6] = -1;
    adj_rel->dy[6] = -1; // Top-Left
    adj_rel->dx[7] = 1;
    adj_rel->dy[7] = 1; // Bottom-Right
    return adj_rel;
}

void freeNodeAdj(NodeAdj **adj_rel)
{
    if (*adj_rel != NULL)
    {
        NodeAdj *tmp;

        tmp = *adj_rel;

        free(tmp->dx);
        free(tmp->dy);
        free(tmp);

        *adj_rel = NULL;
    }
}

//=============================================================================

inline bool areValidNodeCoordsImage(int num_rows, int num_cols, NodeCoords coords)
{
    return (coords.x >= 0 && coords.x < num_cols) &&
           (coords.y >= 0 && coords.y < num_rows);
}

inline int getNodeIndexImage(int num_cols, NodeCoords coords)
{
    return coords.y * num_cols + coords.x;
}

inline int getIndexImage(int num_cols, int row_index, int col_index)
{
    return row_index * num_cols + col_index;
}

inline NodeCoords getAdjacentNodeCoords(NodeAdj *adj_rel, NodeCoords coords, int id)
{
    NodeCoords adj_coords;

    adj_coords.x = coords.x + adj_rel->dx[id];
    adj_coords.y = coords.y + adj_rel->dy[id];

    return adj_coords;
}

inline NodeCoords getNodeCoordsImage(int num_cols, int index)
{
    NodeCoords coords;

    coords.x = index % num_cols;
    coords.y = index / num_cols;

    return coords;
}


