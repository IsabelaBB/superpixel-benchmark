
#include "PrioQueue.h"
#include "Image.h"
#include "Utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include "ift.h"
#include <iostream>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

using namespace cv;
using namespace std;

// #define DEBUG

typedef struct TextsConfig
{
    Vec3b textColor;
    char text[10];
} TextsConfig;

typedef struct Args
{
    char *img_path, *label_path, *label_ext;
    char *imgScoresPath, *imgRecon;
    char *logFile, *dLogFile;
    int buckets, alpha, metric, k;
    bool drawScores;
    double gauss_variance;
    char *saveLabels;
} Args;

void usage()
{
    printf("Usage: main args [optional args] \n");
    printf("----------------------------------\n");
    printf("args:\n");
    printf("--img     - Original image or gt file/directory. Metrics 1,2, and 7: original image. Metrics 3 and 4: image ground-truth \n");
    printf("--eval  - Superpixel metric. Type: int. {1:SIRS, 2:EV, 3:BR, 4:UE, 5:CO, 6:Enforce connectivity, 7:Enforce superpixels' number} \n");
    printf("--label  - A pgm/png file with labeled image \n");
    printf("--ext     - File extension for image labels. Default: pgm \n");

    printf("--log     - txt file with overall results (only used when --img is a directory)\n");
    printf("--dlog    - txt file with evaluated results (for each image) \n");

    printf("--buckets   - Used in metric 1. Default: 16. Type: int \n");
    printf("--alpha     - Used in metric 1. Default: 4. Type: int \n");
    printf("--gaussVar  - Used in metric 1. Default: 0.01. Type: double \n");
    printf("--k         - Used in metric 7. Desired number of superpixels. Type: int \n");

    printf("--imgScores  -  Used in metrics 1 and 2. Save labeled images whose color map indicates SIRS/EV scores. \n");
    printf("--drawScores -  Used in metrics 1 and 2. Boolean option when using --imgScores to show score values \n");
    printf("--recon      -  Used in metrics 1 and 2. Save the reconstructed image \n");
    printf("--save       -  Used in metrics 6 and 7. Save image superpixels after enforce coonectivity/minimum number of superpixels \n");

    printError("main", "Too many/few parameters");
}

bool initArgs(Args *args, int argc, char *argv[])
{
    char *nbucketsChar, *alphaChar, *metricChar, *drawScoresChar, *gauss_varianceChar, *kChar;

    args->img_path = parseArgs(argv, argc, "--img");
    args->label_path = parseArgs(argv, argc, "--label");
    args->label_ext = parseArgs(argv, argc, "--ext");
    args->logFile = parseArgs(argv, argc, "--log");
    args->dLogFile = parseArgs(argv, argc, "--dlog");

    args->imgScoresPath = parseArgs(argv, argc, "--imgScores");
    args->imgRecon = parseArgs(argv, argc, "--recon");
    args->saveLabels = parseArgs(argv, argc, "--save");
    kChar = parseArgs(argv, argc, "--k");
    nbucketsChar = parseArgs(argv, argc, "--buckets");
    alphaChar = parseArgs(argv, argc, "--alpha");
    metricChar = parseArgs(argv, argc, "--metric");
    drawScoresChar = parseArgs(argv, argc, "--drawScores");
    gauss_varianceChar = parseArgs(argv, argc, "--gaussVar");

    if (strcmp(args->saveLabels, "-") == 0)
        args->saveLabels = NULL;

    args->k = strcmp(kChar, "-") != 0 ? atoi(kChar) : 0;
    args->gauss_variance = strcmp(gauss_varianceChar, "-") != 0 ? atof(gauss_varianceChar) : 0.01;
    args->buckets = strcmp(nbucketsChar, "-") != 0 ? atoi(nbucketsChar) : 16;
    args->alpha = strcmp(alphaChar, "-") != 0 ? atoi(alphaChar) : 4;
    args->metric = strcmp(metricChar, "-") != 0 ? atoi(metricChar) : 1;
    args->drawScores = strcmp(drawScoresChar, "-") != 0 ? atoi(drawScoresChar) : false;
    args->gauss_variance = strcmp(gauss_varianceChar, "-") != 0 ? atof(gauss_varianceChar) : 0.01;

    if (strcmp(args->logFile, "-") == 0)
        args->logFile = NULL;
    if (strcmp(args->dLogFile, "-") == 0)
        args->dLogFile = NULL;
    if (strcmp(args->imgScoresPath, "-") == 0)
        args->imgScoresPath = NULL;
    if (strcmp(args->imgRecon, "-") == 0)
        args->imgRecon = NULL;
    if (strcmp(args->label_ext, "-") == 0)
        sprintf(args->label_ext, "pgm");

    if (strcmp(args->img_path, "-") == 0 || strcmp(args->label_path, "-") == 0 || strcmp(args->label_ext, "-") == 0 ||
        (args->metric == 1 && (args->buckets < 1 || args->alpha < 1 || args->alpha > args->buckets * 7)) ||
        args->metric < 1 || args->metric > 7 || (args->k < 1 && args->metric == 7))
        return false;

    return true;
}

Image *loadImage(const char *filepath)
{
    int num_channels, num_rows, num_cols;
    unsigned char *data;
    Image *new_img;

    data = stbi_load(filepath, &num_cols, &num_rows, &num_channels, 0);

    if (data == NULL)
        printError("loadImage", "Could not load the image <%s>", filepath);

    new_img = createImage(num_rows, num_cols, num_channels);

#pragma omp parallel for
    for (int i = 0; i < new_img->num_pixels; i++)
    {
        new_img->val[i] = (int *)calloc(new_img->num_channels, sizeof(int));

        for (int j = 0; j < new_img->num_channels; j++)
            new_img->val[i][j] = data[i * new_img->num_channels + j];
    }

    stbi_image_free(data);

    return new_img;
}

void writeImagePGMInt(int *img, const char *filepath, int num_cols, int num_rows)
{
    int max_val, min_val;
    FILE *fp;

    fp = fopen(filepath, "wb");

    if (fp == NULL)
        printError("writeImagePGM", "Could not open the file <%s>", filepath);

    min_val = max_val = img[0];
    for (int i = 0; i < num_cols * num_rows; i++)
    {
        if (min_val > img[i])
            min_val = img[i];
        if (max_val < img[i])
            max_val = img[i];
    }

    fprintf(fp, "P5\n");
    fprintf(fp, "%d %d\n", num_cols, num_rows);
    fprintf(fp, "%d\n", max_val);

    // 8-bit PGM file
    if (max_val < 256 && min_val >= 0)
    {
        unsigned char *data;

        data = (unsigned char *)calloc(num_cols * num_rows, sizeof(unsigned char));

        for (int i = 0; i < num_cols * num_rows; i++)
            data[i] = (unsigned char)img[i];

        fwrite(data, sizeof(unsigned char), num_cols * num_rows, fp);

        free(data);
    }
    // 16-bit PGM file
    else if (max_val < 65536 && min_val >= 0)
    {
        unsigned short *data;

        data = (unsigned short *)calloc(num_cols * num_rows, sizeof(unsigned short));

        for (int i = 0; i < num_cols * num_rows; i++)
            data[i] = (unsigned short)img[i];

        for (int i = 0; i < num_cols * num_rows; i++)
        {
            int high, low;

            high = ((data[i]) & 0x0000FF00) >> 8;
            low = (data[i]) & 0x000000FF;

            fputc(high, fp);
            fputc(low, fp);
        }

        free(data);
    }
    else
        printError("writeImagePGM", "Invalid min/max spel values <%d,%d>", min_val, max_val);

    fclose(fp);
}

int filterDir(const struct dirent *name)
{
    int pos = 0;
    while (name->d_name[pos] != '.')
    {
        if (name->d_name[pos] == '\0')
            return 0;
        pos++;
    }
    pos++;
    if (name->d_name[pos] == '\0')
        return 0;

    int extSize = 0;
    while (name->d_name[pos + extSize] != '\0')
    {
        extSize++;
    }

    char ext[extSize];
    int pos2 = 0;
    while (pos2 < extSize)
    {
        ext[pos2] = name->d_name[pos + pos2];
        pos2++;
    }

    if ((extSize == 3 && ((ext[0] == 'p' && ext[1] == 'n' && ext[2] == 'g') ||
                          (ext[0] == 'j' && ext[1] == 'p' && ext[2] == 'g') ||
                          (ext[0] == 'p' && ext[1] == 'p' && ext[2] == 'm') ||
                          (ext[0] == 'p' && ext[1] == 'g' && ext[2] == 'm') ||
                          (ext[0] == 'b' && ext[1] == 'm' && ext[2] == 'p'))) ||
        (extSize == 4 && ext[0] == 't' && ext[1] == 'i' && ext[2] == 'f' && ext[2] == 'f'))
        return 1;

    return 0;
}

int getNumSuperpixels(int *L, int num_pixels, int maxLabel)
{
    int tmp[max(num_pixels, maxLabel)], numLabels = 0;

    for (int i = 0; i < max(num_pixels, maxLabel); i++)
        tmp[i] = 0;
    for (int i = 0; i < num_pixels; i++)
        tmp[L[i]] = -1;
    for (int i = 0; i < max(num_pixels, maxLabel); i++)
    {
        if (tmp[i] == -1)
            numLabels++;
    }

#ifdef DEBUG
    printf("Number of superpixels: %d\n", numLabels);
#endif

    return numLabels;
}

double *getImageVariance_channels(Image *image)
{
    double acum_color[image->num_channels];
    double *variance;

    variance = (double *)calloc(image->num_channels, sizeof(double));

    for (int c = 0; c < image->num_channels; c++)
        acum_color[c] = 0.0;

    for (int i = 0; i < image->num_pixels; i++)
    {
        for (int c = 0; c < image->num_channels; c++)
            acum_color[c] += (double)image->val[i][c] / 255.0;
    }

    for (int c = 0; c < image->num_channels; c++)
        acum_color[c] /= (double)image->num_pixels;

    for (int i = 0; i < image->num_pixels; i++)
    {
        for (int c = 0; c < image->num_channels; c++)
            variance[c] += ((double)image->val[i][c] / 255.0 - acum_color[c]) * ((double)image->val[i][c] / 255.0 - acum_color[c]);
    }

    for (int c = 0; c < image->num_channels; c++)
        variance[c] /= (double)image->num_pixels;
    return variance;
}

int relabelSuperpixels(int *labels, int num_rows, int num_cols, int connectivity)
{
    PrioQueue *queue;
    NodeAdj *adj_rel;
    double *visited;
    int *mapLabel, maxLabel = 0;

    visited = (double *)calloc(num_rows * num_cols, sizeof(double));

    for (int i = 0; i < num_rows * num_cols; i++)
    {
        visited[i] = num_rows * num_cols;
        if (labels[i] > maxLabel)
            maxLabel = labels[i];
    }

    mapLabel = (int *)calloc(maxLabel + 1, sizeof(int));

    for (int i = 0; i < maxLabel; i++)
        mapLabel[i] = i;

    queue = createPrioQueue(num_rows * num_cols, visited, MINVAL_POLICY);

    adj_rel = create8NeighAdj();

    int label = labels[0];
    int newLabel = 0;
    visited[0] = 0;
    insertPrioQueue(&queue, 0);
    mapLabel[label] = -1;

    while (!isPrioQueueEmpty(queue))
    {
        int vertex = popPrioQueue(&queue);
        int vertexLabel = labels[vertex];
        NodeCoords vertexCoords;
        vertexCoords = getNodeCoordsImage(num_cols, vertex);

        if (visited[vertex] > 0)
        {
            newLabel++;
            mapLabel[vertexLabel] = -1;
            label = vertexLabel;
        }

        labels[vertex] = newLabel;

        for (int i = 0; i < adj_rel->size; i++)
        {
            NodeCoords adjCoords;

            adjCoords = getAdjacentNodeCoords(adj_rel, vertexCoords, i);

            if (areValidNodeCoordsImage(num_rows, num_cols, adjCoords))
            {
                int adjVertex = getNodeIndexImage(num_cols, adjCoords);

                if (queue->state[adjVertex] != BLACK_STATE)
                {
                    if (labels[adjVertex] == label)
                        visited[adjVertex] = 0;

                    if (queue->state[adjVertex] == WHITE_STATE)
                        insertPrioQueue(&queue, adjVertex);
                    else
                    {
                        if (labels[adjVertex] == label)
                            moveIndexUpPrioQueue(&queue, adjVertex);
                    }
                }
            }
        }
    }

    free(visited);
    free(mapLabel);
    freePrioQueue(&queue);
    freeNodeAdj(&adj_rel);

    newLabel++;

    return getNumSuperpixels(labels, num_cols * num_rows, newLabel);
}

int enforceNumSuperpixel(int *labels, int num_rows, int num_cols, Image *image, int numDesiredSpx)
{
    NodeAdj *adj_rel;

    adj_rel = create8NeighAdj();
    int numSpx = relabelSuperpixels(labels, num_rows, num_cols, 8); // ensure connectivity

    float meanColor[numSpx][3];
    std::vector<std::vector<int>> neighbors(numSpx);
    std::vector<int> ids(numSpx);
    std::vector<int> sizeSpx(numSpx);
    int new_labels[numSpx];

    for (int i = 0; i < numSpx; i++)
    {
        sizeSpx[i] = 0;
        meanColor[i][0] = meanColor[i][1] = meanColor[i][2] = 0;
        neighbors[i] = std::vector<int>(0);
        ids[i] = i;
        new_labels[i] = i;
    }

    for (int node = 0; node < num_rows * num_cols; node++)
    {
        NodeCoords coords;
        int label = labels[node];

        coords = getNodeCoordsImage(num_cols, node);

        sizeSpx[label] += 1;
        meanColor[label][0] += (float)image->val[node][0];
        meanColor[label][1] += (float)image->val[node][1];
        meanColor[label][2] += (float)image->val[node][2];

        for (int j = 0; j < adj_rel->size; j++)
        {
            NodeCoords adjCoords = getAdjacentNodeCoords(adj_rel, coords, j);
            if (areValidNodeCoordsImage(num_rows, num_cols, adjCoords))
            {
                int adjLabel = labels[getNodeIndexImage(num_cols, adjCoords)];
                if (adjLabel != label)
                {
                    if (std::find(neighbors[label].begin(), neighbors[label].end(), adjLabel) == std::end(neighbors[label]))
                    {
                        neighbors[label].push_back(adjLabel);
                    }
                }
            }
        }
    }

    freeNodeAdj(&adj_rel);

    std::sort(ids.begin(), ids.end(), [&sizeSpx](int i, int j)
              { return sizeSpx[i] < sizeSpx[j]; });

    int k = 0;
    while (k <= (int)ids.size() && (int)ids.size() > numDesiredSpx)
    {

        int label = ids[k];
        float minDistance = INFINITY;
        float localMean[3];
        int new_label = -1;
        int parent = new_labels[label];

        while (label != parent)
        {
            label = parent;
            parent = new_labels[label];
        }

        localMean[0] = meanColor[label][0] / (float)sizeSpx[label];
        localMean[1] = meanColor[label][1] / (float)sizeSpx[label];
        localMean[2] = meanColor[label][2] / (float)sizeSpx[label];

        for (int adj = 0; adj < (int)neighbors[label].size(); adj++)
        {

            float adjMean[3];
            int adjLabel = neighbors[label][adj];
            parent = new_labels[adjLabel];

            while (adjLabel != parent)
            {
                adjLabel = parent;
                parent = new_labels[adjLabel];
            }

            if (adjLabel != label)
            {
                adjMean[0] = meanColor[adjLabel][0] / (float)sizeSpx[adjLabel];
                adjMean[1] = meanColor[adjLabel][1] / (float)sizeSpx[adjLabel];
                adjMean[2] = meanColor[adjLabel][2] / (float)sizeSpx[adjLabel];

                float distance = (localMean[0] - adjMean[0]) * (localMean[0] - adjMean[0]) + (localMean[1] - adjMean[1]) * (localMean[1] - adjMean[1]) + (localMean[2] - adjMean[2]) * (localMean[2] - adjMean[2]);

                if (distance < minDistance)
                {
                    minDistance = distance;
                    new_label = adjLabel;
                }
            }
            else
            {
                neighbors[label].erase(neighbors[label].begin() + adj);
            }
        }

        if (new_label != -1)
        {
            new_labels[label] = new_label;

            meanColor[new_label][0] += meanColor[label][0];
            meanColor[new_label][1] += meanColor[label][1];
            meanColor[new_label][2] += meanColor[label][2];
            sizeSpx[new_label] += sizeSpx[label];

            remove(neighbors[new_label].begin(), neighbors[new_label].end(), label);

            for (int adj = 0; adj < (int)neighbors[label].size(); adj++)
            {
                int adjLabel = neighbors[label][adj];

                if (adjLabel != label && adjLabel != new_label && std::find(neighbors[new_label].begin(), neighbors[new_label].end(), adjLabel) == std::end(neighbors[new_label]))
                {
                    neighbors[new_label].push_back(adjLabel);
                }
            }

            ids.erase(ids.begin() + k); // remove merged superpixel
            std::sort(ids.begin(), ids.end(), [&sizeSpx](int i, int j)
                      { return sizeSpx[i] < sizeSpx[j]; }); // sort again because adjLabel grown up
        }
        else
        {
            k++;
        }
    }

    for (k = 0; k < numSpx; k++)
    {
        meanColor[k][0] = -1;
    }

    // atualiza os rótulos que foram unidos
    for (k = 0; k < numSpx; k++)
    {
        int label;
        int parent;
        label = k;
        parent = new_labels[label];

        while (label != parent)
        {
            label = parent;
            parent = new_labels[label];
            new_labels[label] = parent;
        }
        label = k;
        while (label != parent)
        {
            int tmp = new_labels[label];
            new_labels[label] = parent;
            label = tmp;
        }
        meanColor[parent][0] = (float)parent;
    }

    // rotula os pixels com os rótulos novos
    for (int i = 0; i < num_rows * num_cols; i++)
    {
        labels[i] = new_labels[labels[i]];
    }

    return relabelSuperpixels(labels, num_rows, num_cols, 8);
}

//==========================================================

void RBD(Image *image, int *labels, int label, int nbuckets, int *alpha, float **Descriptor)
{
    /* Compute the superpixels descriptors
        image : RGB image
        labels     : Image labels (0,K-1)
        Descriptor : Descriptor[num_channels][alpha]
    */

    PrioQueue *queue;
    int num_histograms = pow(2, image->num_channels) - 1;
    long int ColorHistogram[num_histograms][nbuckets][image->num_channels]; // Descriptor[image->num_channels][nbuckets]
    double V[num_histograms * nbuckets];                                    // buckets priority : V[image->num_channels][nbuckets]
    int superpixel_size;

    superpixel_size = 0;

    queue = createPrioQueue(num_histograms * nbuckets, V, MINVAL_POLICY);

    for (int h = 0; h < num_histograms; h++)
    {
        for (int b = 0; b < nbuckets; b++)
        {
            V[h * nbuckets + b] = 0;
            for (int c = 0; c < image->num_channels; c++)
                ColorHistogram[h][b][c] = 0;
        }
    }

#ifdef DEBUG
    printf("RBD: Compute histogram\n");
#endif

    // compute histograms
    for (int i = 0; i < image->num_pixels; i++)
    {
        if (labels[i] == label)
        {
            superpixel_size++;
            int hist_id = 0;
            int bit = 1;

            for (int c1 = 0; c1 < image->num_channels; c1++)
            {
                for (int c2 = 0; c2 < image->num_channels; c2++)
                {
                    if (image->val[i][c1] < image->val[i][c2])
                    {
                        hist_id = hist_id | bit;
                        c2 = image->num_channels;
                    }
                }
                bit = bit << 1;
            }
            hist_id = ~hist_id;
            hist_id *= -1;

            // find one max channel index
            int max_channel = 0;
            int tmp_hist_id = hist_id;
            while (tmp_hist_id % 2 == 0)
            {
                tmp_hist_id = tmp_hist_id >> 1;
                max_channel++;
            }

            int bin = floor(((float)image->val[i][max_channel] / 255.0) * nbuckets);
            hist_id--;

            if (hist_id < 0 || hist_id > 6)
                printf(">> hist_id:%d \n", hist_id);
            if (bin < 0 || bin > nbuckets)
                printf(">> bin:%d \n", bin);
            if (hist_id * nbuckets + bin < 0 || hist_id * nbuckets + bin > num_histograms * nbuckets)
                printf(">> hist_id * nbuckets + bin:%d \n", hist_id * nbuckets + bin);

            V[hist_id * nbuckets + bin]++;

            for (int c = 0; c < image->num_channels; c++)
                ColorHistogram[hist_id][bin][c] += (long int)image->val[i][c];
        }
    }

#ifdef DEBUG
    printf("getSuperpixelDescriptor: Find the higher alpha buckets\n");
#endif

    for (int c = 0; c < num_histograms; c++)
    {
        for (int b = 0; b < nbuckets; b++)
        {
            if (V[c * nbuckets + b] > 0)
            {
                if (isPrioQueueEmpty(queue) || queue->last_elem_pos < (*alpha) - 1)
                    insertPrioQueue(&queue, c * nbuckets + b); // push (color, frequency) into Q, sorted by V[i]
                else
                {
                    if (!isPrioQueueEmpty(queue) && V[queue->node[0]] < V[c * nbuckets + b])
                    {
                        popPrioQueue(&queue);
                        insertPrioQueue(&queue, c * nbuckets + b); // push (color, frequency) into Q, sorted by V[i]
                    }
                }
            }
        }
    }

    // Get the higher alpha buckets
    if (isPrioQueueEmpty(queue))
    {
        (*alpha) = 0;
        int a = 0;

        for (int c = 0; c < image->num_channels; c++)
            Descriptor[a][c] = 0;
    }
    else
    {
        if (queue->last_elem_pos < (*alpha) - 1)
            (*alpha) = queue->last_elem_pos + 1;

        for (int a = (*alpha) - 1; a >= 0; a--)
        {
            int val = popPrioQueue(&queue);
            int bin = val % nbuckets;
            int hist = (val - bin) / nbuckets;

            for (int c = 0; c < image->num_channels; c++)
            {
                Descriptor[a][c] = (float)ColorHistogram[hist][bin][c] / (float)V[val]; // get the mean color
            }
        }
    }

#ifdef DEBUG
    printf("\nDESCRIPTOR: \n");
    for (int a = 0; a < (*alpha); a++)
    {
        printf("\t alpha %d (", a);
        for (int c = 0; c < image->num_channels; c++)
        {
            printf("%f ", Descriptor[a][c]);
        }
        printf(") \n\n");
    }
#endif

#ifdef DEBUG
    printf("\nFREE STRUCTURES: getSuperpixelDescriptor\n");
#endif

    freePrioQueue(&queue);
}

//==========================================================

void drawtorect(cv::Mat &mat, cv::Rect target, int face, int thickness, cv::Scalar color, const std::string &str)
{
    Size rect = getTextSize(str, face, 1.0, thickness, 0);

    target.width -= 8;
    target.height -= 8;
    target.x += 4;
    target.y += 4;

    if (target.height < 10 || ((float)target.width / (float)target.height < 1.8 && target.width < 36))
        return;
    if (target.width > 60)
    {
        target.x += (target.width - 60) / 2;
        target.width = 60;
    }
    if (target.height > 13)
    {
        target.y += (target.height - 13) / 2;
        target.height = 13;
    }

    double scalex = (double)target.width / (double)rect.width;
    double scaley = (double)target.height / (double)rect.height;
    double scale = min(scalex, scaley);
    int marginx = scale == scalex ? 0 : (int)((double)target.width * (scalex - scale) / scalex * 0.5);
    int marginy = scale == scaley ? 0 : (int)((double)target.height * (scaley - scale) / scaley * 0.5);

    putText(mat, str, Point(target.x + marginx, target.y + target.height - marginy), face, scale, color, thickness, 8, false);
}

Rect findMinRect(const Mat1b &src)
{
    Mat1f W(src.rows, src.cols, float(0));
    Mat1f H(src.rows, src.cols, float(0));

    Rect maxRect(0, 0, 0, 0);
    float maxArea = 0.f;

    for (int r = 0; r < src.rows; ++r)
    {
        for (int c = 0; c < src.cols; ++c)
        {
            if (src(r, c) == 0)
            {
                H(r, c) = 1.f + ((r > 0) ? H(r - 1, c) : 0);
                W(r, c) = 1.f + ((c > 0) ? W(r, c - 1) : 0);
            }

            float minw = W(r, c);
            for (int h = 0; h < H(r, c); ++h)
            {
                minw = min(minw, W(r - h, c));
                float area = (h + 1) * minw;
                if (area > maxArea)
                {
                    maxArea = area;
                    maxRect = Rect(Point(c - minw + 1, r - h), Point(c + 1, r + 1));
                }
            }
        }
    }

    return maxRect;
}

void createImageMetric(int *L, double *colorVariance, int num_rows, int num_cols, const char *filename, bool showScores)
{
    /*
    L : label map
    colorVariance : segmentaton error
    */

    NodeAdj *AdjRel;
    Mat image;
    int K;

    K = 0;
    for (int p = 0; p < num_rows * num_cols; p++)
        K = MAX(K, L[p]);
    K++;

    int superpixel_size[K];
    TextsConfig textsConfig[K];
    vector<vector<Point>> contours(K);
    int count_contours[K];

#ifdef DEBUG
    printf("createImageMetric: Alloc structures\n");
#endif

    AdjRel = create8NeighAdj();
    image = Mat::zeros(num_rows, num_cols, CV_8UC3);

    int max_Label = 0;
    for (int s = 0; s < K; s++)
    {
        superpixel_size[s] = 0;
        count_contours[s] = 0;
    }
    for (int p = 0; p < num_rows * num_cols; p++)
    {
        superpixel_size[L[p]]++;
        max_Label = MAX(max_Label, L[p]);
    }

    for (int s = 0; s < K; s++)
        contours[s] = vector<Point>(superpixel_size[s]);

#ifdef DEBUG
    printf("Iterate over the image pixels.. \n");
#endif
    for (int p = 0; p < num_rows * num_cols; p++)
    {
        Vec3b color, border;
        NodeCoords coords;
        int label = L[p];

        color[0] = color[1] = color[2] = (int)(255 * MIN(1, colorVariance[label])); // get superpixel color according to its error
        border[0] = border[1] = border[2] = (color[0] < 128) ? 255 : 0;

        bool isBorder = false;
        coords = getNodeCoordsImage(num_cols, p);

        for (int j = 0; j < AdjRel->size; j++)
        {
            NodeCoords adj_coords;
            adj_coords = getAdjacentNodeCoords(AdjRel, coords, j);

            if (areValidNodeCoordsImage(num_rows, num_cols, adj_coords))
            {
                int adj_index;
                adj_index = getNodeIndexImage(num_cols, adj_coords);
                if (label != L[adj_index])
                {
                    isBorder = true;
                    if (image.at<Vec3b>(Point(adj_coords.x, adj_coords.y))[0] == 0 && image.at<Vec3b>(Point(adj_coords.x, adj_coords.y))[0] == 255)
                        image.at<Vec3b>(Point(coords.x, coords.y)) = border;
                    else
                        image.at<Vec3b>(Point(coords.x, coords.y)) = image.at<Vec3b>(Point(adj_coords.x, adj_coords.y));
                    break;
                }
            }
        }
        if (!isBorder)
            image.at<Vec3b>(Point(coords.x, coords.y)) = color;

        contours[label][count_contours[label]] = Point(coords.x, coords.y);
        count_contours[label]++;

        gcvt(colorVariance[label], 2, textsConfig[label].text);
        textsConfig[label].textColor[0] = color[0] < 128 ? 255 : 0;
        textsConfig[label].textColor[1] = textsConfig[label].textColor[2] = textsConfig[label].textColor[0];
    }

#ifdef DEBUG
    printf("Apply colormap.. \n");
#endif
    // Apply the colormap
    applyColorMap(image, image, COLORMAP_OCEAN);

    if (showScores)
    {
        vector<vector<Point>> contours_poly(contours.size());
        vector<Rect> boundRect(contours.size());

#ifdef DEBUG
        printf("Draw scores.. \n");
#endif

        // draw scaled text according to region size: https://stackoverflow.com/questions/50353884/calculate-text-size
        // multidimentional vectors in c++: https://www.geeksforgeeks.org/2d-vector-in-cpp-with-user-defined-size/
        // fin min enclosed rectangle: https://stackoverflow.com/questions/34896431/creating-rectangle-within-a-blob-using-opencv

        for (size_t i = 0; i < contours.size(); ++i)
        {
            if (contours[i].size() > 0)
            {
                // Create a mask for each single blob
                Mat1b maskSingleContour(num_rows, num_cols, uchar(0));

                for (size_t j = 0; j < contours[i].size(); j++)
                {
                    maskSingleContour.at<uchar>(contours[i][j]) = 255;
                }

                // Find minimum rect for each blob
                Rect box = findMinRect(~maskSingleContour);

                // Draw rect
                drawtorect(image, box, cv::FONT_HERSHEY_PLAIN, 1, textsConfig[i].textColor, textsConfig[i].text);
            }
        }
    }

#ifdef DEBUG
    printf("Write final image.. \n");
#endif

    imwrite(filename, image);

#ifdef DEBUG
    printf("Free structure.. \n");
#endif

    freeNodeAdj(&AdjRel);
}

//==========================================================
// COLOR HOMOGENEITY MEASURES
//==========================================================

double *SIRS(int *labels, Image *image, int alpha, int nbuckets, char *reconFile, double gauss_variance, double *score)
{
    double *histogramVariation;
    float ***Descriptor;  // Descriptor[numSup][alpha][num_channels]
    int *descriptor_size; // descriptor_size[numSup]
    int emptySuperpixels;
    double **MSE;
    Mat recons;
    double **mean_buckets;         // mean_buckets[superpixels][num_channels];
    double **variation_descriptor; // variation_descriptor[superpixels][num_channels];

    (*score) = 0;
    emptySuperpixels = 0;

    recons = Mat::zeros(image->num_rows, image->num_cols, CV_8UC3);

    // get the higher label = number of superpixels
    int superpixels = 0;
    for (int i = 0; i < image->num_pixels; ++i)
    {
        if (labels[i] > superpixels)
            superpixels = labels[i];
    }
    superpixels++;

    histogramVariation = (double *)calloc(superpixels, sizeof(double));
    descriptor_size = (int *)calloc(superpixels, sizeof(int));
    Descriptor = (float ***)calloc(superpixels, sizeof(float **));
    MSE = (double **)calloc(superpixels, sizeof(double *));
    mean_buckets = (double **)calloc(superpixels, sizeof(double *));
    variation_descriptor = (double **)calloc(superpixels, sizeof(double *));

    std::vector<int> superpixelSize(superpixels, 0);
    for (int i = 0; i < image->num_pixels; ++i)
        superpixelSize[labels[i]]++;

    for (int s = 0; s < superpixels; s++)
    {
        if (superpixelSize[s] == 0)
            emptySuperpixels++;

#ifdef DEBUG
        printf("SUPERPIXEL %d: \n", s);
#endif
        descriptor_size[s] = alpha;
        Descriptor[s] = (float **)calloc(alpha, sizeof(float *));
        MSE[s] = (double *)calloc(image->num_channels, sizeof(double));
        mean_buckets[s] = (double *)calloc(image->num_channels, sizeof(double));
        variation_descriptor[s] = (double *)calloc(image->num_channels, sizeof(double));

        for (int c = 0; c < alpha; c++)
            Descriptor[s][c] = (float *)calloc(image->num_channels, sizeof(float));

#ifdef DEBUG
        printf("call RBD \n");
#endif
        RBD(image, labels, s, nbuckets, &(descriptor_size[s]), Descriptor[s]);

        for (int i = 0; i < image->num_channels; i++)
            MSE[s][i] = 0.0;

        for (int a = 0; a < descriptor_size[s]; a++)
        {
            for (int c = 0; c < image->num_channels; c++)
                mean_buckets[s][c] += ((double)Descriptor[s][a][c] / 255.0);
        }
        for (int c = 0; c < image->num_channels; c++)
            mean_buckets[s][c] /= descriptor_size[s];

        for (int a = 0; a < descriptor_size[s]; a++)
        {
            for (int c = 0; c < image->num_channels; c++)
                variation_descriptor[s][c] = MAX(variation_descriptor[s][c], abs(((double)Descriptor[s][a][c] / 255.0) - mean_buckets[s][c]));
        }
    }

#ifdef DEBUG
    printf("compute variation \n");
#endif

    for (int i = 0; i < image->num_pixels; i++)
    {
        int label;
        NodeCoords coords;
        double minVariance;
        int descIndex;

        label = labels[i];
        coords = getNodeCoordsImage(image->num_cols, i);
        Vec3b &recons_color = recons.at<Vec3b>(coords.y, coords.x);

        descIndex = 0;
        minVariance = 0;

        for (int c = 0; c < image->num_channels; c++)
            minVariance += ((double)image->val[i][c] / 255.0 - (double)Descriptor[label][0][c] / 255.0) * ((double)image->val[i][c] / 255.0 - (double)Descriptor[label][0][c] / 255.0);

        // find the most distance descriptor values
        for (int h1 = 1; h1 < descriptor_size[label]; h1++)
        {
            double val = 0;

            for (int c = 0; c < image->num_channels; c++)
                val += ((double)image->val[i][c] / 255.0 - (double)Descriptor[label][h1][c] / 255.0) * ((double)image->val[i][c] / 255.0 - (double)Descriptor[label][h1][c] / 255.0);

            if (val < minVariance)
            {
                minVariance = val;
                descIndex = h1;
            }
        }

        if (reconFile != NULL)
        {
            for (int c = 0; c < image->num_channels; c++)
                recons_color[image->num_channels - 1 - c] = Descriptor[label][descIndex][c];
        }
        for (int c = 0; c < image->num_channels; c++)
        {
            MSE[label][c] += pow(abs((double)image->val[i][c] / 255.0 - (double)Descriptor[label][descIndex][c] / 255.0), 2 - variation_descriptor[label][c]);
        }
    }

    double *variance = getImageVariance_channels(image);
    double sum_MSE[image->num_channels];

    for (int c = 0; c < image->num_channels; c++)
        sum_MSE[c] = 0.0;

#ifdef DEBUG
    printf("variance: ");
    for (int c = 0; c < image->num_channels; c++)
        printf("%f ", variance[c]);
    printf("\n");
    printf("\n\nScores:\n");
#endif
    for (int s = 0; s < superpixels; s++)
    {
        histogramVariation[s] = 0;
        for (int c = 0; c < image->num_channels; c++)
        {
            histogramVariation[s] += exp(-(MSE[s][c] / superpixelSize[s]) / gauss_variance);
            sum_MSE[c] += MSE[s][c];
        }
        histogramVariation[s] /= image->num_channels;

#ifdef DEBUG
        printf("%f + ", histogramVariation[s]);
#endif
    }

    (*score) = 0;
    for (int c = 0; c < image->num_channels; c++)
    {
        (*score) += sum_MSE[c];
    }
    (*score) /= image->num_channels;
    (*score) = exp(-((*score) / image->num_pixels) / (gauss_variance));

#ifdef DEBUG
    printf("= %f \n", (*score));
#endif

    for (int s = 0; s < superpixels; s++)
    {
        for (int a = 0; a < descriptor_size[s]; a++)
        {
            free(Descriptor[s][a]);
        }
        free(Descriptor[s]);
        free(MSE[s]);
        free(variation_descriptor[s]);
        free(mean_buckets[s]);
    }
    free(Descriptor);
    free(descriptor_size);
    free(MSE);
    free(variance);
    free(variation_descriptor);
    free(mean_buckets);

    if (reconFile != NULL)
        imwrite(reconFile, recons);

    return histogramVariation;
}

double *computeExplainedVariation(int *labels, Image *image, char *reconFile, double *score)
{
    (*score) = 0;
    
    double *supExplainedVariation;
    double *valuesTop, *valuesBottom;
    Mat recons;

    // get the higher label
    int superpixels = 0;
    for (int i = 0; i < image->num_pixels; ++i)
    {
        if (labels[i] > superpixels)
            superpixels = labels[i];
    }
    superpixels++;

    recons = Mat::zeros(image->num_rows, image->num_cols, CV_8UC3);

    supExplainedVariation = (double *)calloc(superpixels, sizeof(double));
    valuesTop = (double *)calloc(superpixels, sizeof(double));
    valuesBottom = (double *)calloc(superpixels, sizeof(double));

    std::vector<cv::Vec3f> mean(superpixels, cv::Vec3f(0, 0, 0));
    std::vector<cv::Vec3f> squared_mean(superpixels, cv::Vec3f(0, 0, 0));
    std::vector<int> count(superpixels, 0);

    cv::Vec3f overall_mean = 0;
    for (int i = 0; i < image->num_pixels; ++i)
    {
        for (int c = 0; c < image->num_channels; ++c)
        {
            mean[labels[i]][c] += image->val[i][c];
            overall_mean[c] += image->val[i][c];
        }
        count[labels[i]]++;
    }

    for (int i = 0; i < superpixels; ++i)
    {
        valuesTop[i] = 0;
        valuesBottom[i] = 0;
        for (int c = 0; c < image->num_channels; ++c)
            mean[i][c] /= count[i];
    }

    overall_mean /= image->num_rows * image->num_cols;

    for (int i = 0; i < image->num_pixels; ++i)
    {
        int label;
        NodeCoords coords;

        label = labels[i];
        coords = getNodeCoordsImage(image->num_cols, i);
        Vec3b &recons_color = recons.at<Vec3b>(coords.y, coords.x);

        for (int c = 0; c < image->num_channels; ++c)
        {
            /*sum_top += (mean[labels[i]][c] - overall_mean[c])
                    *(mean[labels[i]][c] - overall_mean[c]);
            sum_bottom += (image->val[i][c] - overall_mean[c])
                    *(image->val[i][c] - overall_mean[c]);*/
            valuesTop[label] += (mean[label][c] - overall_mean[c]) * (mean[labels[i]][c] - overall_mean[c]);
            valuesBottom[label] += (image->val[i][c] - overall_mean[c]) * (image->val[i][c] - overall_mean[c]);
            if (reconFile != NULL)
                recons_color[image->num_channels - 1 - c] = mean[labels[i]][c];
        }
    }

    double sum_top = 0;
    double sum_bottom = 0;

#ifdef DEBUG
    printf("\n\nScores:\n");
#endif
    for (int s = 0; s < superpixels; s++)
    {
        sum_top += valuesTop[s];
        sum_bottom += valuesBottom[s];

        if (valuesBottom[s] == 0)
            supExplainedVariation[s] = 1;
        else
            supExplainedVariation[s] = valuesTop[s] / valuesBottom[s];

#ifdef DEBUG
        printf("%f + ", supExplainedVariation[s]);
#endif
    }

    (*score) += sum_top / sum_bottom;

#ifdef DEBUG
    printf("= %f \n", (*score));
#endif

    free(valuesTop);
    free(valuesBottom);

    if (reconFile != NULL)
        imwrite(reconFile, recons);

    // return sum_top/sum_bottom;
    return supExplainedVariation;
}

//==========================================================
// DELINEATION MEASURES
//==========================================================

bool is4ConnectedBoundaryPixel(int *labels, int num_rows, int num_cols, int i, int j)
{

    int label_tmp = labels[getIndexImage(num_cols, i, j)];

    if (i > 0)
    {
        if (label_tmp != labels[getIndexImage(num_cols, i - 1, j)])
        {
            return true;
        }
    }

    if (i < num_rows - 1)
    {
        if (label_tmp != labels[getIndexImage(num_cols, i + 1, j)])
        {
            return true;
        }
    }

    if (j > 0)
    {
        if (label_tmp != labels[getIndexImage(num_cols, i, j - 1)])
        {
            return true;
        }
    }

    if (j < num_cols - 1)
    {
        if (label_tmp != labels[getIndexImage(num_cols, i, j + 1)])
        {
            return true;
        }
    }

    return false;
}

bool is4ConnectedBoundaryPixel2(iftImage *gt, int i, int j)
{
    int num_cols = gt->xsize;
    int label_tmp = gt->val[getIndexImage(num_cols, i, j)];

    if (i > 0)
    {
        if (label_tmp != gt->val[getIndexImage(num_cols, i - 1, j)])
        {
            return true;
        }
    }

    if (i < gt->ysize - 1)
    {
        if (label_tmp != gt->val[getIndexImage(num_cols, i + 1, j)])
        {
            return true;
        }
    }

    if (j > 0)
    {
        if (label_tmp != gt->val[getIndexImage(num_cols, i, j - 1)])
        {
            return true;
        }
    }

    if (j < num_cols - 1)
    {
        if (label_tmp != gt->val[getIndexImage(num_cols, i, j + 1)])
        {
            return true;
        }
    }

    return false;
}

void computeIntersectionMatrix(int *labels, iftImage *gt,
                               cv::Mat &intersection_matrix, std::vector<int> &superpixel_sizes, std::vector<int> &gt_sizes)
{
    int superpixels = 0;
    int gt_segments = 0;
    for (int i = 0; i < gt->n; ++i)
    {
        if (labels[i] > superpixels)
        {
            superpixels = labels[i];
        }
        if (gt->val[i] > gt_segments)
        {
            gt_segments = gt->val[i];
        }
    }

    superpixels++;
    gt_segments++;

    superpixel_sizes.resize(superpixels, 0);
    gt_sizes.resize(gt_segments, 0);

    intersection_matrix.create(gt_segments, superpixels, CV_32SC1);

    for (int i = 0; i < intersection_matrix.rows; ++i)
    {
        for (int j = 0; j < intersection_matrix.cols; ++j)
        {
            intersection_matrix.at<int>(i, j) = 0;
        }
    }

    for (int i = 0; i < gt->ysize; ++i)
    {
        for (int j = 0; j < gt->xsize; ++j)
        {
            int index = getIndexImage(gt->xsize, i, j);
            intersection_matrix.at<int>(gt->val[index], labels[index])++;
            superpixel_sizes[labels[index]]++;
            gt_sizes[gt->val[index]]++;
        }
    }
}

double computeBoundaryRecall(int *labels, iftImage *gt, float d)
{
    int H = gt->ysize; // num_rows
    int W = gt->xsize; // num_cols

    int r = std::round(d * std::sqrt(H * H + W * W));

    float tp = 0;
    float fn = 0;

    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            if (is4ConnectedBoundaryPixel2(gt, i, j))
            {

                bool pos = false;
                for (int k = std::max(0, i - r); k < std::min(H - 1, i + r) + 1; k++)
                {
                    for (int l = std::max(0, j - r); l < std::min(W - 1, j + r) + 1; l++)
                    {
                        if (is4ConnectedBoundaryPixel(labels, H, W, k, l))
                        {
                            pos = true;
                        }
                    }
                }

                if (pos)
                {
                    tp++;
                }
                else
                {
                    fn++;
                }
            }
        }
    }

    if (tp + fn > 0)
    {
        return tp / (tp + fn);
    }

    return 0;
}

double computeUndersegmentationError(int *labels, iftImage *gt)
{
    int N = gt->xsize * gt->ysize;

    cv::Mat intersection_matrix;
    std::vector<int> superpixel_sizes;
    std::vector<int> gt_sizes;

    computeIntersectionMatrix(labels, gt, intersection_matrix, superpixel_sizes, gt_sizes);

    float error = 0;
    for (int j = 0; j < intersection_matrix.cols; ++j)
    {

        int min = std::numeric_limits<int>::max();
        for (int i = 0; i < intersection_matrix.rows; ++i)
        {
            int superpixel_j_minus_gt_i = superpixel_sizes[j] - intersection_matrix.at<int>(i, j);

            if (superpixel_j_minus_gt_i < 0)
            {
                printf("Invalid intersection computed, set difference is negative! \n");
                return -1;
            }
            if (superpixel_j_minus_gt_i < min)
            {
                min = superpixel_j_minus_gt_i;
            }
        }

        error += min;
    }

    return error / N;
}

double computeCompactness(int *labels, int num_rows, int num_cols)
{

    int superpixels = 0;
    for (int i = 0; i < num_rows; ++i)
    {
        for (int j = 0; j < num_cols; ++j)
        {
            int label = labels[getIndexImage(num_cols, i, j)];
            if (label > superpixels)
            {
                superpixels = label;
            }
        }
    }
    superpixels++;

    std::vector<float> perimeter(superpixels, 0);
    std::vector<float> area(superpixels, 0); // = number of pixels

    for (int i = 0; i < num_rows; ++i)
    {
        for (int j = 0; j < num_cols; ++j)
        {

            int count = 0;
            int label = labels[getIndexImage(num_cols, i, j)];
            if (i > 0)
            {
                if (label != labels[getIndexImage(num_cols, i - 1, j)])
                {
                    count++;
                }
            }
            else
            {
                count++;
            }

            if (i < num_rows - 1)
            {
                if (label != labels[getIndexImage(num_cols, i + 1, j)])
                {
                    count++;
                }
            }
            else
            {
                count++;
            }

            if (j > 0)
            {
                if (label != labels[getIndexImage(num_cols, i, j - 1)])
                {
                    count++;
                }
            }
            else
            {
                count++;
            }

            if (j < num_cols - 1)
            {
                if (label != labels[getIndexImage(num_cols, i, j + 1)])
                {
                    count++;
                }
            }
            else
            {
                count++;
            }

            perimeter[label] += count;
            area[label] += 1;
        }
    }

    float compactness = 0;

    for (int i = 0; i < superpixels; ++i)
    {
        if (perimeter[i] > 0)
        {
            compactness += area[i] * (4 * M_PI * area[i]) / (perimeter[i] * perimeter[i]);
        }
    }

    compactness /= (num_rows * num_cols);

    return compactness;
}

double computeMetric(char *img_path, char *labels_path, char *saveImg, char *reconFile, Args args, int *numSuperpixels)
{
    iftImage *labels, *gt;
    Image *image;
    int *labeledImage;
    double *explainedVariation, score;
    int num_rows, num_cols;

    score = 0;
    labels = iftReadImageByExt(labels_path);
    num_cols = labels->xsize;
    num_rows = labels->ysize;

    labeledImage = (int *)malloc(labels->n * sizeof(int));
    int maxLabel = 0, minLabel = INT_MAX;
    for (int i = 0; i < labels->n; i++)
    {
        labeledImage[i] = labels->val[i];
        if (maxLabel < labeledImage[i])
            maxLabel = labeledImage[i];
        if (minLabel > labeledImage[i])
            minLabel = labeledImage[i];
    }

    // color homogeneity measures: original img and labels
    if (args.metric < 3 || (args.metric == 7 && saveImg != NULL))
    {
        image = loadImage(img_path);
    }

    // delineation measures: labels and ground-truth
    if (args.metric == 3 || args.metric == 4)
    {
        gt = iftReadImageByExt(img_path);
    }

    iftDestroyImage(&labels);

    maxLabel = relabelSuperpixels(labeledImage, num_rows, num_cols, 8);

    (*numSuperpixels) = getNumSuperpixels(labeledImage, num_cols * num_rows, maxLabel);

#ifdef DEBUG
    printf("Measure segm. error: ");
#endif

    switch (args.metric)
    {
    // ========================================
    // COLOR-BASED MEASURES
    // ========================================
    case 1:
#ifdef DEBUG
        printf("SIRS \n");
#endif
        explainedVariation = SIRS(labeledImage, image, args.alpha, args.buckets, reconFile, args.gauss_variance, &score);
        if (saveImg != NULL)
            createImageMetric(labeledImage, explainedVariation, image->num_rows, image->num_cols, saveImg, args.drawScores);
        break;

    case 2:
#ifdef DEBUG
        printf("EV \n");
#endif
        explainedVariation = computeExplainedVariation(labeledImage, image, reconFile, &score);
        if (saveImg != NULL)
            createImageMetric(labeledImage, explainedVariation, image->num_rows, image->num_cols, saveImg, args.drawScores);
        break;

    // ========================================
    // DELINEATION MEASURES
    // ========================================
    case 3:
#ifdef DEBUG
        printf("BR \n");
#endif
        score = computeBoundaryRecall(labeledImage, gt, 0.0025);
        break;

    case 4:
#ifdef DEBUG
        printf("UE \n");
#endif
        score = computeUndersegmentationError(labeledImage, gt);
        break;

    // ========================================
    // OTHER MEASURES
    // ========================================
    case 5:
#ifdef DEBUG
        printf("CO \n");
#endif
        score = computeCompactness(labeledImage, num_rows, num_cols);
        break;

    case 6:
#ifdef DEBUG
        printf("Connectivity \n");
#endif
        score = (double)relabelSuperpixels(labeledImage, num_rows, num_cols, 8);
        if (saveImg != NULL)
            writeImagePGMInt(labeledImage, saveImg, num_cols, num_rows);
        break;

    case 7:
#ifdef DEBUG
        printf("Enforce superpixels' number \n");
#endif
        if (saveImg != NULL)
        {
            score = enforceNumSuperpixel(labeledImage, num_rows, num_cols, image, (*numSuperpixels));
            if (score != (*numSuperpixels))
                printError("computeMetric", "Computing the wrong number of superpixels.");
            writeImagePGMInt(labeledImage, saveImg, num_cols, num_rows);
        }
        else
            score = 0;
        break;

    default:
        score = -1;
    }

#ifdef DEBUG
    printf("main: Free structures..\n");
#endif

    if (args.metric < 3 || (args.metric == 7 && saveImg != NULL))
    {
        freeImage(&image);
        if (args.metric < 3 && explainedVariation != NULL)
            free(explainedVariation);
    }

    if (args.metric == 3 || args.metric == 4)
    {
        iftDestroyImage(&gt);
    }

    free(labeledImage);
    return score;
}

void getImageName(char *fullImgPath, char *fileName)
{
    char *imgName;
    int length, endFileName;

    imgName = fullImgPath;
    length = strlen(imgName);

    while (imgName[length] != '.')
    {
        length--;
    }
    endFileName = length;
    while (imgName[length] != '/' && length > 0)
    {
        length--;
    }
    if (imgName[length] == '/')
        length++;

    imgName += length;
    strncpy(fileName, imgName, endFileName - length);
    fileName[endFileName - length] = '\0';
}

// return true if the file specified
// by the filename exists
bool file_exists(const char *filename)
{
    struct stat buffer;
    return (stat(filename, &buffer) == 0 && (buffer.st_mode & S_IFMT) == S_IFREG) ? true : false;
}

void runDirectory(Args args)
{
    // determine mode : file or path
    struct stat sb;
    char *label_file, *imgScores_file, *recon_file;
    double score_path = 0;
    int numSuperpixels;

    label_file = (char *)malloc(255 * sizeof(char));
    imgScores_file = (char *)malloc(255 * sizeof(char));
    recon_file = (char *)malloc(255 * sizeof(char));

    if (stat(args.img_path, &sb) == -1)
    {
        perror("stat");
        exit(EXIT_SUCCESS);
    }

    int type;
    switch (sb.st_mode & S_IFMT)
    {
    case S_IFDIR:
        printf("Directory processing: ");
        type = 0;
        break;
    case S_IFREG:
        printf("Single file processing\n");
        type = 1;
        break;
    default:
        type = -1;
        break;
    }

    if (type == -1)
        exit(EXIT_SUCCESS);
    else if (type == 1)
    {
        double score_img = 0;
        char fileName[255];

        getImageName(args.img_path, fileName);

        if (args.label_path != NULL)
        {
            struct stat stats; // determine mode : file or path
            if (stat(args.label_path, &stats) == -1)
                printError("runDirectory", "image/directory %s not found.", args.label_path);

            if ((stats.st_mode & S_IFMT) == S_IFDIR)
                sprintf(label_file, "%s/%s.%s", args.label_path, fileName, args.label_ext);
            else
                sprintf(label_file, "%s", args.label_path);
        }
        else
            label_file = NULL;

        if (args.imgScoresPath != NULL && (args.metric == 1 || args.metric == 2))
        {
            struct stat stats; // determine mode : file or path
            if (stat(args.imgScoresPath, &stats) == -1)
                printError("runDirectory", "image/directory %s not found.", args.imgScoresPath);

            if ((stats.st_mode & S_IFMT) == S_IFDIR)
                sprintf(imgScores_file, "%s/%s.png", args.imgScoresPath, fileName);
            else
                sprintf(imgScores_file, "%s", args.imgScoresPath);
        }
        else
        {
            if (args.saveLabels != NULL && (args.metric == 6 || args.metric == 7))
            {
                struct stat stats; // determine mode : file or path
                if (stat(args.saveLabels, &stats) == -1)
                    printError("runDirectory", "image/directory %s not found.", args.saveLabels);

                if ((stats.st_mode & S_IFMT) == S_IFDIR)
                    sprintf(imgScores_file, "%s/%s.pgm", args.saveLabels, fileName);
                else
                    sprintf(imgScores_file, "%s", args.saveLabels);
            }
            else
                imgScores_file = NULL;
        }

        if (args.imgRecon != NULL)
        {
            struct stat stats; // determine mode : file or path
            if (stat(args.imgRecon, &stats) == -1)
                printError("runDirectory", "image %s not found.", args.imgRecon);

            if ((stats.st_mode & S_IFMT) == S_IFDIR)
                sprintf(recon_file, "%s/%s.png", args.imgRecon, fileName);
            else
                sprintf(recon_file, "%s", args.imgRecon);
        }
        else
            recon_file = NULL;

        score_img = computeMetric(args.img_path, label_file, imgScores_file, recon_file, args, &numSuperpixels);

        // ************************

        if (args.dLogFile != NULL)
        {
            bool file_exist = file_exists(args.dLogFile);
            FILE *fp = fopen(args.dLogFile, "a+");

            if (!file_exist)
            {
                if (args.metric < 6)
                    fprintf(fp, "Image Superpixels Score\n");
                else
                {
                    if (args.metric == 6)
                        fprintf(fp, "Image Superpixels ConnectedSpx\n");
                    else 
                        fprintf(fp, "Image DesiredSpx Superpixels\n");
                }
            }

            if (args.metric == 7)
                fprintf(fp, "%s %d %d\n", fileName, args.k, numSuperpixels);
            else
                fprintf(fp, "%s %d %.5f\n", fileName, numSuperpixels, score_img);
            fclose(fp);
        }

        if (args.metric < 6)
            printf("Score: %.5f , superpixels: %d \n", score_img, numSuperpixels);
        else{
            if (args.metric == 6)
                printf("Superpixels: %d , Connected superpixels: %.5f \n", numSuperpixels, score_img);
            else 
                printf("Desired superpixels: %d , Generated superpixels: %d \n", args.k, numSuperpixels);
        }
    }
    else if (type == 0)
    {
        // get file list
        struct dirent **namelist;
        int n = scandir(args.img_path, &namelist, &filterDir, alphasort);
        if (n == -1)
        {
            printf("No images found.\n");
            perror("scandir");
            exit(EXIT_FAILURE);
        }

        if (n == 0){
            printf("No images found.\n");
            exit(EXIT_SUCCESS);
        }else
            printf("%d Images found.\n", n);

        // process file list
        char *image_name = (char *)malloc(255);

        int numImages = n;
        int sum_num_superpixel = 0;

        while (n--)
        {
            char fileName[255];
            double score_img = 0;

            // get image name
            sprintf(image_name, "%s/%s", args.img_path, namelist[n]->d_name);

            getImageName(image_name, fileName);

            if (args.label_path != NULL)
            {
                struct stat stats; // determine mode : file or path
                if (stat(args.label_path, &stats) == -1)
                    printError("runDirectory", "directory %s not found.", args.label_path);

                if ((stats.st_mode & S_IFMT) == S_IFDIR)
                    sprintf(label_file, "%s/%s.%s", args.label_path, fileName, args.label_ext);
                else
                    printError("runDirectory", "%s must be a directory.", args.label_path);
            }
            else
                label_file = NULL;

            if (args.imgScoresPath != NULL && (args.metric == 1 || args.metric == 2))
            {
                struct stat stats; // determine mode : file or path
                if (stat(args.imgScoresPath, &stats) == -1)
                    printError("runDirectory", "directory %s not found.", args.imgScoresPath);

                if ((stats.st_mode & S_IFMT) == S_IFDIR)
                    sprintf(imgScores_file, "%s/%s.png", args.imgScoresPath, fileName);
                else
                    printError("runDirectory", "%s must be a directory.", args.imgScoresPath);
            }
            else
            {
                if (args.saveLabels != NULL && (args.metric == 6 || args.metric == 7))
                {
                    struct stat stats; // determine mode : file or path
                    if (stat(args.saveLabels, &stats) == -1)
                        printError("runDirectory", "directory %s not found.", args.saveLabels);

                    if ((stats.st_mode & S_IFMT) == S_IFDIR)
                        sprintf(imgScores_file, "%s/%s.pgm", args.saveLabels, fileName);
                    else
                        printError("runDirectory", "%s must be a directory.", args.saveLabels);
                }
                else
                    imgScores_file = NULL;
            }

            if (args.imgRecon != NULL)
            {
                struct stat stats; // determine mode : file or path
                if (stat(args.imgRecon, &stats) == -1)
                    printError("runDirectory", "directory %s not found.", args.imgRecon);

                if ((stats.st_mode & S_IFMT) == S_IFDIR)
                    sprintf(recon_file, "%s/%s.png", args.imgRecon, fileName);
                else
                    printError("runDirectory", "%s must be a directory.", args.imgScoresPath);
            }
            else
                recon_file = NULL;

            // ********
            score_img = computeMetric(image_name, label_file, imgScores_file, recon_file, args, &numSuperpixels);
            // ********

            sum_num_superpixel += numSuperpixels;
            score_path += score_img;

            if (args.dLogFile != NULL)
            {
                bool file_exist = file_exists(args.dLogFile);
                FILE *fp = fopen(args.dLogFile, "a+");

                if (!file_exist)
                {
                    if (args.metric < 6)
                        fprintf(fp, "Image Superpixels Score\n");
                    else
                    {
                        if (args.metric == 6)
                            fprintf(fp, "Image Superpixels ConnectedSpx\n");
                        else 
                            fprintf(fp, "Image DesiredSpx Superpixels\n");
                    }
                }
                if (args.metric == 7)
                    fprintf(fp, "%s %d %d\n", fileName, args.k, numSuperpixels);
                else
                    fprintf(fp, "%s %d %.5f\n", fileName, numSuperpixels, score_img);
                fclose(fp);
            }

            free(namelist[n]);
        }

        free(namelist);

        if (args.logFile != NULL)
        {
            bool file_exist = file_exists(args.logFile);
            FILE *fp = fopen(args.logFile, "a+");

            if (!file_exist)
            {
                if (args.metric < 6)
                    fprintf(fp, "Superpixels Score\n");
                else
                {
                    if (args.metric == 6)
                        fprintf(fp, "Superpixels ConnectedSpx\n");
                    else 
                        fprintf(fp, "DesiredSpx Superpixels\n");
                }
            }

            if (args.metric == 7)
                fprintf(fp, "%d %.5f\n", args.k, (double)sum_num_superpixel / (double)numImages);
            else
                fprintf(fp, "%.5f %.5f\n", (double)sum_num_superpixel / (double)numImages, score_path / (double)numImages);
            fclose(fp);
        }

        free(image_name);
    }

    free(label_file);
    free(imgScores_file);
    free(recon_file);
}

int main(int argc, char *argv[])
{

#ifdef DEBUG
    printf("DEGUB true\n");
#endif

    Args args;
    if (initArgs(&args, argc, argv))
        runDirectory(args);
    else
        usage();

    return 0;
}
