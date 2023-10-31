#include "Eval.h"

int getNumSuperpixels(int *L, int num_pixels, int maxLabel)
{
    int tmp[iftMax(num_pixels, maxLabel)], numLabels = 0;

    for (int i = 0; i < iftMax(num_pixels, maxLabel); i++)
        tmp[i] = 0;
    for (int i = 0; i < num_pixels; i++)
        tmp[L[i]] = -1;
    for (int i = 0; i < iftMax(num_pixels, maxLabel); i++)
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

// return the number of connected components and change labels to have a unique label for each connected component
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
    PrioQueue *queue;

    adj_rel = create8NeighAdj();
    int numSpx = relabelSuperpixels(labels, num_rows, num_cols, 8); // ensure connectivity

    float meanColor[numSpx][3];
    bool **neighbors = (bool**)malloc(numSpx * sizeof(bool*));
    double *sizeSpx = (double *)calloc(numSpx, sizeof(double));

    int new_labels[numSpx];

    for (int i = 0; i < numSpx; i++)
    {
        //sizeSpx[i] = 0;
        meanColor[i][0] = meanColor[i][1] = meanColor[i][2] = 0;
        neighbors[i] = (bool*)calloc(numSpx, sizeof(bool));
        //neighbors[i] = std::vector<int>(0);
        //ids[i] = i;
        new_labels[i] = i;
    }

    // for each superpixel find its adjacent superpixels and compute the mean color and size
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
                    if (neighbors[label][adjLabel] == false)
                    {
                        neighbors[label][adjLabel] = true;
                        neighbors[adjLabel][label] = true;
                    }
                }
            }
        }
    }

    queue = createPrioQueue(numSpx, sizeSpx, MINVAL_POLICY);
    for (int i = 0; i < numSpx; i++) insertPrioQueue(&queue, i);
    
    freeNodeAdj(&adj_rel);

    int k = numSpx - numDesiredSpx; // number of superpixels to merge with other

    // for each superpixel in ascending order of size, merge it with the most similar adjacent superpixel
    while (k > 0 && !isPrioQueueEmpty(queue))
    {
        int label = popPrioQueue(&queue);
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

        bool *ptr_adj = neighbors[label];
        for (int adj = 0; adj < numSpx; adj++, ptr_adj++)
        {
            if(*(ptr_adj) == false) continue;

            float adjMean[3];
            int adjLabel = *(ptr_adj);
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

                float distance = (localMean[0] - adjMean[0]) * (localMean[0] - adjMean[0]) + 
                                (localMean[1] - adjMean[1]) * (localMean[1] - adjMean[1]) + 
                                (localMean[2] - adjMean[2]) * (localMean[2] - adjMean[2]);

                if (distance < minDistance)
                {
                    minDistance = distance;
                    new_label = adjLabel;
                }
            }
            else
            {
                *(ptr_adj) = false;
                neighbors[adj][label] = false;
            }
        }

        // merge with the most similar neighbor superpixel with different label
        if (new_label != -1)
        {
            new_labels[label] = new_label;

            meanColor[new_label][0] += meanColor[label][0];
            meanColor[new_label][1] += meanColor[label][1];
            meanColor[new_label][2] += meanColor[label][2];
            sizeSpx[new_label] += sizeSpx[label];

            neighbors[new_label][label] = false;
            neighbors[label][new_label] = false;
            //remove(neighbors[new_label].begin(), neighbors[new_label].end(), label);

            bool *ptr_adj = neighbors[label];
            bool *ptr_new = neighbors[new_label];
            for (int adj = 0; adj < numSpx; adj++, ptr_adj++, ptr_new++)
            {
                if(adj != label && adj != new_label && *ptr_adj) (*ptr_new) = true;
            }
            moveIndexDownPrioQueue(&queue, new_label);
            k--;
        }
    }

    for (k = 0; k < numSpx; k++)
    {
        free(neighbors[k]);
    }
    free(neighbors);
    free(sizeSpx);

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
    }

    // rotula os pixels com os rótulos novos
    for (int i = 0; i < num_rows * num_cols; i++)
    {
        labels[i] = new_labels[labels[i]];
    }

    return relabelSuperpixels(labels, num_rows, num_cols, 8);
}

//==========================================================

void RBD(Image *image, int *labels, int label, int nbuckets, int *alpha, double **Descriptor)
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
                Descriptor[a][c] = (double)ColorHistogram[hist][bin][c] / (double)V[val]; // get the mean color
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
// COLOR HOMOGENEITY MEASURES
//==========================================================

void RGBtoYCbCr(int* cin, int num_channels, int normalization_value, int *cout)
{
    float a = (16.0/255.0)*(float)normalization_value;
    float b = (128.0/255.0)*(float)normalization_value;

    cout[0]=(int)(0.256789062*(float)cin[0]+
                      0.504128906*(float)cin[1]+
                      0.09790625*(float)cin[2]+a);
    cout[1]=(int)(-0.148222656*(float)cin[0]+
                      -0.290992187*(float)cin[1]+
                      0.439214844*(float)cin[2]+b);
    cout[2]=(int)(0.439214844*(float)cin[0]+
                      -0.367789063*(float)cin[1]+
                      -0.071425781*(float)cin[2]+b);

    for(int i=0; i < 3; i++) {
        if (cout[i] < 0) cout[i] = 0;
        if (cout[i] > normalization_value) cout[i] = normalization_value;
    }
}


double *SIRS(int *labels, Image *image, int alpha, int nbuckets, char *reconFile, double gauss_variance, double *score)
{
    double *histogramVariation;
    double ***Descriptor;  // Descriptor[numSup][alpha][num_channels]
    int *descriptor_size; // descriptor_size[numSup]
    int emptySuperpixels;
    double **MSE;
    //Mat recons;
    iftImage *recons;
    double **mean_buckets;         // mean_buckets[superpixels][num_channels];
    double **variation_descriptor; // variation_descriptor[superpixels][num_channels];

    (*score) = 0;
    emptySuperpixels = 0;

    //recons = Mat::zeros(image->num_rows, image->num_cols, CV_8UC3);
    if(reconFile != NULL) {
        recons = iftCreateColorImage(image->num_cols, image->num_rows, 1, 8); // for RGB colors, depth = 8
        //recons = Mat::zeros(image->num_rows, image->num_cols, CV_8UC3);
    }

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
    Descriptor = (double ***)calloc(superpixels, sizeof(double **));
    MSE = (double **)calloc(superpixels, sizeof(double *));
    mean_buckets = (double **)calloc(superpixels, sizeof(double *));
    variation_descriptor = (double **)calloc(superpixels, sizeof(double *));

    //std::vector<int> superpixelSize(superpixels, 0);
    int *superpixelSize = (int *)calloc(superpixels, sizeof(int));
    for (int s = 0; s < superpixels; s++) superpixelSize[s] = 0;

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
        Descriptor[s] = (double **)calloc(alpha, sizeof(double *));
        MSE[s] = (double *)calloc(image->num_channels, sizeof(double));
        mean_buckets[s] = (double *)calloc(image->num_channels, sizeof(double));
        variation_descriptor[s] = (double *)calloc(image->num_channels, sizeof(double));

        for (int c = 0; c < alpha; c++)
            Descriptor[s][c] = (double *)calloc(image->num_channels, sizeof(double));

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
                variation_descriptor[s][c] = iftMax(variation_descriptor[s][c], abs(((double)Descriptor[s][a][c] / 255.0) - mean_buckets[s][c]));
        }
    }

#ifdef DEBUG
    printf("compute variation \n");
#endif
    int color_in[image->num_channels];
    int color_out[image->num_channels];
    
    for (int i = 0; i < image->num_pixels; i++)
    {
        int label;
        //NodeCoords coords;
        double minVariance;
        int descIndex;

        label = labels[i];
        //coords = getNodeCoordsImage(image->num_cols, i);
        //Vec3b &recons_color = recons.at<Vec3b>(coords.y, coords.x);

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
            for (int c = 0; c < image->num_channels; c++){
                color_in[c] = (int)Descriptor[label][descIndex][c];
            }

            RGBtoYCbCr(color_in, image->num_channels, 255, color_out);
            recons->val[i]=color_out[0];
            recons->Cb[i] =(ushort)color_out[1];
            recons->Cr[i] =(ushort)color_out[2];
            //recons_color[image->num_channels - 1 - c] = Descriptor[label][descIndex][c];
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

    if (reconFile != NULL){
        iftWriteImageByExt(recons, reconFile);
        iftDestroyImage(&recons);
        //imwrite(reconFile, recons);
    }

    return histogramVariation;
}

double *computeExplainedVariation(int *labels, Image *image, char *reconFile, double *score)
{
    (*score) = 0;
    
    double *supExplainedVariation;
    double *valuesTop, *valuesBottom;
    iftImage *recons;
    //Mat recons;

    // get the higher label
    int superpixels = 0;
    for (int i = 0; i < image->num_pixels; ++i)
    {
        if (labels[i] > superpixels)
            superpixels = labels[i];
    }
    superpixels++;

    if(reconFile != NULL) {
        recons = iftCreateColorImage(image->num_cols, image->num_rows, 1, 8); // for RGB colors, depth = 8
        //recons = Mat::zeros(image->num_rows, image->num_cols, CV_8UC3);
    }

    supExplainedVariation = (double *)calloc(superpixels, sizeof(double));
    valuesTop = (double *)calloc(superpixels, sizeof(double));
    valuesBottom = (double *)calloc(superpixels, sizeof(double));

    float **mean = (float**)malloc(superpixels*sizeof(float*));
    int *count = (int*)calloc(superpixels, sizeof(int));
    for(int i=0; i<superpixels; i++){
        mean[i] = (float*)calloc(3, sizeof(float));
    }

    //std::vector<cv::Vec3f> mean(superpixels, cv::Vec3f(0, 0, 0));
    //std::vector<cv::Vec3f> squared_mean(superpixels, cv::Vec3f(0, 0, 0));
    //std::vector<int> count(superpixels, 0);

    float *overall_mean = (float*)calloc(3, sizeof(float));
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
    free(count);

    overall_mean[0] /= image->num_rows * image->num_cols;
    overall_mean[1] /= image->num_rows * image->num_cols;
    overall_mean[2] /= image->num_rows * image->num_cols;

    int color_in[image->num_channels];
    int color_out[image->num_channels];
    for (int i = 0; i < image->num_pixels; ++i)
    {
        int label;
        //NodeCoords coords;

        label = labels[i];
        //coords = getNodeCoordsImage(image->num_cols, i);
        //Vec3b &recons_color = recons.at<Vec3b>(coords.y, coords.x);

        for (int c = 0; c < image->num_channels; ++c)
        {
            valuesTop[label] += (mean[label][c] - overall_mean[c]) * (mean[labels[i]][c] - overall_mean[c]);
            valuesBottom[label] += (image->val[i][c] - overall_mean[c]) * (image->val[i][c] - overall_mean[c]);
            if (reconFile != NULL){
                for (int c = 0; c < image->num_channels; c++){
                    color_in[c] = (int)mean[labels[i]][c];
                }

                RGBtoYCbCr(color_in, image->num_channels, 255, color_out);
                recons->val[i]=color_out[0];
                recons->Cb[i] =(ushort)color_out[1];
                recons->Cr[i] =(ushort)color_out[2];
                //recons_color[image->num_channels - 1 - c] = mean[labels[i]][c];

            }
        }
    }

    free(mean);

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
        iftWriteImageByExt(recons, reconFile);
        //imwrite(reconFile, recons);

    // return sum_top/sum_bottom;
    return supExplainedVariation;
}

//==========================================================
// DELINEATION MEASURES
//==========================================================

bool is4ConnectedBoundaryPixel(int *labels, int num_rows, int num_cols, int i, int j)
{
    int label_tmp = labels[getIndexImage(num_cols, i, j)];

    if (i > 0 && label_tmp != labels[getIndexImage(num_cols, i - 1, j)]) return true;
    if (i < num_rows - 1 && label_tmp != labels[getIndexImage(num_cols, i + 1, j)]) return true;
    if (j > 0 && label_tmp != labels[getIndexImage(num_cols, i, j - 1)]) return true;
    if (j < num_cols - 1 && label_tmp != labels[getIndexImage(num_cols, i, j + 1)]) return true;

    return false;
}

bool is4ConnectedBoundaryPixel2(iftImage *gt, int i, int j)
{
    int num_cols = gt->xsize;
    int label_tmp = gt->val[getIndexImage(num_cols, i, j)];

    if (i > 0 && label_tmp != gt->val[getIndexImage(num_cols, i - 1, j)]) return true;    
    if (i < gt->ysize - 1 && label_tmp != gt->val[getIndexImage(num_cols, i + 1, j)]) return true;
    if (j > 0 && label_tmp != gt->val[getIndexImage(num_cols, i, j - 1)]) return true;
    if (j < num_cols - 1 && label_tmp != gt->val[getIndexImage(num_cols, i, j + 1)]) return true;

    return false;
}

void computeIntersectionMatrix(int *labels, iftImage *gt,
                               int **intersection_matrix, int *superpixel_sizes, int *gt_sizes, int cols, int rows)
{
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j) intersection_matrix[i][j] = 0;
    }

    for (int i = 0; i < gt->ysize; ++i)
    {
        for (int j = 0; j < gt->xsize; ++j)
        {
            int index = getIndexImage(gt->xsize, i, j);
            intersection_matrix[gt->val[index]][labels[index]]++;
            superpixel_sizes[labels[index]]++;
            gt_sizes[gt->val[index]]++;
        }
    }
}

double computeBoundaryRecall(int *labels, iftImage *gt, float d)
{
    int H = gt->ysize; // num_rows
    int W = gt->xsize; // num_cols

    int r = round(d * sqrt(H * H + W * W));

    float tp = 0;
    float fn = 0;

    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            if (is4ConnectedBoundaryPixel2(gt, i, j))
            {
                bool pos = false;
                for (int k = iftMax(0, i - r); k < iftMax(H - 1, i + r) + 1; k++){
                    for (int l = iftMax(0, j - r); l < iftMax(W - 1, j + r) + 1; l++){
                        if (is4ConnectedBoundaryPixel(labels, H, W, k, l)) pos = true;
                    }
                }

                if (pos) tp++;
                else fn++;
            }
        }
    }

    if (tp + fn > 0) return tp / (tp + fn);

    return 0;
}

double computeUndersegmentationError(int *labels, iftImage *gt)
{
    int N = gt->xsize * gt->ysize;
    int **intersection_matrix;
    int *superpixel_sizes;
    int *gt_sizes;

    int superpixels = 0;
    int gt_segments = 0;
    for (int i = 0; i < gt->n; ++i)
    {
        if (labels[i] > superpixels) superpixels = labels[i];
        if (gt->val[i] > gt_segments) gt_segments = gt->val[i];
    }

    superpixels++;
    gt_segments++;

    superpixel_sizes = (int*)calloc(superpixels, sizeof(int));
    gt_sizes = (int*)calloc(superpixels, sizeof(int));

    intersection_matrix = (int**)malloc(gt_segments * sizeof(int*));
    for (int i = 0; i < gt_segments; ++i)
      intersection_matrix[i] = (int*)calloc(superpixels, sizeof(int));
    

    computeIntersectionMatrix(labels, gt, intersection_matrix, superpixel_sizes, gt_sizes, superpixels, gt_segments);

    free(gt_sizes);

    float error = 0;
    for (int j = 0; j < superpixels; ++j)
    {
        int min = 0;
        for (int i = 0; i < gt_segments; ++i)
        {
            int superpixel_j_minus_gt_i = superpixel_sizes[j] - intersection_matrix[i][j];

            if (superpixel_j_minus_gt_i < 0)
            {
                printf("Invalid intersection computed, set difference is negative! \n");
                return -1;
            }
            if (i == 0 || superpixel_j_minus_gt_i < min) min = superpixel_j_minus_gt_i;
        }
        error += min;
    }

    free(superpixel_sizes);
    
    for (int i = 0; i < gt_segments; ++i)
      free(intersection_matrix[i]);
    free(intersection_matrix);

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

    float *perimeter = (float*)calloc(superpixels, sizeof(float));
    float *area = (float*)calloc(superpixels, sizeof(float));

    //std::vector<float> perimeter(superpixels, 0);
    //std::vector<float> area(superpixels, 0); // = number of pixels

    for (int i = 0; i < num_rows; i++)
    {
        for (int j = 0; j < num_cols; j++)
        {
            int count = 0;
            int label = labels[getIndexImage(num_cols, i, j)];
            if (i > 0){
                if (label != labels[getIndexImage(num_cols, i - 1, j)]) count++;
            }
            else count++;

            if (i < num_rows - 1){
                if (label != labels[getIndexImage(num_cols, i + 1, j)]) count++;
            }
            else count++;

            if (j > 0){
                if (label != labels[getIndexImage(num_cols, i, j - 1)]) count++;
            }
            else count++;

            if (j < num_cols - 1){
                if (label != labels[getIndexImage(num_cols, i, j + 1)]) count++;
            }
            else count++;

            perimeter[label] += count;
            area[label] += 1;
        }
    }

    float compactness = 0;

    for (int i = 0; i < superpixels; ++i)
    {
        if (perimeter[i] > 0)
        {
            compactness += area[i] * (4 * IFT_PI * area[i]) / (perimeter[i] * perimeter[i]);
        }
    }

    free(perimeter);
    free(area);

    compactness /= (num_rows * num_cols);

    return compactness;
}


