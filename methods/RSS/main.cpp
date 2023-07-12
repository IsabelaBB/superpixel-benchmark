/**
 * Copyright (c) 2016, Dengfeng Chai
 * Contact: chaidf@zju.edu.cn
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

#include <fstream>
#include <opencv2/opencv.hpp>
#include <bitset>
#include "RSS.h"

using namespace std;
using namespace cv;

/** \brief Command line tool for running RSS.
 * Usage:
 *   $ a.out --help
 *   Allowed options:
 *     --input          arg         the file/folder to process
 *     --label          arg         the file/folder to process
 *     --function       arg         0 diff / 1 range
 *     --k              arg (=400)  number of superpixles
 *     --lambda         arg (=3)    lambda
 *     --perturb-seeds  arg (=0)    perturb seeds: > 0 yes, = 0 no
 *     --parallel       arg (=0)    parallel
 *     --sigma          arg (=1)    sigma used for smoothing (no smoothing if zero)
 *     --time           arg         txt file for time log (mean time for all directory) -- only for a directory of images
 *     --dtime          arg         txt file for time log (time at each image) -- only for a directory of images
 * \endcode
 * \author Dengfeng Chai
 */

typedef struct Args
{
    char *img_path, *label_path;
    char *timeFile, *dTimeFile;
    bool funcRange, perturb_seeds, parallel;
    int superpixels;
    float lambda, sigma;
} Args;

char *parseArgs(char *argv[], int argc, const char *stringKey){
    
    for(int i=1; i < argc; i++){
        // if both strings are identical
        if(strcmp(argv[i], stringKey) == 0){
            return argv[i+1];
        }
    }
    return "-";
}

bool initArgs(Args *args, int argc, char *argv[])
{
    char *superpixels_char, *function_char, *lambda_char, *sigma_char, *perturb_seeds_char, *parallel_char;

    args->img_path = parseArgs(argv, argc, "--img");
    args->label_path = parseArgs(argv, argc, "--label");
    args->timeFile = parseArgs(argv, argc, "--time");
    args->dTimeFile = parseArgs(argv, argc, "--dtime");
    
    superpixels_char = parseArgs(argv, argc, "--k");
    function_char = parseArgs(argv, argc, "--function");
    lambda_char = parseArgs(argv, argc, "--lambda");
    perturb_seeds_char = parseArgs(argv, argc, "--perturb-seeds");
    parallel_char = parseArgs(argv, argc, "--parallel");
    sigma_char = parseArgs(argv, argc, "--sigma");
    
    if (strcmp(args->timeFile, "-") == 0)
        args->timeFile = NULL;

    if (strcmp(args->dTimeFile, "-") == 0)
        args->dTimeFile = NULL;

    if (strcmp(args->label_path, "-") == 0)
        args->label_path = NULL;

    args->superpixels = strcmp(superpixels_char, "-") != 0 ? atoi(superpixels_char) : 200;
    args->lambda = strcmp(lambda_char, "-") != 0 ? atof(lambda_char) : 3;
    args->sigma = strcmp(sigma_char, "-") != 0 ? atof(sigma_char) : 1;
    args->funcRange = strcmp(function_char, "-") != 0 && atoi(function_char) != 0;
    args->perturb_seeds = strcmp(perturb_seeds_char, "-") != 0 && atoi(perturb_seeds_char) != 0;
    args->parallel = strcmp(parallel_char, "-") != 0 && atoi(parallel_char) != 0;

    if (strcmp(args->img_path, "-") == 0 || args->superpixels < 1)
        return false;

    return true;
}

void usage()
{
    printf("Usage: RSS_demo args [optional args] \n");
    printf("----------------------------------\n");
    printf("args:\n");
    printf("--img - Image (STB's supported formats) file/directory\n" );
    printf("--k   - Number of superpixels\n");
    printf("optional args:\n");
    printf("--function      - 0 diff / 1 range\n");
    printf("--lambda        - default:3 \n");
    printf("--perturb-seeds - > 0 yes, = 0 no -- Default:0");
    printf("--parallel      - > 0 yes, = 0 no -- Default:0");
    printf("--sigma         - sigma used for smoothing (no smoothing if zero). Default: 1");
    printf("--label         - File/directory for final segmentation (pgm format)\n");
    printf("--time          - txt file for time log (mean time for all directory) -- only for a directory of images\n");
    printf("--dtime         - txt file for time log (time at each image) -- only for a directory of images\n");
}


void Label2Color(int *label, Mat img_fill)
{
    int num = 16 * 16 * 16;
    vector<Vec3b> rgb(16 * 16 * 16);
    int sr = 0, sg = 0, sb = 0;
    for (int r = 0; r < 16; r++, sr += 9)
        for (int g = 0; g < 16; g++, sg += 9)
            for (int b = 0; b < 16; b++, sb += 9)
            {
                sr = sr % 16;
                sg = sg % 16;
                sb = sb % 16;
                uchar vr = sr * 16;
                uchar vg = sg * 16;
                uchar vb = sb * 16;
                rgb[((r * 16 + g) * 16) + b] = Vec3b(vr, vg, vb);
            }
    int idx = 0;
    std::cout << img_fill.rows << " " << img_fill.cols << "\n";
    for (int i = 0; i < img_fill.rows; i++)
        for (int j = 0; j < img_fill.cols; j++)
        {
            int fij = label[idx++];
            fij = fij % num;
            img_fill.at<Vec3b>(i, j) = rgb[fij];
        }
}

void WriteImageP5(int *labels, int num_cols, int num_rows, const char *filename)
{
    FILE *fp = NULL;
    int p, hi, lo;
    uchar *data8 = NULL;
    ushort *data16 = NULL;

    fp = fopen(filename, "wb");

    fprintf(fp, "P5\n");
    fprintf(fp, "%d %d\n", num_cols, num_rows);

    int maximum = 0, minimum = -1;
    for (int i = 0; i < num_cols * num_rows; i++)
    {
        if (labels[i] > maximum)
            maximum = labels[i];
        if (labels[i] < minimum || minimum == -1)
            minimum = labels[i];
    }

    if ((maximum < 256) && (minimum >= 0))
    {
        fprintf(fp, "%d\n", 255);
        data8 = (uchar *)calloc(num_cols * num_rows, sizeof(uchar));
        for (p = 0; p < num_cols * num_rows; p++)
            data8[p] = (uchar)labels[p];

        fwrite(data8, sizeof(uchar), num_cols * num_rows, fp);
        free(data8);
    }
    else if (maximum < 65536)
    {
        fprintf(fp, "%d\n", 65535);
        data16 = (ushort *)calloc(num_cols * num_rows, sizeof(ushort));
        for (p = 0; p < num_cols * num_rows; p++)
            data16[p] = (ushort)labels[p];

        {
#define HI(num) (((num)&0x0000FF00) >> 8)
#define LO(num) ((num)&0x000000FF)
            for (p = 0; p < num_cols * num_rows; p++)
            {
                hi = HI(data16[p]);
                lo = LO(data16[p]);
                fputc(hi, fp);
                fputc(lo, fp);
            }
        }

        free(data16);
    }
    else
    {
        printf("Cannot write image as P5 (%d/%d)", maximum, minimum);
    }
    fclose(fp);
}

void WriteImageP2(int *labels, int num_cols, int num_rows, const char *filename)
{
    FILE *fp = NULL;
    int p;

    int maximum = 0;
    for (int i = 0; i < num_cols * num_rows; i++)
        if (labels[i] > maximum)
            maximum = labels[i];

    fp = fopen(filename, "w");

    fprintf(fp, "P2\n");
    fprintf(fp, "%d %d\n", num_cols, num_rows);
    fprintf(fp, "%d\n", maximum);

    int k = 0;
    for (int i = 0; i < num_cols; i++)
    {
        for (int j = 0; j < num_rows; j++)
        {
            fprintf(fp, "%d ", labels[k]);
            k++;
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

void WriteLabels(int *labels, int num_cols, int num_rows, const char *filename)
{
    int maximum = 0;
    for (int i = 0; i < num_cols * num_rows; i++)
        if (labels[i] > maximum)
            maximum = labels[i];

    if (maximum > 255)
        WriteImageP2(labels, num_cols, num_rows, filename);
    else
        WriteImageP5(labels, num_cols, num_rows, filename);
}

void Label2Color(int *label, vector<Mat> img_fill)
{
    vector<Vec3b> rgb(16 * 16 * 16);
    int sr = 0, sg = 0, sb = 0;
    for (int r = 0; r < 16; r++, sr += 9)
        for (int g = 0; g < 16; g++, sg += 9)
            for (int b = 0; b < 16; b++, sb += 9)
            {
                sr = sr % 16;
                sg = sg % 16;
                sb = sb % 16;
                uchar vr = sr * 16;
                uchar vg = sg * 16;
                uchar vb = sb * 16;
                rgb[((r * 16 + g) * 16) + b] = Vec3b(vr, vg, vb);
            }
    int idx = 0;
    for (int t = 0; t < img_fill.size(); t++)
        for (int i = 0; i < img_fill[t].rows; i++)
            for (int j = 0; j < img_fill[t].cols; j++)
            {
                int fij = label[idx++];
                img_fill[t].at<Vec3b>(i, j) = rgb[fij];
            }
}

void Label2Boundary(int *label, Mat img, Mat img_boundary)
{
    int rows = img.rows;
    int cols = img.cols;
    img.copyTo(img_boundary);

    Mat istaken = Mat::zeros(img.size(), CV_8U);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
        {
            int np(0);
            int fij = label[i * cols + j];
            // cout << fij << endl;
            bool flag = false;
            int l = max(0, j - 1);
            int r = min(cols - 1, j + 1);
            int u = max(0, i - 1);
            int b = min(rows - 1, i + 1);
            for (int ii = u; ii <= b; ii++)
                for (int jj = l; jj <= r; jj++)
                {
                    int fn = label[ii * cols + jj];
                    if (0 == istaken.at<uchar>(ii, jj))
                    {
                        if (fij != fn)
                            np++;
                    }
                }
            if (np > 1)
            {
                istaken.at<uchar>(i, j) = 255;
                img_boundary.at<Vec3b>(i, j) = Vec3b(0, 255, 0);
            }
        }
}


void SaveVideo(string out_path, vector<string> &fileNames, vector<Mat> &imgsgs)
{
    int status;
    status = mkdir(out_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    for (int i = 0; i < fileNames.size(); i++)
    {
        string fileName = fileNames[i];
        Mat imgsg = imgsgs[i];
        string fileoutput = out_path + fileName;
        imwrite(fileoutput, imgsg);
    }
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
    imgName += length;
    if ((*imgName) == '/')
        imgName++;
    strncpy(fileName, imgName, endFileName - length);
    fileName[endFileName - length] = '\0';
}

int filter(const struct dirent *name)
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

void runDirectory(Args args)
{
    // determine mode : file or path
    struct stat sb;
    clock_t start, end;

    if (stat(args.img_path, &sb) == -1)
    {
        perror("stat");
        exit(EXIT_SUCCESS);
    }

    int type;
    switch (sb.st_mode & S_IFMT)
    {
    case S_IFDIR:
        printf("directory processing\n");
        type = 0;
        break;
    case S_IFREG:
        printf("single file processing\n");
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
        // ALLOC
        Mat image = imread(args.img_path);
        Mat labels;
        double time = 0;

        if (args.sigma > 0.01)
        {
            int size = std::ceil(args.sigma * 4) + 1;
            GaussianBlur(image, image, cv::Size(size, size), args.sigma, args.sigma);
        }
        int nrows = image.rows;
        int ncols = image.cols;
        int nbands = image.channels();
        int nchannels = nbands + 1;
        int total = nrows * ncols * nchannels;

        Mat plane[nbands], res[nchannels];

        split(image, plane);
        res[0] = Mat::zeros(image.size(), CV_8U);
        for (int i = 0; i < nbands; i++)
            res[i + 1] = plane[i];

        Mat dst(image.size(), CV_8UC(nchannels));
        merge(res, nchannels, dst);

        uchar *dataND = new uchar[total];
        for (int i = 0; i < nrows; i++)
            memcpy(dataND + i * ncols * nchannels, dst.data + i * dst.step[0], ncols * nchannels);

        labels = Mat::zeros(image.size(), CV_32SC1);
        int *label = (int *)labels.data;
        RSS sp;

        // ************************
        start = clock();
        sp.Superpixels(dataND, label, nrows, ncols, nchannels, args.superpixels, (int)args.funcRange, args.lambda, args.perturb_seeds, args.parallel);
        end = clock();
        time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
        // ************************
        
        char fileName[255];
        getImageName(args.img_path, fileName);

        if(args.label_path != NULL)
        {
            struct stat label_stats; // determine mode : file or path
            if (stat(args.label_path, &label_stats) == -1)
            {
                perror("stat");
                exit(EXIT_SUCCESS);
            }

            if ((label_stats.st_mode & S_IFMT) == S_IFDIR)
            {
                char output_file[255];
                sprintf(output_file, "%s/%s.pgm", args.label_path, fileName);
                WriteLabels((int *)labels.data, image.cols, image.rows, output_file);
            }
            else
                WriteLabels((int *)labels.data, image.cols, image.rows, args.label_path);
        }

        // write dtime file
        if (args.dTimeFile != NULL)
        {
            FILE *fp = fopen(args.dTimeFile, "a+");
            fprintf(fp, "%s %.5f\n", fileName, time);
            fclose(fp);
        }
    }
    else if (type == 0)
    {
        // get file list
        struct dirent **namelist;
        int n = scandir(args.img_path, &namelist, &filter, alphasort);
        if (n == -1)
        {
            perror("scandir");
            exit(EXIT_FAILURE);
        }

        printf(" %i image(s) found\n", n);
        if (n == 0)
            exit(EXIT_SUCCESS);

        // process file list
        char *image_name = (char *)malloc(255);
        char *output_file;
        int numImages = n;
        double sum_time = 0;
        if(args.label_path != NULL) output_file = (char *)malloc(255);
        
        while (n--)
        {
            double time = 0;

            // get image name
            sprintf(image_name, "%s/%s", args.img_path, namelist[n]->d_name);

            // alloc structures
            Mat image = imread(image_name);
            Mat labels;

            if (args.sigma > 0.01)
            {
                int size = std::ceil(args.sigma * 4) + 1;
                GaussianBlur(image, image, cv::Size(size, size), args.sigma, args.sigma);
            }
            int nrows = image.rows;
            int ncols = image.cols;
            int nbands = image.channels();
            int nchannels = nbands + 1;
            int total = nrows * ncols * nchannels;

            Mat plane[nbands], res[nchannels];

            split(image, plane);
            res[0] = Mat::zeros(image.size(), CV_8U);
            for (int i = 0; i < nbands; i++)
                res[i + 1] = plane[i];

            Mat dst(image.size(), CV_8UC(nchannels));
            merge(res, nchannels, dst);

            uchar *dataND = new uchar[total];
            for (int i = 0; i < nrows; i++)
                memcpy(dataND + i * ncols * nchannels, dst.data + i * dst.step[0], ncols * nchannels);

            labels = Mat::zeros(image.size(), CV_32SC1);
            int *label = (int *)labels.data;
            RSS sp;
            // ********

            // run code
            start = clock();
            sp.Superpixels(dataND, label, nrows, ncols, nchannels, args.superpixels, (int)args.funcRange, args.lambda, args.perturb_seeds, args.parallel);
            end = clock();
            time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
            sum_time += time;
            // ********

            // get image name 
            char fileName[255];
            getImageName(namelist[n]->d_name, fileName);

            // write segmentation image
            if(args.label_path != NULL){
                sprintf(output_file, "%s/%s.pgm", args.label_path, fileName);
                WriteLabels((int *)labels.data, image.cols, image.rows, output_file);
            }

            //write dtime file
            if(args.dTimeFile != NULL){
                FILE *fp = fopen(args.dTimeFile, "a+");
                fprintf(fp, "%s %.5f\n", fileName, time);
                fclose(fp);
            }
            // ********

            free(namelist[n]);
        }

        //write time file
        if(args.timeFile != NULL){
            FILE *fp = fopen(args.timeFile, "a+");
            fprintf(fp, "%d %.6f\n", args.superpixels, sum_time / (double)numImages);
            fclose(fp);
        }

        // Free the rest of structures
        free(image_name);
        if(args.label_path != NULL) free(output_file);
        free(namelist);
        // ************************
    }
}

int main(int argc, char *argv[])
{
    Args args;
    if(initArgs(&args, argc, argv)) runDirectory(args);
    else usage(); 

    return 0;
}
