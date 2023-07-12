/* -- Serge Bobbia : serge.bobbia@u-bourgogne.fr -- Le2i 2018
 * This work is distributed for non commercial use only,
 * it implements the IBIS method as described in the CVPM 2018 paper.
 * Read the ibis.h file for options and benchmark instructions
 *
 * This file show how to instanciate the IBIS class
 * You can either provide a file, or a directory, path to segment images
 */

#include <iostream>
#include "ibis.h"
#include <opencv2/opencv.hpp>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <opencv2/highgui.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include <time.h>

using namespace std;
//=================================================================================
/// DrawContoursAroundSegments
///
/// Internal contour drawing option exists. One only needs to comment the if
/// statement inside the loop that looks at neighbourhood.
//=================================================================================
void DrawContoursAroundSegments(
    unsigned int *&ubuff,
    int *&labels,
    const int &width,
    const int &height,
    const unsigned int &color)
{
    const int dx8[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
    const int dy8[8] = {0, -1, -1, -1, 0, 1, 1, 1};

    int sz = width * height;
    vector<bool> istaken(sz, false);
    vector<int> contourx(sz);
    vector<int> contoury(sz);
    int mainindex(0);
    int cind(0);

    for (int j = 0; j < height; j++)
    {
        for (int k = 0; k < width; k++)
        {
            int np(0);
            for (int i = 0; i < 8; i++)
            {
                int x = k + dx8[i];
                int y = j + dy8[i];

                if ((x >= 0 && x < width) && (y >= 0 && y < height))
                {
                    int index = y * width + x;

                    // if( false == istaken[index] )//comment this to obtain internal contours
                    {
                        if (labels[mainindex] != labels[index])
                            np++;
                    }
                }
            }
            if (np > 1)
            {
                contourx[cind] = k;
                contoury[cind] = j;
                istaken[mainindex] = true;
                // img[mainindex] = color;
                cind++;
            }
            mainindex++;
        }
    }

    int numboundpix = cind; // int(contourx.size());
    for (int j = 0; j < numboundpix; j++)
    {
        int ii = contoury[j] * width + contourx[j];
        ubuff[ii] = 0xffffff;

        for (int n = 0; n < 8; n++)
        {
            int x = contourx[j] + dx8[n];
            int y = contoury[j] + dy8[n];
            if ((x >= 0 && x < width) && (y >= 0 && y < height))
            {
                int ind = y * width + x;
                if (!istaken[ind])
                    ubuff[ind] = 0;
            }
        }
    }
}

typedef struct Args
{
    char *img_path, *label_path;
    char *timeFile, *dTimeFile;
    int compacity, superpixels;
} Args;

char *parseArgs(char *argv[], int argc, const char *stringKey)
{

    for (int i = 1; i < argc; i++)
    {
        // if both strings are identical
        if (strcmp(argv[i], stringKey) == 0)
        {
            return argv[i + 1];
        }
    }
    return "-";
}

bool initArgs(Args *args, int argc, char *argv[])
{
    char *cmp_car, *spx_char;

    args->img_path = parseArgs(argv, argc, "--img");
    args->label_path = parseArgs(argv, argc, "--label");
    args->timeFile = parseArgs(argv, argc, "--time");
    args->dTimeFile = parseArgs(argv, argc, "--dtime");

    cmp_car = parseArgs(argv, argc, "--cmp");
    spx_char = parseArgs(argv, argc, "--k");

    if (strcmp(args->timeFile, "-") == 0)
        args->timeFile = NULL;

    if (strcmp(args->dTimeFile, "-") == 0)
        args->dTimeFile = NULL;

    if (strcmp(args->label_path, "-") == 0)
        args->label_path = NULL;

    args->compacity = strcmp(cmp_car, "-") != 0 ? atoi(cmp_car) : 20;
    args->superpixels = strcmp(spx_char, "-") != 0 ? atoi(spx_char) : 200;

    if (strcmp(args->img_path, "-") == 0 || args->compacity < 0 || args->superpixels < 1)
        return false;

    return true;
}

void usage()
{
    printf("Usage: ./ibis args [optional args] \n");
    printf("----------------------------------\n");
    printf("args:\n");
    printf("--img   - Image (STB's supported formats) file/directory\n");
    printf(" --k     - user fixed number of superpixels, > 0\n");
    printf(" --cmp   - Factor of compacity, set to 20 for benchmark, > 0\n");
    printf("optional args:\n");
    printf("--label  - File/directory for final segmentation (pgm format)\n");
    printf("--time   - txt file for time log (mean time for all directory) -- only for a directory of images\n");
    printf("--dtime  - txt file for time log (time at each image) -- only for a directory of images\n");
    exit(EXIT_SUCCESS);
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

std::string get_name(const std::string &path_with_ext)
{
    int deb = path_with_ext.find_last_of("/");
    int fin = path_with_ext.find_last_of(".");
    return path_with_ext.substr(deb + 1, fin - deb - 1);
}

void execute_IBIS(int K, int compa, IBIS *Super_Pixel, cv::Mat *img, char *output_file, double *time)
{

    clock_t start, end;

    start = clock();
    // process IBIS
    Super_Pixel->process(img);
    end = clock();

    (*time) = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);

    // convert int* labels to Mat* labels in gray scale
    int *labels = Super_Pixel->getLabels();

    const int width = img->cols;
    const int height = img->rows;
    const int color = 0xFFFFFFFF;

    if (output_file != NULL)
        WriteLabels(labels, width, height, output_file);
}

int filter(const struct dirent *name)
{
    std::string file_name = std::string(name->d_name);
    std::size_t found = file_name.find(".png");
    if (found != std::string::npos)
    {
        return 1;
    }

    found = file_name.find(".jpg");
    if (found != std::string::npos)
    {
        return 1;
    }

    found = file_name.find(".ppm");
    if (found != std::string::npos)
    {
        return 1;
    }

    found = file_name.find(".bmp");
    if (found != std::string::npos)
    {
        return 1;
    }

    found = file_name.find(".tiff");
    if (found != std::string::npos)
    {
        return 1;
    }

    return 0;
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

int main(int argc, char *argv[])
{
    Args args;

    if (argc == 1)
        usage();

    if (!initArgs(&args, argc, argv))
        usage();

    // determine mode : file or path
    struct stat sb;

    if (stat(args.img_path, &sb) == -1)
    {
        perror("stat");
        exit(EXIT_SUCCESS);
    }

    int type;
    switch (sb.st_mode & S_IFMT)
    {
    case S_IFDIR:
        //printf("directory processing\n");
        type = 0;
        break;
    case S_IFREG:
        //printf("single file processing\n");
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
        double time;

        // IBIS
        IBIS Super_Pixel(args.superpixels, args.compacity);

        // get picture
        cv::Mat img = cv::imread(args.img_path);

        // execute IBIS
        char fileName[255];
        getImageName(args.img_path, fileName);

        if (args.label_path == NULL)
            execute_IBIS(args.superpixels, args.compacity, &Super_Pixel, &img, NULL, &time);
        else
        {
            struct stat label_stats; // determine mode : file or path
            if (stat(args.label_path, &label_stats) == -1)
            {
                perror("stat");
                exit(EXIT_SUCCESS);
            }

            char output_file[255];
            if ((label_stats.st_mode & S_IFMT) == S_IFDIR)
                sprintf(output_file, "%s/%s.pgm", args.label_path, fileName);
            else
                printf(output_file, "%s", args.label_path);
            
            execute_IBIS(args.superpixels, args.compacity, &Super_Pixel, &img, output_file, &time);
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

        //printf(" %i image(s) found\n", n);
        if (n == 0)
            exit(EXIT_SUCCESS);

        // process file list
        int width = 0;
        int height = 0;
        IBIS *Super_Pixel;
        char *image_name = (char *)malloc(255);
        int numImages = n;
        char *output_file;
        double sum_time = 0;

        if (args.label_path != NULL)
            output_file = (char *)malloc(255);

        while (n--)
        {
            double time = 0;

            // get picture
            sprintf(image_name, "%s/%s", args.img_path, namelist[n]->d_name);
            cv::Mat img = cv::imread(image_name);

            // execute IBIS
            if (width == 0)
            {
                width = img.cols;
                height = img.rows;

                // IBIS
                Super_Pixel = new IBIS(args.superpixels, args.compacity);
            }
            else
            {
                if (width != img.cols)
                {
                    delete Super_Pixel;
                    Super_Pixel = new IBIS(args.superpixels, args.compacity);
                    width = img.cols;
                    height = img.rows;
                }
            }

            // get image name
            char fileName[255];
            getImageName(namelist[n]->d_name, fileName);

            if (args.label_path == NULL)
                execute_IBIS(args.superpixels, args.compacity, Super_Pixel, &img, NULL, &time);
            else
            {
                sprintf(output_file, "%s/%s.pgm", args.label_path, fileName);
                execute_IBIS(args.superpixels, args.compacity, Super_Pixel, &img, output_file, &time);
            }
            sum_time += time;

            // write dtime file
            if (args.dTimeFile != NULL)
            {
                FILE *fp = fopen(args.dTimeFile, "a+");
                fprintf(fp, "%s %.5f\n", fileName, time);
                fclose(fp);
            }

            free(namelist[n]);
        }

        // write time file
        if (args.timeFile != NULL)
        {
            FILE *fp = fopen(args.timeFile, "a+");
            fprintf(fp, "%d %.6f\n", args.superpixels, sum_time / (double)numImages);
            fclose(fp);
        }

        // Free the rest of structures
        free(image_name);
        if (args.label_path != NULL)
            free(output_file);
        delete Super_Pixel;
        free(namelist);
    }

    exit(EXIT_SUCCESS);
}
