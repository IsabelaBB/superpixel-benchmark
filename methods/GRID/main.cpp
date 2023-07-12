
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

typedef struct
{
    int num_cols, num_rows, num_channels, num_pixels;
    int **val; // Access by val[i < num_pixels][f < num_channels]
} Image;

void usage();
Image *loadImage(const char *filepath);
void writeImagePGM(Image *img, char *filepath);

typedef struct Args
{
    char *img_path, *label_path;
    char *timeFile, *dTimeFile;
    int k;
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

void printError(const char *function_name, const char *message, ...)
{
    va_list args;
    char full_msg[4096];

    va_start(args, message);
    vsprintf(full_msg, message, args);
    va_end(args);

    fprintf(stderr, "\nError in %s:\n%s!\n", function_name, full_msg);
    fflush(stdout);
    exit(-1);
}

Image *createImage(int num_rows, int num_cols, int num_channels)
{
    Image *new_img;

    new_img = (Image *)calloc(1, sizeof(Image));

    new_img->num_rows = num_rows;
    new_img->num_cols = num_cols;
    new_img->num_pixels = num_rows * num_cols;
    new_img->num_channels = num_channels;

    new_img->val = (int **)calloc(new_img->num_pixels, sizeof(int *));
    for (int i = 0; i < new_img->num_pixels; i++)
        new_img->val[i] = (int *)calloc(num_channels, sizeof(int));

    return new_img;
}

bool initArgs(Args *args, int argc, char *argv[])
{
    char *k_char;

    args->img_path = parseArgs(argv, argc, "--img");
    args->label_path = parseArgs(argv, argc, "--label");
    args->timeFile = parseArgs(argv, argc, "--time");
    args->dTimeFile = parseArgs(argv, argc, "--dtime");

    k_char = parseArgs(argv, argc, "--k");

    if (strcmp(args->timeFile, "-") == 0)
        args->timeFile = NULL;

    if (strcmp(args->dTimeFile, "-") == 0)
        args->dTimeFile = NULL;

    if (strcmp(args->label_path, "-") == 0)
        args->label_path = NULL;

    args->k = strcmp(k_char, "-") != 0 ? atoi(k_char) : 500;

    if (strcmp(args->img_path, "-") == 0 || args->k < 1)
        return false;

    return true;
}

void usage()
{
    printf("Usage: GRID args [optional args] \n");
    printf("----------------------------------\n");
    printf("args:\n");
    printf("--img - Image (STB's supported formats) file/directory\n");
    printf("--k  - Number of superpixels\n");
    printf("optional args:\n");
    printf("--label  - File/directory for final segmentation (pgm format)\n");
    printf("--time   - txt file for time log (mean time for all directory) -- only for a directory of images\n");
    printf("--dtime  - txt file for time log (time at each image) -- only for a directory of images\n");
    printError("main", "Too many/few parameters");
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

Image *loadImage(const char *filepath)
{
    int num_channels, num_rows, num_cols;
    unsigned char *data;
    Image *new_img;

    data = stbi_load(filepath, &num_cols, &num_rows, &num_channels, 0);

    if (data == NULL)
        printError("loadImage", "Could not load the image <%s>", filepath);

    new_img = createImage(num_rows, num_cols, num_channels);

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

//==========================================================

int *GridSegmentation(Image *image, int num_superpixels)
{
    double width;
    double height;

    int *labels = (int *)malloc(image->num_pixels * sizeof(int));

    // find the superpixel size
    width = sqrt((double)(num_superpixels * image->num_cols) / (double)image->num_rows);
    height = ((double)image->num_rows / (double)image->num_cols) * width;

    int final_height, final_width;

    if ((int)((image->num_cols / ceil(width)) * (image->num_rows / ceil(height))) <= num_superpixels)
    {
        final_width = (int)(image->num_cols / ceil(width));
        final_height = (int)(image->num_rows / ceil(height));
    }
    else
    {
        if ((int)((image->num_cols / ceil(width)) * (image->num_rows / floor(height))) <= num_superpixels)
        {
            final_width = (int)(image->num_cols / ceil(width));
            final_height = (int)(image->num_rows / floor(height));
        }
        else
        {
            if ((int)((image->num_cols / floor(width)) * (image->num_rows / ceil(height))) <= num_superpixels)
            {
                final_width = (int)(image->num_cols / floor(width));
                final_height = (int)(image->num_rows / ceil(height));
            }
            else
            {
                final_width = (int)(image->num_cols / floor(width));
                final_height = (int)(image->num_rows / floor(height));
            }
        }
    }

    int last_sup_width = (image->num_cols % final_width);
    int last_sup_height = (image->num_rows % final_height);

    int numSupWeight = image->num_cols / final_width;
    for (int row = 0; row < image->num_rows; row++)
    {
        for (int col = 0; col < image->num_cols; col++)
        {
            if (image->num_cols - col <= last_sup_width)
                labels[image->num_cols * row + col] = labels[image->num_cols * row + col - 1];
            else
            {
                if (image->num_rows - row <= last_sup_height)
                    labels[image->num_cols * row + col] = labels[image->num_cols * (row - 1) + col];
                else
                    labels[image->num_cols * row + col] = (image->num_cols / final_width) * (row / final_height) + (col / final_width);
            }
        }
    }

    return labels;
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
        // ALLOC
        Image *img;
        img = loadImage(args.img_path);
        int *labels;
        double time = 0;

        // ************************
        start = clock();
        labels = GridSegmentation(img, args.k);
        end = clock();
        time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
        // ************************

        char fileName[255];
        getImageName(args.img_path, fileName);

        if (args.label_path != NULL)
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
                writeImagePGMInt(labels, output_file, img->num_cols, img->num_rows);
            }
            else
                writeImagePGMInt(labels, args.label_path, img->num_cols, img->num_rows);
        }

        // write dtime file
        if (args.dTimeFile != NULL)
        {
            FILE *fp = fopen(args.dTimeFile, "a+");
            fprintf(fp, "%s %.5f\n", fileName, time);
            fclose(fp);
        }

        freeImage(&img);
        free(labels);
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
        char *image_name = (char *)malloc(255);
        char *output_file;
        int numImages = n;
        double sum_time = 0;
        if (args.label_path != NULL)
            output_file = (char *)malloc(255);

        while (n--)
        {
            double time = 0;

            // get image name
            sprintf(image_name, "%s/%s", args.img_path, namelist[n]->d_name);

            // alloc structures
            Image *img;
            img = loadImage(image_name);
            int *labels;
            // ********

            // run code
            start = clock();
            labels = GridSegmentation(img, args.k);
            end = clock();
            time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
            sum_time += time;
            // ********

            // get image name
            char fileName[255];
            getImageName(namelist[n]->d_name, fileName);

            // write segmentation image
            if (args.label_path != NULL)
            {
                sprintf(output_file, "%s/%s.pgm", args.label_path, fileName);
                writeImagePGMInt(labels, output_file, img->num_cols, img->num_rows);
            }
            
            // free structures
            free(labels);
            freeImage(&img);

            // write dtime file
            if (args.dTimeFile != NULL)
            {
                FILE *fp = fopen(args.dTimeFile, "a+");
                fprintf(fp, "%s %.5f\n", fileName, time);
                fclose(fp);
            }
            // ********

            free(namelist[n]);
        }

        // write time file
        if (args.timeFile != NULL)
        {
            FILE *fp = fopen(args.timeFile, "a+");
            fprintf(fp, "%d %.6f\n", args.k, sum_time / (double)numImages);
            fclose(fp);
        }

        // Free the rest of structures
        free(image_name);
        if (args.label_path != NULL)
            free(output_file);
        free(namelist);
        // ************************
    }
}

//==========================================================

int main(int argc, char *argv[])
{
    /*
    args:
        image       :   RGB image / image path
        labels      :   Image name (String)
        k           :   Number of superpixels
    */

    Args args;
    if (initArgs(&args, argc, argv))
        runDirectory(args);
    else
        usage();
    return 0;
}
