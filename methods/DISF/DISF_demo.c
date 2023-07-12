/**
* Dynamic and Iterative Spanning Forest (C)
*
* @date September, 2019
*/

//=============================================================================
// Includes
//=============================================================================
#include "Image.h"
#include "DISF.h"
#include "Utils.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

//=============================================================================
// Prototypes
//=============================================================================
void usage();
Image *loadImage(const char* filepath);
void writeImagePGM(Image *img, char* filepath);

typedef struct Args
{
    char *img_path, *label_path;
    char *timeFile, *dTimeFile;
    int n0, nf;
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
    char *n0_char, *nf_char;

    args->img_path = parseArgs(argv, argc, "--img");
    args->label_path = parseArgs(argv, argc, "--label");
    args->timeFile = parseArgs(argv, argc, "--time");
    args->dTimeFile = parseArgs(argv, argc, "--dtime");
    
    n0_char = parseArgs(argv, argc, "--n0");
    nf_char = parseArgs(argv, argc, "--nf");
    
    if (strcmp(args->timeFile, "-") == 0)
        args->timeFile = NULL;

    if (strcmp(args->dTimeFile, "-") == 0)
        args->dTimeFile = NULL;

    if (strcmp(args->label_path, "-") == 0)
        args->label_path = NULL;

    args->n0 = strcmp(n0_char, "-") != 0 ? atoi(n0_char) : 8000;
    args->nf = strcmp(nf_char, "-") != 0 ? atoi(nf_char) : 200;
    
    if (strcmp(args->img_path, "-") == 0 || args->n0 < 1 || args->nf < 1 || args->n0 < args->nf)
        return false;

    return true;
}

void usage()
{
    printf("Usage: DISF_demo args [optional args] \n");
    printf("----------------------------------\n");
    printf("args:\n");
    printf("--img - Image (STB's supported formats) file/directory\n" );
    printf("--n0  - Initial number of seeds (e.g., N0 = 8000)\n");
    printf("--nf  - Final number of superpixels (e.g., Nf = 50)\n");
    printf("optional args:\n");
    printf("--label  - File/directory for final segmentation (pgm format)\n");
    printf("--time   - txt file for time log (mean time for all directory) -- only for a directory of images\n");
    printf("--dtime  - txt file for time log (time at each image) -- only for a directory of images\n");
    printError("main", "Too many/few parameters");
}


int filter(const struct dirent *name)
{
    int pos=0;
    while(name->d_name[pos] != '.'){
        if(name->d_name[pos] == '\0') return 0;
        pos++;
    }
    pos++;
    if(name->d_name[pos] == '\0') return 0;

    int extSize = 0;
    while(name->d_name[pos+extSize] != '\0'){
        extSize++;
    }

    char ext[extSize];
    int pos2=0;
    while(pos2 < extSize){
        ext[pos2] = name->d_name[pos+pos2];
        pos2++;
    }

    if((extSize == 3 && ((ext[0]=='p' && ext[1]=='n' && ext[2]=='g') || 
        (ext[0]=='j' && ext[1]=='p' && ext[2]=='g') || 
        (ext[0]=='p' && ext[1]=='p' && ext[2]=='m') || 
        (ext[0]=='p' && ext[1]=='g' && ext[2]=='m') || 
        (ext[0]=='b' && ext[1]=='m' && ext[2]=='p') )) || 
        (extSize == 4 && ext[0]=='t' && ext[1]=='i' && ext[2]=='f' && ext[2]=='f'))
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
        type = 0;
        break;
    case S_IFREG:
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
        Image *img, *border_img, *label_img;
        Graph *graph;
        double time = 0;

        img = loadImage(args.img_path);
        border_img = createImage(img->num_rows, img->num_cols, 1);
        graph = createGraph(img);
        freeImage(&img);

        start = clock();
        label_img = runDISF(graph, args.n0, args.nf, &border_img);
        end = clock();
        time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
        // ************************
        
        freeGraph(&graph);
        freeImage(&border_img);

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
                writeImagePGM(label_img, output_file);
            }
            else
                writeImagePGM(label_img, args.label_path);
        }

        // write dtime file
        if (args.dTimeFile != NULL)
        {
            FILE *fp = fopen(args.dTimeFile, "a+");
            fprintf(fp, "%s %.5f\n", fileName, time);
            fclose(fp);
        }

        freeImage(&label_img);
        
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
        if(args.label_path != NULL) output_file = (char *)malloc(255);
        
        while (n--)
        {
            double time = 0;

            // get image name
            sprintf(image_name, "%s/%s", args.img_path, namelist[n]->d_name);

            // alloc structures
            Image *img, *border_img, *label_img;
            Graph *graph;
            
            img = loadImage(image_name);
            border_img = createImage(img->num_rows, img->num_cols, 1);
            graph = createGraph(img);
            freeImage(&img);
            // ********

            // run code
            start = clock();
            label_img = runDISF(graph, args.n0, args.nf, &border_img);
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
                writeImagePGM(label_img, output_file);
            }

            // free structures
            freeGraph(&graph);
            freeImage(&label_img);
            freeImage(&border_img);

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
            fprintf(fp, "%d %.6f\n", args.nf, sum_time / (double)numImages);
            fclose(fp);
        }

        // Free the rest of structures
        free(image_name);
        if(args.label_path != NULL) free(output_file);
        free(namelist);
        // ************************
        
    }
}


//=============================================================================
// Main
//=============================================================================
int main(int argc, char* argv[])
{
    Args args;
    if(initArgs(&args, argc, argv)) runDirectory(args);
    else usage();
}

//=============================================================================
// Methods
//=============================================================================

Image *loadImage(const char* filepath)
{
    int num_channels, num_rows, num_cols;
    unsigned char *data;    
    Image *new_img;
    
    data = stbi_load(filepath, &num_cols, &num_rows, &num_channels, 0);

    if(data == NULL)
        printError("loadImage", "Could not load the image <%s>", filepath);

    new_img = createImage(num_rows, num_cols, num_channels);

    #pragma omp parallel for
    for(int i = 0; i < new_img->num_pixels; i++)
    {
        new_img->val[i] = (int*)calloc(new_img->num_channels, sizeof(int));

        for(int j = 0; j < new_img->num_channels; j++) 
        {
            new_img->val[i][j] = data[i * new_img->num_channels + j];
        }
    }

    stbi_image_free(data);

    return new_img;
}

void writeImagePGM(Image *img, char* filepath)
{
    int max_val, min_val;
    FILE *fp;

    fp = fopen(filepath, "wb");

    if(fp == NULL)
        printError("writeImagePGM", "Could not open the file <%s>", filepath);

    max_val = getMaximumValue(img, -1);
    min_val = getMinimumValue(img, -1);

    fprintf(fp, "P5\n");
    fprintf(fp, "%d %d\n", img->num_cols, img->num_rows);
    fprintf(fp, "%d\n", max_val);

    // 8-bit PGM file
    if(max_val < 256 && min_val >= 0)
    {
        unsigned char* data;

        data = (unsigned char*)calloc(img->num_pixels, sizeof(unsigned char));

        for(int i = 0; i < img->num_pixels; i++)
            data[i] = (unsigned char)img->val[i][0];

        fwrite(data, sizeof(unsigned char), img->num_pixels, fp);

        free(data);
    }
    // 16-bit PGM file
    else if(max_val < 65536 && min_val >= 0) 
    {
        unsigned short* data;

        data = (unsigned short*)calloc(img->num_pixels, sizeof(unsigned short));

        for(int i = 0; i < img->num_pixels; i++)
            data[i] = (unsigned short)img->val[i][0];
        
        for(int i = 0; i < img->num_pixels; i++)
        {
            int high, low;

            high = ((data[i]) & 0x0000FF00) >> 8;
            low = (data[i]) & 0x000000FF;

            fputc(high,fp);
            fputc(low, fp);
        }

        free(data);   
    }
    else
        printError("writeImagePGM", "Invalid min/max spel values <%d,%d>", min_val, max_val);

    fclose(fp);
}
