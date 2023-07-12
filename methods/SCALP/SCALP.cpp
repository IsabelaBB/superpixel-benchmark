#include"./include/SCALP.h"
#include "./include/CImg.h"
#include <iostream>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>

using namespace cimg_library;

typedef struct Args
{
    char *img_path, *label_path, *contour_path;
    char *timeFile, *dTimeFile;
    int superpixels;
    float compactness;
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
    char *superpixels_char, *compactness_char;

    args->img_path = parseArgs(argv, argc, "--img");
    args->label_path = parseArgs(argv, argc, "--label");
    args->timeFile = parseArgs(argv, argc, "--time");
    args->dTimeFile = parseArgs(argv, argc, "--dtime");
    args->contour_path = parseArgs(argv, argc, "--contour");
    
    superpixels_char = parseArgs(argv, argc, "--k");
    compactness_char = parseArgs(argv, argc, "--compactness");
    
    if (strcmp(args->timeFile, "-") == 0)
        args->timeFile = NULL;

    if (strcmp(args->dTimeFile, "-") == 0)
        args->dTimeFile = NULL;

    if (strcmp(args->label_path, "-") == 0)
        args->label_path = NULL;

    if (strcmp(args->contour_path, "-") == 0)
        args->contour_path = NULL;

    args->superpixels = strcmp(superpixels_char, "-") != 0 ? atoi(superpixels_char) : 450;
    args->compactness = strcmp(compactness_char, "-") != 0 ? atof(compactness_char) : 0.075;
    
    if (strcmp(args->img_path, "-") == 0 || args->superpixels < 1)
        return false;

    return true;
}

void usage()
{
    printf("Usage: SCALP args [optional args] \n");
    printf("----------------------------------\n");
    printf("args:\n");
    printf("--img          - Image file/directory\n" );
    printf("--k            - Number of superpixlels. Default: 450 \n");
    printf("--compactness  - Compactness. Default: 0.075\n");
    printf("--contour      - Input contour image file\n");
    printf("optional args:\n");
    printf("--label  - File/directory for final segmentation (pgm format)\n");
    printf("--time   - txt file for time log (mean time for all directory) -- only for a directory of images\n");
    printf("--dtime  - txt file for time log (time at each image) -- only for a directory of images\n");
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

    //Default parameters
    float * kernel = NULL;
    int bres_color = 1;
    int pw = 3;
    float gamma = 50;

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
        CImg<int> img_in(args.img_path);
        CImg<int> img = img_in;
        int nCols = img.width();
        int nRows = img.height();
        double time = 0;

        //Contour loading
        float *contour = NULL;
        int bres_contour = 0;

        char fileName[255];
        getImageName(args.img_path, fileName);

        if(args.contour_path != NULL){
            bres_contour = 1;
            char contour_file[255];

            struct stat contour_stats; // determine mode : file or path
            if (stat(args.contour_path, &contour_stats) == -1)
            {
                perror("stat");
                exit(EXIT_SUCCESS);
            }

            if ((contour_stats.st_mode & S_IFMT) == S_IFDIR){
                sprintf(contour_file, "%s/%s.png", args.contour_path, fileName);
            }else
                sprintf(contour_file, "%s", args.contour_path);

            CImg<int> contour_img(contour_file);
            if ( (contour_img.width() != nCols) || (contour_img.height() != nRows) ) {
                printf("Contour dimensions do not match image size!\n");
                exit(EXIT_SUCCESS);
            }
        
            contour = (float *) malloc(nCols*nRows*sizeof(float));
            for (int j=0; j<nCols; j++) {
                for (int i=0; i<nRows; i++)
                    contour[i+j*nRows] = 1 + gamma*(1-expf(-contour_img(j,i)*contour_img(j,i)/(255*255)));
            }
        }

        unsigned char* R,* G,* B;
        unsigned short* label;
        int pixel=nRows*nCols;
        R=new unsigned char[pixel];
        G=new unsigned char[pixel];
        B=new unsigned char[pixel];
        label=new unsigned short[pixel];

        for (int j=0; j<nCols; j++) {
            for (int i=0; i<nRows; i++) {
                if(img.spectrum() == 3){
                    R[i+j*nRows] = img(j,i,0,0);
                    G[i+j*nRows] = img(j,i,0,1);
                    B[i+j*nRows] = img(j,i,0,2);
                }else{
                    R[i+j*nRows] = img(j,i);
                    G[i+j*nRows] = img(j,i);
                    B[i+j*nRows] = img(j,i);
                }    
            }
        }

        // ************************
        start = clock();
        SCALP(R,G,B,nCols,nRows,args.superpixels,args.compactness,label,kernel,bres_color,bres_contour,contour,pw);
        end = clock();
        time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
        // ************************
        
        if(args.label_path != NULL)
        {
            struct stat label_stats; // determine mode : file or path
            if (stat(args.label_path, &label_stats) == -1)
            {
                perror("stat");
                exit(EXIT_SUCCESS);
            }

            int max_sp = 0;
            CImg<int> output(nCols,nRows);
            for (int i=0; i<nCols; i++){
                for (int j=0; j<nRows; j++) {
                    output(i,j) = (int) label[j+i*nRows];
                    if (output(i,j) > max_sp)
                        max_sp = output(i,j);
                }
            }

            if ((label_stats.st_mode & S_IFMT) == S_IFDIR)
            {
                char output_file[255];
                sprintf(output_file, "%s/%s.pgm", args.label_path, fileName);
                output.save(output_file);
            }
            else
                output.save(args.label_path);
        }

        delete [] R;
        delete [] G;
        delete [] B;
        delete [] label;

        if(args.contour_path != NULL) free(contour);

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
        char *output_file, *contour_file;
        int numImages = n;
        double sum_time = 0;
        if(args.label_path != NULL) output_file = (char *)malloc(255);
        if(args.contour_path != NULL) contour_file = (char *)malloc(255);
        
        while (n--)
        {
            double time = 0;

            // get image name
            sprintf(image_name, "%s/%s", args.img_path, namelist[n]->d_name);

            // alloc structures
            CImg<int> img_in(image_name);
            CImg<int> img = img_in;
            int nCols = img.width();
            int nRows = img.height();
            
            // get image name 
            char fileName[255];
            getImageName(namelist[n]->d_name, fileName);

            //Contour loading
            float *contour = NULL;
            int bres_contour = 0;

            if(args.contour_path != NULL){
                bres_contour = 1;
                sprintf(contour_file, "%s/%s.png", args.contour_path, fileName);

                CImg<int> contour_img(contour_file);

                if ( (contour_img.width() != nCols) || (contour_img.height() != nRows) ) {
                    printf("Contour dimensions do not match image size!\n");
                    continue;
                }
        
                contour = (float *) malloc(nCols*nRows*sizeof(float));
                for (int j=0; j<nCols; j++) {
                    for (int i=0; i<nRows; i++)
                        contour[i+j*nRows] = 1 + gamma*(1-expf(-contour_img(j,i)*contour_img(j,i)/(255*255)));
                }
            }
            // ********

            unsigned char* R,* G,* B;
            unsigned short* label;
            int pixel=nRows*nCols;
            R=new unsigned char[pixel];
            G=new unsigned char[pixel];
            B=new unsigned char[pixel];
            label=new unsigned short[pixel];

            for (int j=0; j<nCols; j++) {
                for (int i=0; i<nRows; i++) {
                    if(img.spectrum() == 3){
                        R[i+j*nRows] = img(j,i,0,0);
                        G[i+j*nRows] = img(j,i,0,1);
                        B[i+j*nRows] = img(j,i,0,2);
                    }else{
                        R[i+j*nRows] = img(j,i);
                        G[i+j*nRows] = img(j,i);
                        B[i+j*nRows] = img(j,i);
                    }    
                }
            }

            // run code
            start = clock();
            SCALP(R,G,B,nCols,nRows,args.superpixels,args.compactness,label,kernel,bres_color,bres_contour,contour,pw);
            end = clock();
            time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
            sum_time += time;
            // ********
            
            // write segmentation image
            if(args.label_path != NULL)
            {
                int max_sp = 0;
                CImg<int> output(nCols,nRows);
                for (int i=0; i<nCols; i++){
                    for (int j=0; j<nRows; j++) {
                        output(i,j) = (int) label[j+i*nRows];
                        if (output(i,j) > max_sp)
                            max_sp = output(i,j);
                    }
                }
                sprintf(output_file, "%s/%s.png", args.label_path, fileName);
                output.save(output_file);
            }

            delete [] R;
            delete [] G;
            delete [] B;
            delete [] label;

            if(args.contour_path != NULL) free(contour);

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


int main(int argc, char* argv[])
{
    Args args;
    if(initArgs(&args, argc, argv)) runDirectory(args);
    else usage();

    return 0;
}
