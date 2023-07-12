/*
    ===== DESCRIPTION =====

    This is an implementation of the algorithms in

    Real-Time Coarse-to-fine Topologically Preserving Segmentation
    by Jian Yao, Marko Boben, Sanja Fidler, Raquel Urtasun

    published in CVPR 2015.

    http://www.cs.toronto.edu/~urtasun/publications/yao_etal_cvpr15.pdf

    Code is available from https://bitbucket.org/mboben/spixel

    ===== LICENSE =====

    Copyright(C) 2015  Jian Yao, Marko Boben, Sanja Fidler and Raquel Urtasun

    This program is free software : you can redistribute it and / or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "stdafx.h"
#include "segengine.h"
#include "functions.h"
#include "utils.h"
#include <fstream>
#include <cstdlib>
#include "contrib/SGMStereo.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>

using namespace cv;
using namespace std;

typedef struct Args
{
    char *img_path, *label_path;
    char *timeFile, *dTimeFile, *configFile;
    int superpixels;
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
    char *superpixels_char;

    args->img_path = parseArgs(argv, argc, "--img");
    args->label_path = parseArgs(argv, argc, "--label");
    args->timeFile = parseArgs(argv, argc, "--time");
    args->dTimeFile = parseArgs(argv, argc, "--dtime");
    args->configFile = parseArgs(argv, argc, "--config");

    superpixels_char = parseArgs(argv, argc, "--k");

    if (strcmp(args->timeFile, "-") == 0)
        args->timeFile = NULL;

    if (strcmp(args->dTimeFile, "-") == 0)
        args->dTimeFile = NULL;

    if (strcmp(args->label_path, "-") == 0)
        args->label_path = NULL;

    args->superpixels = strcmp(superpixels_char, "-") != 0 ? atoi(superpixels_char) : 200;

    if (strcmp(args->img_path, "-") == 0 || args->superpixels < 1 || strcmp(args->configFile, "-") == 0)
        return false;

    return true;
}

void usage()
{
    printf("Usage: spixel args [optional args] \n");
    printf("----------------------------------\n");
    printf("args:\n");
    printf("--img     - Image file/directory\n");
    printf("--k       - Number of superpixels \n");
    printf("--config  - Configure file. Example: ./ETPS/examples/basic/config.yml\n");
    printf("optional args:\n");
    printf("--label   - File/directory for final segmentation (pgm format)\n");
    printf("--time    - txt file for time log (mean time for all directory) -- only for a directory of images\n");
    printf("--dtime   - txt file for time log (time at each image) -- only for a directory of images\n");
}

void ConvertOCVToPNG(const Mat &ocvImg, png::image<png::rgb_pixel> &pngImg)
{
    pngImg.resize(ocvImg.cols, ocvImg.rows);

    for (int r = 0; r < ocvImg.rows; r++)
    {
        for (int c = 0; c < ocvImg.cols; c++)
        {
            const Vec3b &bgrColor = ocvImg.at<Vec3b>(r, c);
            pngImg.set_pixel(c, r, png::rgb_pixel(bgrColor[2], bgrColor[1], bgrColor[0]));
        }
    }
}

Mat ConvertFloatToOCV(int width, int height, const float *data)
{
    Mat1w result(height, width);

    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            result(y, x) = (ushort)(data[width * y + x] * 256.0f + 0.5);
        }
    }
    return result;
}

double SGMPreprocessing(const Mat &leftImage, const Mat &rightImage, Mat &dispImage)
{
    png::image<png::rgb_pixel> leftImageSGM, rightImageSGM;

    ConvertOCVToPNG(leftImage, leftImageSGM);
    ConvertOCVToPNG(rightImage, rightImageSGM);

    size_t width = leftImageSGM.get_width();
    size_t height = leftImageSGM.get_height();

    if (width != rightImageSGM.get_width() || height != rightImageSGM.get_height())
    {
        dispImage = Mat1w();
        return 0.0;
    }

    float *dispImageFloat = (float *)malloc(width * height * sizeof(float));

    SGMStereo sgm;
    Timer t;

    sgm.compute(leftImageSGM, rightImageSGM, dispImageFloat);
    t.Stop();
    dispImage = ConvertFloatToOCV(width, height, dispImageFloat);

    free(dispImageFloat);

    return t.GetTimeInSec();
}

void ProcessFilesBatch(SPSegmentationParameters &params, const vector<string> &files, const string &fileDir)
{
    MkDir(fileDir + "out");
    MkDir(fileDir + "seg");

    int nProcessed = 0;
    double totalTime = 0.0;

    for (const string &f : files)
    {
        string fileName = fileDir + f;
        Mat image = imread(fileName, IMREAD_COLOR);

        if (image.rows == 0 || image.cols == 0)
        {
            cout << "Failed reading image '" << fileName << "'" << endl;
            continue;
        }

        cout << "Processing: " << fileName << endl;

        SPSegmentationEngine engine(params, image);

        engine.ProcessImage();

        engine.PrintPerformanceInfo();
        totalTime += engine.ProcessingTime();

        string outImage = ChangeExtension(fileDir + "out/" + f, "_sp.png");
        string outImageSeg = ChangeExtension(fileDir + "seg/" + f, ".png");

        imwrite(outImage, engine.GetSegmentedImage());
        imwrite(outImageSeg, engine.GetSegmentation());

        nProcessed++;
    }

    if (nProcessed > 1 && params.timingOutput)
    {
        cout << "Processed " << nProcessed << " files in " << totalTime << " sec. ";
        cout << "Average per image " << (totalTime / nProcessed) << " sec." << endl;
    }
}

// If dispPattern is empty, dispDir is the whole file name (in case we call this
// function to process only one file)
void ProcessFilesStereoBatch(SPSegmentationParameters &params, const vector<string> &files, const string &fileDir,
                             const string &dispDir, const string &dispPattern)
{
    MkDir(fileDir + "out");
    MkDir(fileDir + "seg");
    MkDir(fileDir + "disp");

    int nProcessed = 0;
    double totalTime = 0.0;

    for (const string &f : files)
    {
        string fileName = fileDir + f;
        string dispFileName = dispPattern.empty() ? dispDir : ChangeExtension(dispDir + f, dispPattern);
        Mat image = imread(fileName, IMREAD_COLOR);
        Mat dispImage = imread(dispFileName, IMREAD_ANYDEPTH);

        if (image.rows == 0 || image.cols == 0)
        {
            cout << "Failed reading image '" << fileName << "'" << endl;
            continue;
        }
        if (dispImage.rows == 0 || dispImage.cols == 0)
        {
            cout << "Failed reading disparity image '" << dispFileName << "'" << endl;
            continue;
        }
        cout << "Processing: " << fileName << endl;

        SPSegmentationEngine engine(params, image, dispImage);

        engine.ProcessImageStereo();
        engine.PrintDebugInfoStereo();
        engine.PrintPerformanceInfo();
        totalTime += engine.ProcessingTime();

        string outImage = ChangeExtension(fileDir + "out/" + f, "_sp.png");
        string outImageSeg = ChangeExtension(fileDir + "seg/" + f, ".png");
        string outImageDisp = ChangeExtension(fileDir + "disp/" + f, ".png");

        imwrite(outImage, engine.GetSegmentedImage());
        imwrite(outImageSeg, engine.GetSegmentation());
        imwrite(outImageDisp, engine.GetDisparity());

        nProcessed++;
    }

    if (nProcessed > 1 && params.timingOutput)
    {
        cout << "Processed " << nProcessed << " files in " << totalTime << " sec. ";
        cout << "Average per image " << (totalTime / nProcessed) << " sec." << endl;
    }
}

void ProcessFilesStereoBatchSGM(SPSegmentationParameters &params, const vector<string> &files, const string &leftFileDir,
                                const string &rightFileDir, bool rightIsName)
{
    MkDir(leftFileDir + "out");
    MkDir(leftFileDir + "seg");
    MkDir(leftFileDir + "disp");

    int nProcessed = 0;
    double totalTime = 0.0;

    for (const string &f : files)
    {
        string leftFileName = leftFileDir + f;
        string rightFileName = rightIsName ? rightFileDir : rightFileDir + f;
        Mat leftImage = imread(leftFileName, IMREAD_COLOR);
        Mat rightImage = imread(rightFileName, IMREAD_COLOR);

        if (leftImage.empty())
        {
            cout << "Failed reading left image '" << leftImage << "'" << endl;
            continue;
        }
        if (rightImage.empty())
        {
            cout << "Failed reading right image '" << rightImage << "'" << endl;
            continue;
        }
        cout << "Processing: " << leftFileName << "/" << rightFileName << endl;

        Mat dispImage;
        double sgmTime;

        sgmTime = SGMPreprocessing(leftImage, rightImage, dispImage);
        totalTime += sgmTime;

        if (dispImage.empty())
        {
            cout << "Failed creating SGM image for '" << leftImage << "'/'" << rightImage << "' pair" << endl;
            continue;
        }

        if (params.timingOutput)
        {
            cout << "SGM processing time: " << sgmTime << " sec." << endl;
        }

        SPSegmentationEngine engine(params, leftImage, dispImage);

        engine.ProcessImageStereo();
        engine.PrintDebugInfoStereo();
        engine.PrintPerformanceInfo();
        totalTime += engine.ProcessingTime();

        string outImage = ChangeExtension(leftFileDir + "out/" + f, "_sp.png");
        string outImageSeg = ChangeExtension(leftFileDir + "seg/" + f, ".png");
        string outImageDisp = ChangeExtension(leftFileDir + "disp/" + f, ".png");

        imwrite(outImage, engine.GetSegmentedImage());
        imwrite(outImageSeg, engine.GetSegmentation());
        imwrite(outImageDisp, engine.GetDisparity());

        nProcessed++;
    }

    if (nProcessed > 1 && params.timingOutput)
    {
        cout << "Processed " << nProcessed << " files in " << totalTime << " sec. ";
        cout << "Average per image " << (totalTime / nProcessed) << " sec." << endl;
    }
}

void ProcessFilesStereoBatchSGM(SPSegmentationParameters &params, string leftDir, string rightDir,
                                const string &filePattern)
{
    vector<string> files;

    FindFiles(leftDir, filePattern, files, false);
    EndDir(leftDir);
    EndDir(rightDir);
    ProcessFilesStereoBatchSGM(params, files, leftDir, rightDir, false);
}

void ProcessFileStereoSGM(SPSegmentationParameters &params, string leftFile, string rightFile,
                          const string &filePattern)
{
    vector<string> files;
    string leftDir = FilePath(leftFile);
    string leftName = FileName(leftFile);

    files.push_back(leftName);
    ProcessFilesStereoBatchSGM(params, files, leftDir, rightFile, true);
}

void ProcessFilesBatch(SPSegmentationParameters &params, const string &dirName, const string &pattern)
{
    vector<string> files;
    string fileDir = dirName;

    FindFiles(fileDir, pattern, files, false);
    EndDir(fileDir);
    ProcessFilesBatch(params, files, fileDir);
}

void ProcessFile(SPSegmentationParameters &params, const string &file, const string &file_out)
{
    int nProcessed = 0;
    double totalTime = 0.0;

    Mat image = imread(file, IMREAD_COLOR);

    if (image.rows == 0 || image.cols == 0)
    {
        cout << "Failed reading image '" << file << "'" << endl;
        return;
    }
    cout << "Processing: " << file << endl;

    SPSegmentationEngine engine(params, image);

    engine.ProcessImage();

    engine.PrintPerformanceInfo();
    totalTime += engine.ProcessingTime();

    imwrite(file_out, engine.GetSegmentation());

    nProcessed++;
}

void ProcessFilesStereoBatch(SPSegmentationParameters &params, const string &dirName, const string &pattern,
                             const string &dispPattern)
{
    vector<string> files;
    string fileDir = dirName;

    FindFiles(fileDir, pattern, files, false);
    EndDir(fileDir);
    ProcessFilesStereoBatch(params, files, fileDir, fileDir, dispPattern);
}

void ProcessFileStereo(SPSegmentationParameters &params, const string &file, const string &dispFile)
{
    vector<string> files;
    string fileDir = FilePath(file);
    string fileName = FileName(file);

    files.push_back(fileName);
    ProcessFilesStereoBatch(params, files, fileDir, dispFile, "");
}

void ProcessFiles(const string &paramFile, const string &name1, const string &name2,
                  const string &name3)
{
    SPSegmentationParameters params = ReadParameters(paramFile);

    if (params.randomSeed > 0)
    {
        srand(params.randomSeed);
    }

    if (params.stereo)
    {
        if (params.batchProcessing)
        {
            if (params.computeSGM)
                ProcessFilesStereoBatchSGM(params, name1, name2, name3);
            else
                ProcessFilesStereoBatch(params, name1, name2, name3);
        }
        else
        { // !batchProcessing
            if (params.computeSGM)
                ProcessFileStereoSGM(params, name1, name2, name3);
            else
                ProcessFileStereo(params, name1, name2);
        }
    }
    else
    { // !stereo
        if (params.batchProcessing)
            ProcessFilesBatch(params, name1, name2);
        else
        {
            params.superpixelNum = atoi(name3.c_str());
            ProcessFile(params, name1, name2);
        }
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
        SPSegmentationParameters params = ReadParameters(args.configFile);
        double time = 0;

        if (params.randomSeed > 0)
        {
            srand(params.randomSeed);
        }

        Mat image = imread(args.img_path, IMREAD_COLOR);
        params.superpixelNum = args.superpixels;

        if (params.randomSeed > 0)
        {
            srand(params.randomSeed);
        }

        if (image.rows == 0 || image.cols == 0)
        {
            cout << "Failed reading image '" << args.img_path << "'" << endl;
            exit(EXIT_SUCCESS);
        }

        SPSegmentationEngine engine(params, image);

        // ************************
        start = clock();
        engine.ProcessImage();
        end = clock();
        time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
        // engine.PrintPerformanceInfo();
        //  ************************

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
                imwrite(output_file, engine.GetSegmentation());
            }
            else
                imwrite(args.label_path, engine.GetSegmentation());
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
        char *image_name = (char *)malloc(255);
        char *output_file;
        int numImages = n;
        clock_t start, end;
        double sum_time = 0;
        if (args.label_path != NULL)
            output_file = (char *)malloc(255);

        while (n--)
        {
            double time = 0;
            SPSegmentationParameters params = ReadParameters(args.configFile);
            if (params.randomSeed > 0)
                srand(params.randomSeed);
            
            params.superpixelNum = args.superpixels;

            // get image name
            sprintf(image_name, "%s/%s", args.img_path, namelist[n]->d_name);

            // alloc structures
            Mat image = imread(image_name, IMREAD_COLOR);

            if (image.rows == 0 || image.cols == 0)
            {
                cout << "Failed reading image '" << image_name << "'" << endl;
                continue;
            }

            SPSegmentationEngine engine(params, image);
            // ********

            // run code
            start = clock();
            engine.ProcessImage();
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
                sprintf(output_file, "%s/%s.png", args.label_path, fileName);
                imwrite(output_file, engine.GetSegmentation());
            }

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
            fprintf(fp, "%d %.6f\n", args.superpixels, sum_time / (double)numImages);
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

int main(int argc, char *argv[])
{
    Args args;
    if (initArgs(&args, argc, argv))
        runDirectory(args);
    else
        usage();
}
