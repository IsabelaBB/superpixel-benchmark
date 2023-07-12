
#include "gft.h"

#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>
#include <strings.h>

#include <cstring>

#define IMG_COLOR 1
#define IMG_GRAYSCALE 0

typedef struct Args
{
  char *img_path, *label_path;
  char *timeFile, *dTimeFile;
  int superpixels, type;
  bool isMixed, isRoot;
  double alpha;
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
  char *spx_char, *sampling_char, *cost_char, *alpha_char, *type_char;

  args->img_path = parseArgs(argv, argc, "--img");
  args->label_path = parseArgs(argv, argc, "--label");
  args->timeFile = parseArgs(argv, argc, "--time");
  args->dTimeFile = parseArgs(argv, argc, "--dtime");

  spx_char = parseArgs(argv, argc, "--k");
  type_char = parseArgs(argv, argc, "--type");
  sampling_char = parseArgs(argv, argc, "--sampling");
  cost_char = parseArgs(argv, argc, "--cost");
  alpha_char = parseArgs(argv, argc, "--alpha");

  if (strcmp(args->timeFile, "-") == 0)
    args->timeFile = NULL;

  if (strcmp(args->dTimeFile, "-") == 0)
    args->dTimeFile = NULL;

  if (strcmp(args->label_path, "-") == 0)
    args->label_path = NULL;

  args->superpixels = strcmp(spx_char, "-") != 0 ? atoi(spx_char) : 0;
  args->type = strcmp(type_char, "-") != 0 ? atoi(type_char) : 1;
  args->isMixed = strcmp(sampling_char, "-") != 0 ? atoi(sampling_char) : 1;
  args->isRoot = strcmp(cost_char, "-") != 0 ? atoi(cost_char) : 0;
  args->alpha = strcmp(alpha_char, "-") != 0 ? atof(alpha_char) : 0.5;

  if (strcmp(args->img_path, "-") == 0 || args->superpixels < 1 || args->type < 0 || args->type > 1 || args->alpha < 0.01 || args->alpha > 0.5)
    return false;

  return true;
}

void usage()
{
  printf("Usage: superpixels args [optional args] \n");
  printf("----------------------------------\n");
  printf("args:\n");
  printf("--img   - Image (STB's supported formats) file/directory\n");
  printf("--k     - Number of superpixels\n");
  printf("optional args:\n");
  printf("--type  - image color type {0: grayscale, 1: color} -- default:1\n");
  printf("--alpha  - The relative importance between color similarity and spatial proximity. It can be in the range [0.01, 0.50]. Default:0.5\n");
  printf("--sampling - {0: grid, 1: mixed}. Default:1\n");
  printf("--cost - {0: mean, 1: root}. Default:0\n");
  printf("--label  - File/directory for final segmentation (pgm format)\n");
  printf("--time   - txt file for time log (mean time for all directory) -- only for a directory of images\n");
  printf("--dtime  - txt file for time log (time at each image) -- only for a directory of images\n");
}

gft::CImage::CImage *ReadAnyCImage(char *file)
{
  gft::CImage::CImage *cimg;
  char command[512];
  int s;

  s = strlen(file);
  if (strcasecmp(&file[s - 3], "ppm") == 0)
  {
    cimg = gft::CImage::Read(file);
  }
  else
  {
    sprintf(command, "convert %s cimage_tmp.ppm", file);
    system(command);

    cimg = gft::CImage::Read((char *)"cimage_tmp.ppm");
    system("rm cimage_tmp.ppm");
  }
  return cimg;
}

gft::Image32::Image32 *ReadAnyImage(char *file)
{
  gft::Image32::Image32 *img;
  char command[512];
  int s;

  s = strlen(file);
  if (strcasecmp(&file[s - 3], "pgm") == 0)
  {
    img = gft::Image32::Read(file);
  }
  else
  {
    sprintf(command, "convert %s image_tmp.pgm", file);
    system(command);
    img = gft::Image32::Read((char *)"image_tmp.pgm");
    system("rm image_tmp.pgm");
  }
  return img;
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
    length--;
  
  endFileName = length;
  while (imgName[length] != '/' && length > 0)
    length--;
  
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
  int color = 0xFF0000, black = 0x000000;
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
    gft::CImage::CImage *cimg = NULL, *spixels, *tmp;
    gft::Image32::Image32 *img = NULL, *label;
    double time = 0;

    if (args.type == IMG_GRAYSCALE)
    {
      img = ReadAnyImage(args.img_path);
      cimg = gft::CImage::Clone(img);
    }
    else
      cimg = ReadAnyCImage(args.img_path);

    // ************************
    if (args.type == IMG_GRAYSCALE)
    {
      start = clock();
      label = gft::Superpixels::ISF(img, args.superpixels, args.alpha, // 12.0,
                                    5.0, 2.0, 10, args.isMixed, args.isRoot);
      end = clock();
      time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
    }
    else
    {
      start = clock();
      label = gft::Superpixels::ISF(cimg, args.superpixels, args.alpha, // 12.0,
                                    5.0, 2.0, 10, args.isMixed, args.isRoot);
      end = clock();
      time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
    }
    // ************************

    char fileName[255];
    getImageName(args.img_path, fileName);

    tmp = gft::Highlight::CWideLabels(cimg, label, 2.0, &black, 0.0,
                                      true, true);
    spixels = gft::Highlight::CWideLabels(tmp, label, 1.0, &color, 0.0,
                                          true, true);
    gft::CImage::Destroy(&tmp);

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
        gft::Image32::Write(label, output_file);
      }
      else
        gft::Image32::Write(label, args.label_path);
    }

    // write dtime file
    if (args.dTimeFile != NULL)
    {
      FILE *fp = fopen(args.dTimeFile, "a+");
      fprintf(fp, "%s %.5f\n", fileName, time);
      fclose(fp);
    }

    gft::Image32::Destroy(&label);
    gft::CImage::Destroy(&spixels);
    if (cimg != NULL)
      gft::CImage::Destroy(&cimg);
    if (img != NULL)
      gft::Image32::Destroy(&img);
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
      gft::CImage::CImage *cimg = NULL, *spixels, *tmp;
      gft::Image32::Image32 *img = NULL, *label;

      if (args.type == IMG_GRAYSCALE)
      {
        img = ReadAnyImage(image_name);
        cimg = gft::CImage::Clone(img);
      }
      else
        cimg = ReadAnyCImage(image_name);
      // ********

      // run code
      if (args.type == IMG_GRAYSCALE)
      {
        start = clock();
        label = gft::Superpixels::ISF(img, args.superpixels, args.alpha, // 12.0,
                                      5.0, 2.0, 10, args.isMixed, args.isRoot);
        end = clock();
      }
      else
      {
        start = clock();
        label = gft::Superpixels::ISF(cimg, args.superpixels, args.alpha, // 12.0,
                                      5.0, 2.0, 10, args.isMixed, args.isRoot);
        end = clock();
      }
      time = ((double)(end - start)) / CLOCKS_PER_SEC;
      sum_time += time;
      // ********

      // free structures
      tmp = gft::Highlight::CWideLabels(cimg, label, 2.0, &black, 0.0,
                                        true, true);
      spixels = gft::Highlight::CWideLabels(tmp, label, 1.0, &color, 0.0,
                                            true, true);
      gft::CImage::Destroy(&tmp);

      // get image name
      char fileName[255];
      getImageName(namelist[n]->d_name, fileName);

      // write segmentation image
      if (args.label_path != NULL)
      {
        sprintf(output_file, "%s/%s.pgm", args.label_path, fileName);
        gft::Image32::Write(label, output_file);
      }

      // write dtime file
      if (args.dTimeFile != NULL)
      {
        FILE *fp = fopen(args.dTimeFile, "a+");
        fprintf(fp, "%s %.5f\n", fileName, time);
        fclose(fp);
      }

      gft::Image32::Destroy(&label);
      gft::CImage::Destroy(&spixels);
      if (cimg != NULL)
        gft::CImage::Destroy(&cimg);
      if (img != NULL)
        gft::Image32::Destroy(&img);
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

int main(int argc, char **argv)
{

  Args args;
  if (initArgs(&args, argc, argv))
    runDirectory(args);
  else
    usage();
  //-------------------------------------------------

  return 0;
}
