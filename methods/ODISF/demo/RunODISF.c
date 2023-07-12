/*****************************************************************************\
* RunODISF.c
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-05-15
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#include "ift.h"
#include "iftArgs.h"
#include "iftODISF.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>

void usage();
void readImgInputs(const iftArgs *args, iftImage **img, iftImage **mask, iftImage **objsm,
                   const char **path);
void setODISFParams(iftODISF **odisf, const iftArgs *args);
void setODISFSampl(iftODISF **odisf, const iftArgs *args);

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

void runDirectory(iftArgs *args)
{
  // determine mode : file or path
  struct stat sb;
  clock_t start, end;

  const char *path;
  if (iftHasArgVal(args, "img") == true)
    path = iftGetArg(args, "img");

  if (stat(path, &sb) == -1)
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
    const char *LABEL_PATH = NULL;
    iftImage *img, *mask, *objsm, *labels;
    iftODISF *odisf;
    double time = 0;

    readImgInputs(args, &img, &mask, &objsm, &LABEL_PATH);
    odisf = iftCreateODISF(img, mask, objsm);
    iftDestroyImage(&img);
    if (mask != NULL)
      iftDestroyImage(&mask);
    if (objsm != NULL)
      iftDestroyImage(&objsm);
    setODISFParams(&odisf, args);

    // ************************
    start = clock();
    setODISFSampl(&odisf, args);
    labels = iftRunODISF(odisf);
    end = clock();
    time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
    // ************************

    char fileName[255];
    getImageName(path, fileName);

    iftDestroyODISF(&odisf);

    if (iftExistArg(args, "out") == true)
    {
      struct stat label_stats; // determine mode : file or path
      if (stat(LABEL_PATH, &label_stats) == -1)
      {
        perror("stat");
        exit(EXIT_SUCCESS);
      }

      if ((label_stats.st_mode & S_IFMT) == S_IFDIR)
      {
        char output_file[255];
        sprintf(output_file, "%s/%s.pgm", LABEL_PATH, fileName);
        iftWriteImageByExt(labels, output_file);
      }
      else
        iftWriteImageByExt(labels, LABEL_PATH);
    }
    iftDestroyImage(&labels);

    // write dtime file
    if (iftExistArg(args, "dtime") == true)
    {
      const char *dfile_time = iftGetArg(args, "dtime");
      FILE *fp = fopen(dfile_time, "a+");
      fprintf(fp, "%s %.5f\n", fileName, time);
      fclose(fp);
    }
  }
  else if (type == 0)
  {
    // get file list
    struct dirent **namelist;
    int n = scandir(path, &namelist, &filter, alphasort);
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

    int numImages = n;
    double sum_time = 0;
    const char *out_dir = NULL;
    const char *mask_dir;
    const char *objsm_dir;
    const char *dfile_time;
    int n0, nf;
    const char *sampling;

    if (iftExistArg(args, "out") == true)
      out_dir = iftGetArg(args, "out");
    else
      out_dir = NULL;

    if (iftExistArg(args, "objsm") == true)
      objsm_dir = iftGetArg(args, "objsm");
    else
      objsm_dir = NULL;

    if (iftExistArg(args, "sampl-op") == true)
      sampling = iftGetArg(args, "sampl-op");

    if (iftExistArg(args, "mask") == true)
      mask_dir = iftGetArg(args, "mask");
    else
      mask_dir = NULL;

    if (iftExistArg(args, "dtime") == true)
      dfile_time = iftGetArg(args, "dtime");
    else
      dfile_time = NULL;

    if (iftExistArg(args, "n0") == true)
      n0 = atoi(iftGetArg(args, "n0"));

    if (iftExistArg(args, "nf") == true)
      nf = atoi(iftGetArg(args, "nf"));

    while (n--)
    {
      // get image name
      sprintf(image_name, "%s/%s", path, namelist[n]->d_name);

      // alloc structures
      iftArgs *args_local;
      const char *LABEL_PATH;
      iftImage *img, *mask, *objsm, *labels;
      iftODISF *odisf;
      double time = 0;

      // get image name
      char fileName[255];
      getImageName(namelist[n]->d_name, fileName);

      // get args
      char **argv = (char **)malloc(14 * sizeof(char *));
      int idx = 0;

      argv[idx] = (char *)malloc(255 * sizeof(char));
      strcpy(argv[idx], "--img");
      idx++;
      argv[idx] = (char *)malloc(255 * sizeof(char));
      sprintf(argv[idx], "%s", image_name);

      if (out_dir != NULL)
      {
        idx++;
        argv[idx] = (char *)malloc(255 * sizeof(char));
        strcpy(argv[idx], "--out");
        idx++;
        argv[idx] = (char *)malloc(255 * sizeof(char));
        sprintf(argv[idx], "%s/%s.pgm", out_dir, fileName);
      }

      if (objsm_dir != NULL)
      {
        idx++;
        argv[idx] = (char *)malloc(255 * sizeof(char));
        strcpy(argv[idx], "--objsm");
        idx++;
        argv[idx] = (char *)malloc(255 * sizeof(char));
        sprintf(argv[idx], "%s/%s.png", objsm_dir, fileName);
      }

      if (mask_dir != NULL)
      {
        idx++;
        argv[idx] = (char *)malloc(255 * sizeof(char));
        strcpy(argv[idx], "--mask");
        idx++;
        argv[idx] = (char *)malloc(255 * sizeof(char));
        sprintf(argv[idx], "%s/%s.png", mask_dir, fileName);
      }

      idx++;
      argv[idx] = (char *)malloc(255 * sizeof(char));
      strcpy(argv[idx], "--sampl-op");
      idx++;
      argv[idx] = (char *)malloc(255 * sizeof(char));
      sprintf(argv[idx], "%s", sampling);

      idx++;
      argv[idx] = (char *)malloc(255 * sizeof(char));
      strcpy(argv[idx], "--n0");
      idx++;
      argv[idx] = (char *)malloc(255 * sizeof(char));
      sprintf(argv[idx], "%d", n0);

      idx++;
      argv[idx] = (char *)malloc(255 * sizeof(char));
      strcpy(argv[idx], "--nf");
      idx++;
      argv[idx] = (char *)malloc(255 * sizeof(char));
      sprintf(argv[idx], "%d", nf);
      idx++;

      const int argc = idx;
      args_local = iftCreateArgs(argc, argv);

      // init structures
      readImgInputs(args_local, &img, &mask, &objsm, &LABEL_PATH);
      odisf = iftCreateODISF(img, mask, objsm);
      iftDestroyImage(&img);
      if (mask != NULL)
        iftDestroyImage(&mask);
      if (objsm != NULL)
        iftDestroyImage(&objsm);
      setODISFParams(&odisf, args_local);
      // ********

      // run code
      start = clock();
      setODISFSampl(&odisf, args_local);
      labels = iftRunODISF(odisf);
      end = clock();
      time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
      sum_time += time;
      // ********

      // write segmentation image
      if (out_dir != NULL)
        iftWriteImageByExt(labels, LABEL_PATH);

      // write dtime file
      if (dfile_time != NULL)
      {
        FILE *fp = fopen(dfile_time, "a+");
        fprintf(fp, "%s %.5f\n", fileName, time);
        fclose(fp);
      }

      // free structures
      iftDestroyODISF(&odisf);
      iftDestroyImage(&labels);
      iftDestroyArgs(&args_local);
      for (int i = 0; i < argc; i++)
        free(argv[i]);
      free(argv);
      // ********

      free(namelist[n]);
    }

    if (iftExistArg(args, "time") == true)
    {
      char *timeFile = iftGetArg(args, "time");
      int superpixels = atoi(iftGetArg(args, "nf"));
      FILE *fp = fopen(timeFile, "a+");

      fprintf(fp, "%d %.6f\n", superpixels, sum_time / (double)numImages);
      fclose(fp);
    }
    iftDestroyArgs(&args);

    free(image_name);
    free(namelist);
    // ************************
  }
}

int main(int argc, char const *argv[])
{
  //-------------------------------------------------------------------------//
  bool has_req, has_help;
  iftArgs *args;

  args = iftCreateArgs(argc, argv);

  has_req = iftExistArg(args, "img");
  has_help = iftExistArg(args, "help");

  if (has_req == false || has_help == true)
  {
    usage();
    iftDestroyArgs(&args);
    return EXIT_FAILURE;
  }
  //-------------------------------------------------------------------------//

  runDirectory(args);
  iftDestroyArgs(&args);

  return EXIT_SUCCESS;
}

void usage()
{
  const int SKIP_IND = 15; // For indentation purposes
  printf("\nThe required parameters are:\n");
  printf("%-*s %s\n", SKIP_IND, "--img",
         "Input 2D image");

  printf("\nThe optional parameters are:\n");
  printf("%-*s %s\n", SKIP_IND, "--out",
         "Output 2D label image");
  printf("%-*s %s\n", SKIP_IND, "--no-diag-adj",
         "Disable search scope to consider 8-adjacency.");
  printf("%-*s %s\n", SKIP_IND, "--mask",
         "Mask image indicating the region of interest.");
  printf("%-*s %s\n", SKIP_IND, "--n0",
         "Desired initial number of seeds. Default: 8000");
  printf("%-*s %s\n", SKIP_IND, "--nf",
         "Desired final number of superpixels. Default: 200");
  printf("%-*s %s\n", SKIP_IND, "--sampl-op",
         "Seed sampling algorithm. Options: "
         "grid, rnd. Default: grid");
  printf("%-*s %s\n", SKIP_IND, "--objsm",
         "Grayscale object saliency map.");
  printf("%-*s %s\n", SKIP_IND, "--dtime",
         "txt file for time log (mean time for all directory) -- only for a directory of images");
  printf("%-*s %s\n", SKIP_IND, "--time",
         "txt file for time log (time at each image) -- only for a directory of images.");
  printf("%-*s %s\n", SKIP_IND, "--help",
         "Prints this message");
  printf("\n");
}

void readImgInputs(const iftArgs *args, iftImage **img, iftImage **mask, iftImage **objsm,
                   const char **path)
{
#if IFT_DEBUG //-----------------------------------------------------------//
  assert(args != NULL);
  assert(img != NULL);
  assert(mask != NULL);
  assert(objsm != NULL);
  assert(path != NULL);
#endif //------------------------------------------------------------------//
  const char *PATH;

  if (iftHasArgVal(args, "img") == true)
    (*img) = iftReadImageByExt(iftGetArg(args, "img"));
  else
    iftError("No image path was given", "readImgInputs");

  if (iftHasArgVal(args, "out") == true)
  {
    PATH = iftGetArg(args, "out");

    if (iftIsImagePathnameValid(PATH) == true)
      (*path) = PATH;
    else
      iftError("Unknown image type", "readImgInputs");
  }
  // else iftError("No output label path was given", "readImgInputs");

  if (iftExistArg(args, "mask") == true)
  {
    if (iftHasArgVal(args, "mask") == true)
    {
      (*mask) = iftReadImageByExt(iftGetArg(args, "mask"));

      iftVerifyImageDomains(*img, *mask, "readImgInputs");
    }
    else
      iftError("No mask path was given", "readImgInputs");
  }
  else
    (*mask) = NULL;

  if (iftExistArg(args, "objsm") == true)
  {
    if (iftHasArgVal(args, "objsm") == true)
    {
      (*objsm) = iftReadImageByExt(iftGetArg(args, "objsm"));

      iftVerifyImageDomains(*img, *objsm, "readImgInputs");
    }
    else
      iftError("No object saliency map path was given", "readImgInputs");
  }
  else
    (*objsm) = NULL;
}

void setODISFParams(iftODISF **odisf, const iftArgs *args)
{
#if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(*odisf != NULL);
  assert(args != NULL);
#endif //------------------------------------------------------------------//
  int n0, nf;

  iftODISFUseDiagAdj(odisf, !iftExistArg(args, "no-diag-adj"));

  if (iftExistArg(args, "n0") == true)
  {
    if (iftHasArgVal(args, "n0") == true)
    {
      n0 = atoi(iftGetArg(args, "n0"));
      iftODISFSetN0(odisf, n0);
    }
    else
      iftError("No initial number of seeds was given", "setODISFParams");
  }

  if (iftExistArg(args, "nf") == true)
  {
    if (iftHasArgVal(args, "nf") == true)
    {
      nf = atoi(iftGetArg(args, "nf"));
      iftODISFSetNf(odisf, nf);
    }
    else
      iftError("No desired quantity of superpixels was given",
               "setODISFParams");
  }
}

void setODISFSampl(iftODISF **odisf, const iftArgs *args)
{
#if IFT_DEBUG //-----------------------------------------------------------//
  assert(odisf != NULL);
  assert(*odisf != NULL);
  assert(args != NULL);
#endif //------------------------------------------------------------------//
  if (iftExistArg(args, "sampl-op") == true)
  {
    if (iftHasArgVal(args, "sampl-op") == true)
    {
      const char *OPT;

      OPT = iftGetArg(args, "sampl-op");

      if (iftCompareStrings(OPT, "grid"))
        iftODISFSetSamplOpt(odisf, IFT_ODISF_SAMPL_GRID);
      else if (iftCompareStrings(OPT, "rnd"))
        iftODISFSetSamplOpt(odisf, IFT_ODISF_SAMPL_RND);
      else
        iftError("Unknown seed sampling algorithm", "setODISFSampl");
    }
    else
      iftError("No sampling algorithm was given", "setODISFSampl");
  }
}