/*****************************************************************************\
* RunSICLE.c
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-07-08
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#include "ift.h"
#include "iftArgs.h"
#include "iftSICLE.h"

#include <sys/stat.h>
#include <dirent.h>
#include <time.h>
#include <string.h>

void usage();
void readImgInputs
(const iftArgs *args, iftImage **img, iftImage **mask, iftImage **objsm,
 const char **path);
void setSICLEParams
(iftSICLE **sicle, const iftArgs *args);
void setSICLESampl
(iftSICLE **sicle, const iftArgs *args);
void setSICLEArcCost
(iftSICLE **sicle, const iftArgs *args);
void setSICLERem
(iftSICLE **sicle, const iftArgs *args);
void writeLabels
(const iftImage **labels, const char *path);

// from: https://stackoverflow.com/questions/2736753/how-to-remove-extension-from-file-name
char *remove_ext (const char* myStr, char extSep, char pathSep) {
    char *retStr, *lastExt, *lastPath;

    // Error checks and allocate string.

    if (myStr == NULL) return NULL;
    if ((retStr = malloc (strlen (myStr) + 1)) == NULL) return NULL;

    // Make a copy and find the relevant characters.

    strcpy (retStr, myStr);
    lastExt = strrchr (retStr, extSep);
    lastPath = (pathSep == 0) ? NULL : strrchr (retStr, pathSep);

    // If it has an extension separator.

    if (lastExt != NULL) {
        // and it's to the right of the path separator.

        if (lastPath != NULL) {
            if (lastPath < lastExt) {
                // then remove it.

                *lastExt = '\0';
            }
        } else {
            // Has extension separator with no path separator.

            *lastExt = '\0';
        }
    }

    // Return the modified string.

    return retStr;
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
  for(int i=0; i < extSize; i++)
    ext[i] = name->d_name[pos + i];
  
  if (extSize > 0 && (extSize == 3 && (strcmp(ext, "jpg") || strcmp(ext, "png") || strcmp(ext, "ppm") || strcmp(ext, "pgm") || strcmp(ext, "bmp"))) ||
      (extSize == 4 && (strcmp(ext, "tiff") || strcmp(ext, "jpeg"))))
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

void copyArg(char** argv, int *idx, char* key, char *value){
  argv[*idx] = (char *)malloc(255 * sizeof(char));
  strcpy(argv[*idx], key);
  (*idx)++;
  argv[*idx] = (char *)malloc(255 * sizeof(char));
  sprintf(argv[*idx], "%s", value);
  (*idx)++;
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
    const char *LABEL_PATH = NULL;
    iftImage *img, *mask, *objsm, **labels;
    iftSICLE *sicle;
    double time = 0;

    readImgInputs(args, &img, &mask, &objsm, &LABEL_PATH);
    sicle = iftCreateSICLE(img, mask, objsm);
    iftDestroyImage(&img);
    if (mask != NULL) iftDestroyImage(&mask);
    if (objsm != NULL) iftDestroyImage(&objsm);

    setSICLEParams(&sicle, args);
    setSICLESampl(&sicle, args);
    setSICLEArcCost(&sicle, args);
    setSICLERem(&sicle, args);
    
    // ************************
    start = clock();
    labels = iftRunSICLE(sicle);
    end = clock();
    time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
    // ************************

    char fileName[255];
    getImageName((char*)path, fileName);

    iftDestroySICLE(&sicle);

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
        iftWriteImageByExt(labels[0], output_file);
      }
      else
        iftWriteImageByExt(labels[0], LABEL_PATH);
    }
    iftDestroyImage(&(labels[0]));
    free(labels);

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

    //printf(" %i image(s) found\n", n);
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

    if (iftExistArg(args, "out") == true)
      out_dir = iftGetArg(args, "out");
    else
      out_dir = NULL;

    if (iftExistArg(args, "objsm") == true)
      objsm_dir = iftGetArg(args, "objsm");
    else
      objsm_dir = NULL;

    if (iftExistArg(args, "mask") == true)
      mask_dir = iftGetArg(args, "mask");
    else
      mask_dir = NULL;

    if (iftExistArg(args, "dtime") == true)
      dfile_time = iftGetArg(args, "dtime");
    else
      dfile_time = NULL;

    while (n--)
    {
      // get image name
      sprintf(image_name, "%s/%s", path, namelist[n]->d_name);

      // alloc structures
      iftArgs *args_local;
      const char *LABEL_PATH;
      const iftImage *img, *mask, *objsm, **labels;
      iftSICLE *sicle;
      double time = 0;

      // get image name
      char fileName[255];

      getImageName(namelist[n]->d_name, fileName);
      
      // get args
      char **argv = (char **)malloc(24 * sizeof(char *));
      int idx = 0;

      copyArg(argv, &idx, "--img", image_name);

      if (out_dir != NULL)
      {
        char value[255];
        sprintf(value, "%s/%s.pgm", out_dir, fileName);
        copyArg(argv, &idx, "--out", value);
      }

      if (objsm_dir != NULL)
      {
        char value[255];
        sprintf(value, "%s/%s.png", objsm_dir, fileName);
        copyArg(argv, &idx, "--objsm", value);
      }

      if (mask_dir != NULL)
      {
        char value[255];
        sprintf(value, "%s/%s.png", mask_dir, fileName);
        copyArg(argv, &idx, "--mask", value);
      }

      if (iftExistArg(args, "no-diag-adj")){
          char *value = (char *)iftGetArg(args, "no-diag-adj");
          copyArg(argv, &idx, "--no-diag-adj", value);
      }

      if (iftExistArg(args, "boost")){
          char *value = (char *)iftGetArg(args, "boost");
          copyArg(argv, &idx, "--boost", value);
      }

      if (iftExistArg(args, "max-iters")){
          char *value = (char *)iftGetArg(args, "max-iters");
          copyArg(argv, &idx, "--max-iters", value);
      }

      if (iftExistArg(args, "n0")){
          char *value = (char *)iftGetArg(args, "n0");
          copyArg(argv, &idx, "--n0", value);
      }

      if (iftExistArg(args, "nf")){
          char *value = (char *)iftGetArg(args, "nf");
          copyArg(argv, &idx, "--nf", value);
      }

      if (iftExistArg(args, "sampl-op")){
          char *value = (char *)iftGetArg(args, "sampl-op");
          copyArg(argv, &idx, "--sampl-op", value);
      }

      if (iftExistArg(args, "arc-op")){
          char *value = (char *)iftGetArg(args, "arc-op");
          copyArg(argv, &idx, "--arc-op", value);
      }

      if (iftExistArg(args, "rem-op")){
          char *value = (char *)iftGetArg(args, "rem-op");
          copyArg(argv, &idx, "--rem-op", value);
      }

      const int argc = idx;
      args_local = iftCreateArgs(argc, argv);

      readImgInputs(args_local, &img, &mask, &objsm, &LABEL_PATH);

      sicle = iftCreateSICLE(img, mask, objsm);

      setSICLEParams(&sicle, args_local);
      setSICLESampl(&sicle, args_local);
      setSICLEArcCost(&sicle, args_local);
      setSICLERem(&sicle, args_local);

      // ********

      // run code
      start = clock();
      labels = iftRunSICLE(sicle);
      end = clock();
      time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
      sum_time += time;
      // ********

      iftDestroyImage(&img);
      if (mask != NULL)
        iftDestroyImage(&mask);
      if (objsm != NULL)
        iftDestroyImage(&objsm);

      // write segmentation image
      if (out_dir != NULL)
        iftWriteImageByExt(labels[0], LABEL_PATH);

      // write dtime file
      if (dfile_time != NULL)
      {
        FILE *fp = fopen(dfile_time, "a+");
        fprintf(fp, "%s %.5f\n", fileName, time);
        fclose(fp);
      }

      // free structures
      iftDestroySICLE(&sicle);
      iftDestroyImage(&(labels[0]));
      free(labels);
      iftDestroyArgs(&args_local);
      for (int i = 0; i < argc; i++)
        free(argv[i]);
      free(argv);
      // ********

      free(namelist[n]);
    }

    if (iftExistArg(args, "time") == true)
    {
      const char *timeFile = iftGetArg(args, "time");
      int superpixels = atoi(iftGetArg(args, "nf"));
      FILE *fp = fopen(timeFile, "a+");

      fprintf(fp, "%d %.5f\n", superpixels, sum_time / (double)numImages);
      fclose(fp);
    }

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

  has_req = iftExistArg(args, "img")  && iftExistArg(args, "out");
  has_help = iftExistArg(args, "help");

  if(has_req == false || has_help == true)
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
         "Input image");
  printf("%-*s %s\n", SKIP_IND, "--out",
         "Output label image");

  printf("\nThe optional parameters are:\n");
  printf("%-*s %s\n", SKIP_IND, "--no-diag-adj",
         "Disable search scope to consider 8-adjacency.");
  printf("%-*s %s\n", SKIP_IND, "--mask",
         "Mask image indicating the region of interest.");
  printf("%-*s %s\n", SKIP_IND, "--max-iters",
         "Maximum number of iterations for segmentation. Default: 5");
  printf("%-*s %s\n", SKIP_IND, "--n0",
         "Desired initial number of seeds. Default: 3000");
  printf("%-*s %s\n", SKIP_IND, "--nf",
         "Desired final number of superpixels. Default: 200");
  printf("%-*s %s\n", SKIP_IND, "--scales",
         "Comma-separated list of superpixel scales. Default: None");
 printf("%-*s %s\n", SKIP_IND, "--boost",
        "Enables boosting performance for good saliency maps.");
  printf("%-*s %s\n", SKIP_IND, "--objsm",
         "Grayscale object saliency map.");
  printf("%-*s %s\n", SKIP_IND, "--help",
         "Prints this message");

  printf("\nThe SICLE configuration options are:\n");
  printf("%-*s %s\n", SKIP_IND, "--arc-op",
         "IFT arc cost function. Options: "
         "root, dyn. Default: root");
  printf("%-*s %s\n", SKIP_IND, "--sampl-op",
         "Seed sampling algorithm. Options: "
         "grid, rnd. Default: rnd");
  printf("%-*s %s\n", SKIP_IND, "--rem-op",
         "Seed removal criterion. Options: "
         "max-contr, min-contr, size, rnd, "
         "max-sc, min-sc. Default: min-sc");

  printf("\n");
}

void readImgInputs
(const iftArgs *args, iftImage **img, iftImage **mask, iftImage **objsm,
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

  if(iftHasArgVal(args, "img") == true)
  {
    PATH = iftGetArg(args, "img");

    (*img) = iftReadImageByExt(PATH);
  }
  else iftError("No image path was given", "readImgInputs");
  //else img = NULL;

  if(iftHasArgVal(args, "out") == true)
  {
    PATH = iftGetArg(args, "out");

    if(iftIsImagePathnameValid(PATH) == true) (*path) = PATH;
    else iftError("Unknown image type", "readImgInputs");
  }
  //else iftError("No output label path was given", "readImgInputs");
  
  if(iftExistArg(args, "mask") == true)
  {
    if(iftHasArgVal(args, "mask") == true)
    {
      PATH = iftGetArg(args, "mask");

      (*mask) = iftReadImageByExt(PATH);

      iftVerifyImageDomains(*img, *mask, "readImgInputs");
    }
    else iftError("No mask path was given", "readImgInputs");
  }
  else (*mask) = NULL;

  if(iftExistArg(args, "objsm") == true)
  {
    if(iftHasArgVal(args, "objsm") == true)
    {
      PATH = iftGetArg(args, "objsm");

      (*objsm) = iftReadImageByExt(PATH);

      iftVerifyImageDomains(*img, *objsm, "readImgInputs");
    }
    else iftError("No object saliency map path was given", "readImgInputs");
  }
  else (*objsm) = NULL;
}

void setSICLEParams
(iftSICLE **sicle, const iftArgs *args)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  assert(args != NULL);
  #endif //------------------------------------------------------------------//
  int n0, nf, max_iters;

  iftSICLEUseDiagAdj(sicle, !iftExistArg(args, "no-diag-adj"));

  iftSICLEEnableBoost(sicle, iftExistArg(args, "boost"));

  if(iftExistArg(args, "max-iters") == true)
  {
    if(iftHasArgVal(args, "max-iters") == true)
    {
      max_iters = atoi(iftGetArg(args, "max-iters"));
      if(max_iters <= 1)
        iftError("The maximum number of iterations must be greater than 1",
                 "setSICLEParams");
      iftSICLESetMaxIters(sicle, max_iters);
    }
    else iftError("No maximum number of iterations was given",
                  "setSICLEParams");
  }

  if(iftExistArg(args, "n0") == true)
  {
    if(iftHasArgVal(args, "n0") == true)
    {
      n0 = atoi(iftGetArg(args, "n0"));
      iftSICLESetN0(sicle, n0);
    }
    else iftError("No initial number of seeds was given", "setSICLEParams");
  }

  if(iftExistArg(args, "nf") == true)
  {
    if(iftHasArgVal(args, "nf") == true)
    {
      nf = atoi(iftGetArg(args, "nf"));

	  if(nf < iftSICLEGetN0(*sicle)) iftSICLESetNf(sicle, nf);
      else iftError("The number of superpixels must be greater than N0",
                    "setSICLEParams");
    }
    else iftError("No desired quantity of superpixels was given",
                  "setSICLEParams");
  }

  if(iftExistArg(args, "scales") == true)
  {
    if(iftHasArgVal(args, "scales") == true)
    {
      const char *VAL;
      int i, prev;
      char *tmp, *tok;
	  int num_scales;
	  int *scales;
	  iftSet *list_scales;

	  list_scales = NULL;

      VAL = iftGetArg(args, "scales");
      tmp = iftCopyString(VAL);
      tok = strtok(tmp, ",");

      i = 0;
	  prev = iftSICLEGetN0(*sicle);
      while(tok != NULL)
      {
        int curr;

        curr = atoi(tok);
        if(curr < prev && curr > iftSICLEGetNf(*sicle))
        {
		  iftInsertSet(&list_scales,curr);
		  prev = curr;
        }
        else iftError("The order of scales must be statically decreasing and within ]N0,Nf[",
        		   	  "setSICLEParams");

        tok = strtok(NULL, ",");
        i++;
      }
      free(tmp);

      if(i == 0) iftError("No list of scales was given", "setSICLEParams");

	  num_scales = i;
	  scales = calloc(num_scales, sizeof(int));
	  while(i > 0)
		scales[--i] = iftRemoveSet(&list_scales);

      iftSICLESetScales(sicle, num_scales, scales);
	  free(scales);
    }
    else iftError("No list of superpixel scales was given",
                  "setSICLEParams");
  }


}

void setSICLESampl
(iftSICLE **sicle, const iftArgs *args)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  assert(args != NULL);
  #endif //------------------------------------------------------------------//
  if(iftExistArg(args, "sampl-op") == true)
  {
    if(iftHasArgVal(args, "sampl-op") == true)
    {
      const char *OPT;

      OPT = iftGetArg(args, "sampl-op");

      if(iftCompareStrings(OPT, "grid"))
        iftSICLESetSamplOpt(sicle, IFT_SICLE_SAMPL_GRID);
      else if(iftCompareStrings(OPT, "rnd"))
        iftSICLESetSamplOpt(sicle, IFT_SICLE_SAMPL_RND);
      else iftError("Unknown seed sampling algorithm", "setSICLESampl");
    }
    else iftError("No sampling algorithm was given", "setSICLESampl");
  }
}

void setSICLEArcCost
(iftSICLE **sicle, const iftArgs *args)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  assert(args != NULL);
  #endif //------------------------------------------------------------------//
  if(iftExistArg(args, "arc-op") == true)
  {
    if(iftHasArgVal(args, "arc-op") == true)
    {
      const char *OPT;

      OPT = iftGetArg(args, "arc-op");

      if(iftCompareStrings(OPT, "root"))
        iftSICLESetArcCostOpt(sicle, IFT_SICLE_ARCCOST_ROOT);
      else if(iftCompareStrings(OPT, "dyn"))
        iftSICLESetArcCostOpt(sicle, IFT_SICLE_ARCCOST_DYN);
      else iftError("Unknown arc-cost function", "setSICLEArcCost");
    }
    else iftError("No arc-cost function was given", "setSICLEArcCost");
  }
}

void setSICLERem
(iftSICLE **sicle, const iftArgs *args)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(sicle != NULL);
  assert(*sicle != NULL);
  assert(args != NULL);
  #endif //------------------------------------------------------------------//
  if(iftExistArg(args, "rem-op") == true)
  {
    if(iftHasArgVal(args, "rem-op") == true)
    {
      const char *OPT;

      OPT = iftGetArg(args, "rem-op");

      if(iftCompareStrings(OPT, "min-contr"))
        iftSICLESetRemOpt(sicle, IFT_SICLE_REM_MINCONTR);
      else if(iftCompareStrings(OPT, "max-contr"))
        iftSICLESetRemOpt(sicle, IFT_SICLE_REM_MAXCONTR);
      else if(iftCompareStrings(OPT, "max-sc"))
        iftSICLESetRemOpt(sicle, IFT_SICLE_REM_MAXSC);
      else if(iftCompareStrings(OPT, "min-sc"))
        iftSICLESetRemOpt(sicle, IFT_SICLE_REM_MINSC);
      else if(iftCompareStrings(OPT, "size"))
        iftSICLESetRemOpt(sicle, IFT_SICLE_REM_SIZE);
      else if(iftCompareStrings(OPT, "rnd"))
        iftSICLESetRemOpt(sicle, IFT_SICLE_REM_RND);
      else iftError("Unknown seed removal criterion", "setSICLERem");
    }
    else iftError("No seed removal criterion was given", "setSICLERem");
  }
}
