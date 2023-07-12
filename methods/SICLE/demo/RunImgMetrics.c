/*****************************************************************************\
* RunImgMetrics.c
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-06-15
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#include "ift.h"
#include "iftArgs.h"
#include "iftMetrics.h"

void usage();
void readImgInputs
(const iftArgs *args, iftImage **label_img,  iftImage **img, iftImage **gt_img);

int main(int argc, char const *argv[])
{
  //-------------------------------------------------------------------------//
  bool has_req, has_help;
  iftArgs *args;

  args = iftCreateArgs(argc, argv);

  has_req = iftExistArg(args, "labels") && iftExistArg(args, "gt");
  has_help = iftExistArg(args, "help");

  if(has_req == false || has_help == true)
  {
    usage(); 
    iftDestroyArgs(&args);
    return EXIT_FAILURE;
  }
  //-------------------------------------------------------------------------//
  bool as_csv;
  int min_label, max_label, num_labels;
  char out_str[IFT_STR_DEFAULT_SIZE], tmp_str[IFT_STR_DEFAULT_SIZE];
  iftImage *label_img, *img, *gt_img;

  readImgInputs(args, &label_img, &img, &gt_img);

  iftMinMaxValues(label_img, &min_label, &max_label);
  num_labels = max_label - min_label + 1;

  as_csv = iftExistArg(args, "csv");
  if(as_csv == false) 
    sprintf(tmp_str, "Number of Superspels: %d\n", num_labels);
  else sprintf(tmp_str, "%d", num_labels);
  strcat(out_str, tmp_str);

  if(gt_img != NULL)
  {
    float br, ue;

    br = iftBoundRecall(label_img, gt_img);
    ue = iftUnderSegmError(label_img, gt_img);

    if(as_csv == false) 
    {
      sprintf(tmp_str, "BR: %.3f\n", br);
      strcat(out_str, tmp_str);
      sprintf(tmp_str, "UE: %.3f\n", ue);
      strcat(out_str, tmp_str);
    }
    else 
    {
      sprintf(tmp_str, ",%f", br);
      strcat(out_str, tmp_str);
      sprintf(tmp_str, ",%f", ue);
      strcat(out_str, tmp_str);
    }
  }
  
  iftDestroyArgs(&args);
  iftDestroyImage(&label_img);
  if(gt_img != NULL) iftDestroyImage(&gt_img);
  if(img != NULL) iftDestroyImage(&img);

  puts(out_str);
  return EXIT_SUCCESS;
}

void usage()
{
  const int SKIP_IND = 15; // For indentation purposes
  printf("\nThe required parameters are:\n");
  printf("%-*s %s\n", SKIP_IND, "--labels", 
         "Input label image");
  printf("%-*s %s\n", SKIP_IND, "--gt", 
         "Input ground-truth image.");

  printf("\nThe optional parameters are:\n");
  printf("%-*s %s\n", SKIP_IND, "--csv", 
         "Flag for printing the output as comma-separated values (CSV).");
  printf("%-*s %s\n", SKIP_IND, "--help", 
         "Prints this message");

  printf("\n");
}

void readImgInputs
(const iftArgs *args, iftImage **label_img,  iftImage **img, iftImage **gt_img)
{
  #if IFT_DEBUG //-----------------------------------------------------------//
  assert(args != NULL);
  assert(img != NULL);
  assert(label_img != NULL);
  assert(gt_img != NULL);
  #endif //------------------------------------------------------------------//
  const char *PATH;

  if(iftHasArgVal(args, "labels") == true) 
  {
    PATH = iftGetArg(args, "labels");

    (*label_img) = iftReadImageByExt(PATH);
  }
  else iftError("No label image path was given", "readImgInputs");

  
  if(iftHasArgVal(args, "gt") == true) 
  {
    PATH = iftGetArg(args, "gt");

    (*gt_img) = iftReadImageByExt(PATH);
  }
  else iftError("No ground-truth image path was given", "readImgInputs");

}