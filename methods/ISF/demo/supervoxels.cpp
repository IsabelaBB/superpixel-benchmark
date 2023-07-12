
#include "gft.h"

#include <string>

int main(int argc, char **argv){
  gft::Scene32::Scene32 *scn, *label;
  char filename[512];
  int k;
  double alpha;
  bool isMixed, isRoot;

  if(argc != 7){
    fprintf(stdout,"usage:\n");
    fprintf(stdout,"supervoxels <image> <k> <alpha> <sampling> <cost-function> <output>\n");
    fprintf(stdout,"\t k........ The desired number of supervoxels\n");
    fprintf(stdout,"\t alpha.... The relative importance between color similarity and spatial proximity.\n");
    fprintf(stdout,"\t           It can be in the range [0.01, 0.50]\n");
    fprintf(stdout,"\t sampling. 0: grid\n");
    fprintf(stdout,"\t           1: mixed\n");
    fprintf(stdout,"\t cost..... 0: mean\n");
    fprintf(stdout,"\t           1: root\n");
    fprintf(stdout,"\t output... The target path of the output sans extension.\n");
    exit(0);
  }

  strcpy(filename, argv[1]);
  k = atoi(argv[2]);
  alpha = atof(argv[3]);
  isMixed = atoi(argv[4]);
  isRoot = atoi(argv[5]);

  scn = gft::Scene32::Read(filename);

  //-------------------------------------------------
  // IFT-SLIC:
  float real_proc_time = 0.0;
  label = gft::Superpixels::ISF(scn, k, alpha, 5.0, 2.0, 10, isMixed, isRoot/*, &real_proc_time*/);

  sprintf(filename, "%s.scn", argv[6]);
  gft::Scene32::Write(label, filename);

  gft::Scene32::Destroy(&label);
  gft::Scene32::Destroy(&scn);
  return 0;
}



