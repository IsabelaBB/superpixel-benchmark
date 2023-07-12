/**
 * @file
 * @brief Relabels each connected component to guarantee a consistent map.
 * @author Felipe Lemes Galvao
 */

#include <ift.h>
#include "iftSuperpixelFeatures.h"
#include "iftExperimentUtility.h"

int main(int argc, char* argv[]) {
  if (argc < 3)
    iftError("Usage: %s <input img> <output_img> [min_region_size]", "main", argv[0]);

  iftImage *labelMap = iftReadImageByExt(argv[1]);
  int minSize = (argc > 3) ? atol(argv[3]) : 0;
  iftImage *res = iftForceLabelMapConnectivity(labelMap, minSize);
  iftWriteImageByExt(res, argv[2]);

  iftDestroyImage(&labelMap);
  iftDestroyImage(&res);

  return 0;
}
