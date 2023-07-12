/**
 * @file
 * @brief Overlay borders from one or more segmentations.
 * @author Felipe Lemes Galvao
 * @date November, 2017
 */

#include <ift.h>
#include "iftExperimentUtility.h"

int main(int argc, char* argv[])
{
  if (argc < 4) 
    iftError("Usage: iftOverlayBorders <image> <output> [seg1, seg2...]", "main");

  // Load Image
  iftImage *img = iftReadImageByExt(argv[1]);

  // Prepare border colors
  iftColorTable *colorTb = iftCreateColorTable(argc - 3);

  // Overlay borders, one segmentation at a time
  for (int i = 3; i < argc; ++i) {
    // Load label map
    iftImage *superpixelLabelMap = iftReadImageByExt(argv[i]);

    // Actual overlay
    iftImage *tmpImg = iftOverlaySegmentationBorders(img, superpixelLabelMap, colorTb->color[i-3]);
    iftSwap(img, tmpImg);

    // Clean up
    iftDestroyImage(&superpixelLabelMap);
    iftDestroyImage(&tmpImg);
  }

  iftWriteImageByExt(img, "%s", argv[2]);
  iftDestroyImage(&img);

  return 0;
}
