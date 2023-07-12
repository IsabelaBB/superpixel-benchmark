#include "ift.h"
#include "iftExperimentUtility.h"

typedef enum connectivity_status {
  LABEL_MAP_UNCONNECTED = 0,
  LABEL_MAP_8_CONNECTED = 1,
  LABEL_MAP_4_CONNECTED = 2
} ConnectivityStatus;

int main(int argc, char *argv[])
{
  if (argc < 2)
    iftError("Usage: %s <label map> [<orig img> <output visualization prefix>]", "main", argv[0]);

  iftImage *labelMap = iftReadImageByExt(argv[1]);
  iftImage *img = argc > 2 ? iftReadImageByExt(argv[2]) : NULL;
  char *outPrefix = argc > 3 ? argv[3] : NULL;

  if (iftMinimumValue(labelMap) < 0)
    iftError("Label should only have non-negative values.", "main");

  int nLabels = iftMaximumValue(labelMap) + 1;
  ConnectivityStatus status = LABEL_MAP_4_CONNECTED;

  iftAdjRel *A[2];
  A[0] = iftIs3DImage(labelMap) ? iftSpheric(1.0) : iftCircular(1.0f);
  A[1] = iftIs3DImage(labelMap) ? iftSpheric(sqrt(3.0f)) : iftCircular(sqrt(2.0f));

  for (int adj = 0; adj < 2; ++adj) {
    // Init
    iftImage *res = iftCreateImage(labelMap->xsize, labelMap->ysize, labelMap->zsize);
    iftFIFO *Q = iftCreateFIFO(res->n);
    int *visited = iftAllocIntArray(nLabels);

    int label = 1;
    for (int pTop = 0; pTop < res->n; ++pTop) {
      // Connected component already processed
      if (res->val[pTop] > 0)
        continue;

      visited[labelMap->val[pTop]] += 1;

      // Initialization
      iftResetFIFO(Q);
      iftInsertFIFO(Q, pTop);

      // Cover entire component with BFS
      res->val[pTop] = label;
      while (!iftEmptyFIFO(Q)) {
        // Process voxel
        int p = iftRemoveFIFO(Q);

        // Visit neighbors
        iftVoxel u = iftGetVoxelCoord(res, p);
        for (int i = 1; i < A[adj]->n; ++i) {
          iftVoxel v = iftGetAdjacentVoxel(A[adj], u, i);
          if (!iftValidVoxel(res, v))
            continue;
          int q = iftGetVoxelIndex(res, v);
          if (res->val[q] == 0) {
            if (labelMap->val[p] == labelMap->val[q]) {
              res->val[q] = label;
              iftInsertFIFO(Q, q);
            }
          }
        }
      }

      label += 1;
    }


    iftColorTable *colorTb = iftCreateColorTable(1);
    for (int i = 0; i < nLabels; ++i) {
      if (visited[i] > 1) {
        status = (adj == 0 ? LABEL_MAP_8_CONNECTED : LABEL_MAP_UNCONNECTED);
        if (outPrefix != NULL) {
          char visPath[256];
          sprintf(visPath, "%s_%dconn_%d.png", outPrefix, adj == 0 ? 4 : 8, i);
          printf("Label %d is not %d-connected! Saving visualization to %s\n", i, adj == 0 ? 4 : 8, visPath);
          iftImage *tmp = iftCopyImage(img);
          iftDrawBordersSingleLabel(tmp, labelMap, i, colorTb->color[0]);
          iftWriteImageByExt(tmp, visPath);
          iftDestroyImage(&tmp);
        }
      }
    }

    iftDestroyFIFO(&Q);
    free(visited);
    iftDestroyColorTable(&colorTb);
    iftDestroyImage(&res);
    iftDestroyAdjRel(&(A[adj]));
  }

  iftDestroyImage(&labelMap);
  iftDestroyImage(&img);

  switch (status) {
    case LABEL_MAP_4_CONNECTED:
      printf("Label map is 4-connected.\n");
      return 0;
    case LABEL_MAP_8_CONNECTED:
      printf("Label map is 8-connected.\n");
      return 1;
    case LABEL_MAP_UNCONNECTED:
      printf("Label map is not connected.\n");
      return -1;
  };
}
