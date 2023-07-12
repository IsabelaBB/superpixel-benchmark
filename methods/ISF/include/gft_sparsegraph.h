
#ifndef _GFT_SPARSEGRAPH_H_
#define _GFT_SPARSEGRAPH_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_cimage.h"
#include "gft_adjrel.h"

namespace gft{
  namespace SparseGraph{

#define DISSIMILARITY 0
#define CAPACITY      1

    // Image Sparse Graph:
    // n-links connect pairs of neighboring pixels.
    typedef struct _sparseGraph {
      int type;
      int Wmax;
      int **n_link;
      int ncols, nrows;
      AdjRel::AdjRel *A;
    } SparseGraph;
    
    
    SparseGraph *ByAbsDiff(Image32::Image32 *img, float r);
    SparseGraph *ByAccAbsDiff(Image32::Image32 *img, float r, float R);
    SparseGraph *ByAccAbsDiff(CImage::CImage *cimg, float r, float R);
    SparseGraph *ByWeightImage(Image32::Image32 *W, float r);

    SparseGraph *ByHomogeneityAffinity(Image32::Image32 *img, float r);
    SparseGraph *ByCHomogeneityAffinity(CImage::CImage *cimg, float r);


    SparseGraph   *WeightedMean(SparseGraph *G1, 
				SparseGraph *G2,
				float w2);
    SparseGraph   *ByLevel(Image32::Image32 *obj, int T, float r);
    SparseGraph   *ByExclusion(Image32::Image32 *obj,
			       Image32::Image32 *bkg, float r);

    void LinearStretch(SparseGraph *sg, 
		       int f1, int f2, 
		       int g1, int g2);
    
    //------ SparseGraph Functions ------------
    
    int            get_edge_index(int p, int q, 
				  SparseGraph *g);
    
    void           Destroy(SparseGraph **g);
    SparseGraph   *Create(int ncols, int nrows, 
			  AdjRel::AdjRel *A);
    SparseGraph   *Clone(SparseGraph *g);

    SparseGraph   *ReadFromTxt(char *filename);
    void           Write2Txt(SparseGraph *sg, 
			     char *filename);
    
    void           ChangeType(SparseGraph *g, 
			      int type);
    
    void           Pow(SparseGraph *sg, int power, int max);
    SparseGraph   *Convert2Symmetric(SparseGraph *g,
				     int method);
    
    void  SuppressZeroWeightedArcs(SparseGraph *sg);
    
    void Transpose(SparseGraph *sg);
    
    Image32::Image32 *ArcWeightImage(SparseGraph *sg);
    
    void Orient2Digraph(SparseGraph *sg, 
			Image32::Image32 *img, float per);
    
    void Orient2DigraphInner(SparseGraph *sg, 
			     Image32::Image32 *P_sum);
    void Orient2DigraphOuter(SparseGraph *sg, 
			     Image32::Image32 *P_sum);

  } //end SparseGraph namespace
} //end gft namespace

#endif

