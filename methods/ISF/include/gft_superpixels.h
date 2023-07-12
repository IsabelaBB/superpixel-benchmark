
#ifndef _GFT_SUPERPIXELS_H_
#define _GFT_SUPERPIXELS_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_cimage.h"
#include "gft_cimage32f.h"
#include "gft_color.h"
#include "gft_adjrel.h"
#include "gft_filtering.h"
#include "gft_queue.h"
#include "gft_heap.h"
#include "gft_scene32.h"
#include "gft_adjrel3.h"


namespace gft{
  namespace Superpixels{

    typedef struct _feature{
      int j;
      int i;
      float l;
      float a;
      float b;
    } Feature;

    typedef struct _spixelseed{
      Feature actual;
      Feature last_computed;
      bool recompute;
      bool updateseeds;
      int n;
    } SPixelSeed;

    //--------------------------

    float ColorSquaredDistance(CImage32f::CImage32f *cimg_lab, 
			       const SPixelSeed& seed, int q);
    
    float LastColorSquaredDistance(const SPixelSeed& seed);
    float LastSpatialSquaredDistance(const SPixelSeed& seed);
    
    void RunSPixelsByIFT(CImage32f::CImage32f *cimg_lab,
			 Image32::Image32 *label,
			 Heap::Heap *Q,
			 AdjRel::AdjRel *A,
			 AdjRel::AdjPxl *P,
			 Image32::Image32 *pred,
			 float *cost,
			 float alpha,
			 //float beta,
			 SPixelSeed *seeds, int nseeds);

    void RunSPixelsByDIFT(CImage32f::CImage32f *cimg_lab,
			  Image32::Image32 *label,
			  Heap::Heap *Q,
			  AdjRel::AdjRel *A,
			  AdjRel::AdjPxl *P,
			  Image32::Image32 *pred,
			  float *cost,
			  float alpha,
			  //float beta,
			  SPixelSeed *seeds, int nseeds);
    
    Image32::Image32 *ISF(CImage::CImage *cimg,
			       int k, float alpha, //float beta,
			       float err_dc, float err_ds, int itMax,
                               bool isMixed, bool isRoot/*, float *real_proc_time*/);
    Image32::Image32 *ISF(Image32::Image32 *img,
			       int k, float alpha, //float beta,
			       float err_dc, float err_ds, int itMax,
                               bool isMixed, bool isRoot/*, float *real_proc_time*/);
    //--------------------------
    float WeightedDistanceMeasure(CImage32f::CImage32f *cimg_lab, 
				  SPixelSeed seed, int i, int j,
				  float m, float S);

    Image32::Image32 *mySLIC(CImage::CImage *cimg, 
			     int k, float m);
    Image32::Image32 *mySLIC(Image32::Image32 *img, 
			     int k, float m);
    //--------------------------
    SPixelSeed *GetSPixelsSeeds(CImage32f::CImage32f *cimg_lab,
				Image32::Image32 *grad,
				int k,
				int *nseeds);

    SPixelSeed *GetMixedSamplingSecondStageSPixelsSeeds(CImage32f::CImage32f *cimg,
				Image32::Image32 *grad,
				int nsamples,
				int *final_nseeds);

    SPixelSeed *GetMixedSamplingSPixelsSeeds(CImage32f::CImage32f *cimg,
				Image32::Image32 *grad,
				int nsamples,
				int *final_nseeds);

    void Move2LowerGradient(Image32::Image32 *grad, 
			    int *j, int *i);

    void UpdateSPixelsSeeds(CImage32f::CImage32f *cimg_lab,
			    Image32::Image32 *label,
			    Image32::Image32 *pred,
			    AdjRel::AdjPxl *P,
			    SPixelSeed *seeds, int nseeds,
			    int inside, bool isRoot);

    int GetNumberOfSuperPixels(Image32::Image32 *label);


    void Postprocessing_CC(Image32::Image32 *label,
			   SPixelSeed *seeds, int nseeds);
    void Postprocessing(Image32::Image32 *label,
			SPixelSeed *seeds, int nseeds);


    //--------------------3D--------------------------------
    typedef struct _feature3D{
      int j;
      int i;
      int k;
      float l;
    } Feature3D;

    typedef struct _svoxelseed{
      Feature3D actual;
      Feature3D last_computed;
      bool recompute;
      bool updateseeds;
      int n;
    } SVoxelSeed;


    float FeatureSquaredDistance(Scene32::Scene32 *scn, 
				 const SVoxelSeed& seed, int q);
    
    float LastFeatureSquaredDistance(const SVoxelSeed& seed);
    float LastSpatialSquaredDistance(const SVoxelSeed& seed);
    
    void RunSVoxelsByIFT(Scene32::Scene32 *scn,
			 Scene32::Scene32 *label,
			 Heap::Heap *Q,
			 AdjRel3::AdjRel3 *A,
			 AdjRel3::AdjVxl *V,
			 Scene32::Scene32 *pred,
			 float *cost,
			 float alpha,
			 //float beta,
			 SVoxelSeed *seeds, int nseeds);

    void RunSVoxelsByDIFT(Scene32::Scene32 *scn,
			  Scene32::Scene32 *label,
			  Heap::Heap *Q,
			  AdjRel3::AdjRel3 *A,
			  AdjRel3::AdjVxl *V,
			  Scene32::Scene32 *pred,
			  float *cost,
			  float alpha,
			  //float beta,
			  SVoxelSeed *seeds, int nseeds);
    
    Scene32::Scene32 *ISF(Scene32::Scene32 *scn,
			       int k, float alpha, //float beta,
			       float err_df, float err_ds, int itMax,
                               bool isMixed, bool isRoot/*, float *real_proc_time*/);

    SVoxelSeed *GetSVoxelsSeeds(Scene32::Scene32 *scn,
				int k,
				int *nseeds);

    SVoxelSeed *GetMixedSamplingSecondStageSVoxelsSeeds(Scene32::Scene32 *scn,
				                        int nsamples,
			 	                        int *final_nseeds);

    SVoxelSeed *GetMixedSamplingSVoxelsSeeds(Scene32::Scene32 *scn,
				             int nsamples,
                                             int *final_nseeds);
    
    void UpdateSVoxelsSeeds(Scene32::Scene32 *scn,
			    Scene32::Scene32 *label,
			    Scene32::Scene32 *pred,
			    AdjRel3::AdjVxl *V,
			    SVoxelSeed *seeds, int nseeds,
			    int inside, bool isRoot);
    
  } //end Superpixels namespace
} //end gft namespace


#endif

