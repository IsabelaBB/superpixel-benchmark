
#ifndef _GFT_FEATURES_H_
#define _GFT_FEATURES_H_

#include "gft_common.h"

/* It supports the following styles: */

#define gftFeature_ALLOCLABEL  0x000001
#define gftFeature_ALLOCINDEX  0x000002
#define gftFeature_ALLOCFV     0x000004
#define gftFeature_ALLOCDIST   0x000008
#define gftFeature_ALLOCALL    0x00000F


namespace gft{

    namespace Features{


      typedef struct _Features {
	int n;
	int nfeats;
	int *label; /* For supervised learning from labeled training data */
	int *index; /* Feature vector position in the image */
	float **fv; /* fv[0..n-1][0..nfeats-1]*/

	float **dist; /* dist[0..n-1][0..n-1] */
      } Features;


      Features *Create(int n, int nfeats, int style);
      Features *Clone(Features *f);
      void     *Destroy(Features **f);

      /*It shuffles the samples at random.*/
      void Randomize(Features *f);

      /*Returns a subset with samples of all classes in the 
	proportion given by "rate". This subset is removed from the
	original dataset. You should call "Randomize"
	first to ensure the selection of random samples.*/
      Features *RemoveSamples(Features **f,
			      float rate);
      
      /*Concatenate two datasets.*/
      Features *Merge(Features *f1, 
		      Features *f2);
      
      
    } //end Features namespace
} //end gft namespace



#endif

