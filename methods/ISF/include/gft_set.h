
#ifndef _GFT_SET_H_
#define _GFT_SET_H_

#include "gft_common.h"
#include "gft_bmap.h"

#include "gft_image32.h"
#include "gft_adjrel.h"

namespace gft{
  namespace Set{

    typedef struct _set {
      int elem;
      struct _set *next;
    } Set;

    
    Set *Create();
    Set *Create(gft::Image32::Image32 *bin,
		gft::AdjRel::AdjRel *A);
    Set *Create(gft::Image32::Image32 *img);
    
    void Destroy(Set **S);
    Set *Clone(Set *S);

    void Insert(Set **S, int elem);
    int  Remove(Set **S);
    void RemoveElem(Set **S, int elem);
    bool IsInSet(Set *S, int elem);
    int  MinimumValue(Set *S);
    int  MaximumValue(Set *S);
    void Convert2DisjointSets(Set **S1,
			      Set **S2);
    int  GetNElems(Set *S);
    
    /**
     * \brief Merge two sets. 
     *
     * The next field of the last element of set S 
     * points to the first element of set T. 
     * T does not change.
     */
    void Merge(Set **S, Set **T);

  } //end Set namespace
} //end gft namespace



namespace gft{
  namespace Image32{

    void DrawSet(Image32 *img,
		 gft::Set::Set *S, 
		 int value);
    

  } //end Image32 namespace
} //end gft namespace


#include "gft_cimage.h"

namespace gft{
  namespace CImage{
    
    void DrawSet(CImage *img,
		 gft::Set::Set *S, 
		 int color);

    
  } //end CImage namespace
} //end gft namespace



#endif

