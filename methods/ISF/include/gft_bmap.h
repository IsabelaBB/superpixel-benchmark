
#ifndef _GFT_BMAP_H_
#define _GFT_BMAP_H_

#include "gft_common.h"


namespace gft{
  /**
   * \brief Common definitions and functions to manipulate a vector of booleans.
   */
  namespace BMap{

    /**
     * \brief Vector of booleans. 
     *
     * It uses one bit per boolean (i.e., size = ceil (n / 8)).
     */
    typedef struct _bmap {
      char *data;
      int N, VN;
    } BMap;
    
    BMap *Create(int n);
    void  Destroy(BMap **b);
    void  Copy(BMap *dest, BMap *src);
    void  Fill(BMap *b, int value);

    inline int   Get(BMap *b, int p){
      return ((b->data[p>>3]&(1<<(p&0x07)))!=0);
    }

    inline void  Set(BMap *b, int p, int value){
      if(value) b->data[p>>3]|=(1<<(p&0x07));
      else      b->data[p>>3]&=((~0)^(1<<(p&0x07)));
    }

    inline void Set0(BMap *b, int p){
      b->data[p>>3]&=((~0)^(1<<(p&0x07)));
    }

    inline void Set1(BMap *b, int p){
      b->data[p>>3]|=(1<<(p&0x07));
    }

    inline void  Toggle(BMap *b, int p){
      b->data[p>>3]^=(1<<(p&0x07));
    }


  } //end BMap namespace
} //end gft namespace

#endif

