// version 00.00.02

#ifndef _GFT_ARRAYLIST_H_
#define _GFT_ARRAYLIST_H_

#include "gft_common.h"

namespace gft{
  namespace ArrayList{

    typedef struct _ArrayList {
      void **array;
      int cap; //Capacity.
      int n;   //Number of objects added.
      void (*clean)(void**); //Clean function
    } ArrayList;

    ArrayList *Create(int cap);
    void       Destroy(ArrayList **A);
    
    void       SetCleanFunc(ArrayList *A,
			    void (*clean)(void**));

    void       AddElement(ArrayList *A, 
			  void *elem);
    void      *GetElement(ArrayList *A, 
			  int index);
    void       DelElement(ArrayList *A, 
			  int index);
    void       DelElement(ArrayList *A,
			  void **elem);
    
    void       Resize(ArrayList *A, int n);
    
    //Trims the capacity of this ArrayList instance 
    //to be the list's current size.
    void       Trim2Size(ArrayList *A);
    
  } /*end ArrayList namespace*/
} /*end gft namespace*/


#endif

