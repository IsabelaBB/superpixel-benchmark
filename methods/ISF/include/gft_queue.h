
#ifndef _GFT_QUEUE_H_
#define _GFT_QUEUE_H_

#include "gft_common.h"

namespace gft{
  namespace Queue{

    /**
     * \brief FIFO Queue with circular and growing features.
     */
    typedef struct _queue {
      int *data;
      int get, put;
      int nbuckets;
      int nadded;   //!< Number of elements added.
    } Queue;

    Queue  *Create(int nbuckets);
    void    Destroy(Queue **Q);
    
    void    Push(Queue *Q, int p);

    /**
     * @return Returns NIL if empty.
     */
    int         Pop(Queue *Q);
    
    void        Reset(Queue *Q);
    bool        IsEmpty(Queue *Q);
    bool        IsFull(Queue *Q);

  } //end Queue namespace
} //end gft namespace

#endif


