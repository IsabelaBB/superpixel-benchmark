#ifndef _GFT_PQUEUE32_H_
#define _GFT_PQUEUE32_H_

#include "gft_common.h"

namespace gft{
  namespace PQueue32{

    typedef struct _pqnode { 
      int  next;  //!< Next node.
      int  prev;  //!< Previous node.
      char color; //!< WHITE=0, GRAY=1, BLACK=2.
    } PQNode;
    
    typedef struct _pqdoublylinkedlists {
      PQNode *elem; //!< All possible doubly-linked lists of the circular queue.
      int nelems;   //!< Total number of elements.
      int *value;   //!< The value of the nodes in the graph.
    } PQDoublyLinkedLists; 

    typedef struct _pqcircularqueue { 
      int  *first;   //!< List of the first elements of each doubly-linked list.
      int  *last;    //!< List of the last  elements of each doubly-linked list.
      int  nbuckets; //!< Number of buckets in the circular queue.
      int  minvalue; //!< Minimum value of a node in queue.
      int  maxvalue; //!< Maximum value of a node in queue.
    } PQCircularQueue;


    /**
     * \brief Priority queue by Dial implemented as proposed by 
     * A.X. Falcao with circular and growing features.
     */
    typedef struct _priorityqueue { 
      PQCircularQueue C;
      PQDoublyLinkedLists L;
      int nadded;      //!< Number of elements added.
    } PQueue32;


    PQueue32   *Create(int nbuckets, int nelems, int *value);
    void        Destroy(PQueue32 **Q);
    PQueue32   *Grow(PQueue32 **Q, int nbuckets);
    void        Reset(PQueue32 *Q);
    inline bool IsEmpty(PQueue32 *Q){ 
      return (Q->nadded==0); 
    }
    inline bool IsFull(PQueue32 *Q){ 
      return (Q->nadded==(Q->L).nelems); 
    }
    

    /**
     * Generic version with circular and growing features.
     */
    void   InsertElem(PQueue32 **Q, int elem);
    /**
     * Generic version with circular and growing features.
     */
    void   RemoveElem(PQueue32 *Q, int elem);
    /**
     * Generic version with circular and growing features.
     */
    void   UpdateElem(PQueue32 **Q, int elem, int newvalue);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMinFIFO(PQueue32 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMinLIFO(PQueue32 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMaxFIFO(PQueue32 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMaxLIFO(PQueue32 *Q);


    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastInsertElem(PQueue32 *Q, int elem);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastRemoveElem(PQueue32 *Q, int elem);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastUpdateElem(PQueue32 *Q, int elem, int newvalue);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMinFIFO(PQueue32 *Q);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMinLIFO(PQueue32 *Q);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMaxFIFO(PQueue32 *Q);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMaxLIFO(PQueue32 *Q);

    void   FastInsertElemAsFirst(PQueue32 *Q, int elem);
    
  } //end PQueue32 namespace
} //end gft namespace

#endif


