
#ifndef _GFT_PQUEUE16_H_
#define _GFT_PQUEUE16_H_

#include "gft_common.h"

namespace gft{
  namespace PQueue16{

    typedef struct _pqnode { 
      int  next;  //!< Next node.
      int  prev;  //!< Previous node.
      char color; //!< WHITE=0, GRAY=1, BLACK=2.
    } PQNode;
    
    typedef struct _pqdoublylinkedlists {
      PQNode *elem;  //!< All possible doubly-linked lists of the circular queue.
      int nelems;    //!< Total number of elements.
      ushort *value; //!< The value of the nodes in the graph.
    } PQDoublyLinkedLists;

    typedef struct _pqcircularqueue { 
      int  *first;   //!< List of the first elements of each doubly-linked list.
      int  *last;    //!< List of the last  elements of each doubly-linked list.
      int  nbuckets; //!< Number of buckets in the circular queue.
      ushort minvalue; //!< Minimum value of a node in queue.
      ushort maxvalue; //!< Maximum value of a node in queue.
    } PQCircularQueue;


    /**
     * \brief Priority queue by Dial implemented as proposed by 
     * A.X. Falcao with circular and growing features.
     */
    typedef struct _priorityqueue { 
      PQCircularQueue C;
      PQDoublyLinkedLists L;
      int nadded;      //!< Number of elements added.
    } PQueue16;


    PQueue16   *Create(int nbuckets, int nelems, ushort *value);
    void        Destroy(PQueue16 **Q);
    PQueue16   *Grow(PQueue16 **Q, int nbuckets);
    void        Reset(PQueue16 *Q);
    inline bool IsEmpty(PQueue16 *Q){ 
      return (Q->nadded==0); 
    }
    inline bool IsFull(PQueue16 *Q){ 
      return (Q->nadded==(Q->L).nelems); 
    }
    
    /**
     * Generic version with circular and growing features.
     */
    void   InsertElem(PQueue16 **Q, int elem);
    /**
     * Generic version with circular and growing features.
     */
    void   RemoveElem(PQueue16 *Q, int elem);
    /**
     * Generic version with circular and growing features.
     */
    void   UpdateElem(PQueue16 **Q, int elem, ushort newvalue);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMinFIFO(PQueue16 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMinLIFO(PQueue16 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMaxFIFO(PQueue16 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMaxLIFO(PQueue16 *Q);


    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastInsertElem(PQueue16 *Q, int elem);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastRemoveElem(PQueue16 *Q, int elem);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastUpdateElem(PQueue16 *Q, int elem, ushort newvalue);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMinFIFO(PQueue16 *Q);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMinLIFO(PQueue16 *Q);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMaxFIFO(PQueue16 *Q);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMaxLIFO(PQueue16 *Q);

    
  } //end PQueue16 namespace
} //end gft namespace

#endif


