#ifndef _GFT_HASHTABLE_H_
#define _GFT_HASHTABLE_H_

#include "gft_common.h"

namespace gft{
  namespace HashTable{

    /**
     * \brief A hashing table node.
     */
    typedef struct _hashnode {
      char *key;               //!< Search key.
      void *value;             //!< Associated value.
      struct _hashnode *next;  //!< Pointer to next node.
    } HashNode;


    typedef struct _hashtable {
      HashNode **data;
      int size;
    } HashTable;


    /**
     * \brief A constructor.
     */ 
    HashTable *Create(int size);

    /**
     * \brief A destructor.
     */
    void       Destroy(HashTable **ht);
    
    /**
     * \brief Insert the pair (key, value) in the hash table.
     */
    void Insert(HashTable *ht, char *key, void *value);

    /**
     * \brief Search for the value associated with a key.
     */
    void *Search(HashTable *ht, char *key);


  } //end HashTable namespace
} //end gft namespace

#endif
