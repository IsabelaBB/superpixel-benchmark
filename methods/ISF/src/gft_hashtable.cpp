
#include "gft_hashtable.h"

namespace gft{
  namespace HashTable{


    HashTable *Create(int size){
      HashTable *ht;
      int i;
      if(size<=0)
	gft::Error((char *)"Invalid size value",
		   (char *)"CreateHashTable");
      ht = (HashTable *)malloc(sizeof(HashTable));
      if(ht==NULL)
	gft::Error((char *)MSG1,
		   (char *)"CreateHashTable");
      ht->size = size;
      ht->data = (HashNode **)malloc(size*sizeof(HashNode *));
      if(ht->data==NULL)
	gft::Error((char *)MSG1,
		   (char *)"CreateHashTable");
      for(i=0; i<size; i++)
	ht->data[i]=NULL;
      
      return ht;
    }


    void Destroy(HashTable **ht){
      HashNode *node,*tmp;
      int i;
      
      if(*ht!=NULL){
	for(i=0; i<(*ht)->size; i++){
	  node = (*ht)->data[i];
	  while(node!=NULL){
	    tmp = node;
	    node = node->next;
	    free(tmp);
	  }
	}
	free((*ht)->data);
	free(*ht);
	*ht = NULL;
      }
    }


    // Function to calculate the position in the hash table
    int HashPosition(HashTable *ht, char *key){
      int h = 0, a = 127;
      
      for (; *key != '\0'; key++)
	h = (a * h + *key) % ht->size;
      
      return h;
    }


    void Insert(HashTable *ht, char *key, void *value){
      HashNode *node;
      int i;
      
      i = HashPosition(ht, key);
      node = (HashNode *)malloc(sizeof(HashNode));
      node->key = key;
      node->value = value;
      node->next = ht->data[i];
      ht->data[i] = node;
    }


    void *Search(HashTable *ht, char *key){
      HashNode *node;
      int i;
      
      i = HashPosition(ht, key);
      for(node=ht->data[i]; node!=NULL; node=node->next)
	if(strcmp(key, node->key)==0)
	  return node->value;
      
      return NULL;
    }


  } //end HashTable namespace
} //end gft namespace





