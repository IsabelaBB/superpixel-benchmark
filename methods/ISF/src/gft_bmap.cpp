
#include "gft_bmap.h"

namespace gft{
  namespace BMap{

    BMap *Create(int n) {
      BMap *x;
      x= (BMap *) malloc(sizeof(BMap));
      if (!x) return 0;
      x->N    = n;
      x->VN   = n/8;
      if (n%8) x->VN++;
      x->data = (char *) malloc(sizeof(char) * x->VN);
      if (!x->data) return 0;
      Fill(x,0);
      return x;
    }


    void  Destroy(BMap **b) {
      BMap *aux;

      aux = *b;
      if(aux != NULL){
	if(aux->data != NULL){
	  free(aux->data);
	  aux->data = NULL;
	}
	free(aux);
	*b = NULL;
      }
    }


    void   Fill(BMap *b, int value) {
      memset(b->data, value?0xff:0, b->VN);
    }


    void   Copy(BMap *dest, BMap *src) {
      int n;
      Fill(dest, 0);
      n = dest->VN;
      if (n > src->VN) n = src->VN;
      memcpy(dest->data,src->data,sizeof(char) * src->VN);
    }

    /*
    static char bmap_set[8]   = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };
    static char bmap_reset[8] = { 0xfe, 0xfd, 0xfb, 0xf7, 0xef, 0xdf, 0xbf, 0x7f };
    
    int    Get(BMap *b, int n) {
      int thebyte, thebit, value;
      
      thebyte = n >> 3;
      thebit  = n & 0x07;
      
      value = b->data[thebyte] & bmap_set[thebit];
      if (value) value=1;
      return value;
    }

    void Set(BMap *b, int n, int value) {
      int thebyte, thebit;
      
      thebyte = n >> 3;
      thebit  = n & 0x07;
      
      if (value)
	b->data[thebyte] |= bmap_set[thebit];
      else
	b->data[thebyte] &= bmap_reset[thebit];
    }

    void   Toggle(BMap *b, int n) {
      int thebyte, thebit;
      
      thebyte = n >> 3;
      thebit  = n & 0x07;

      b->data[thebyte] ^= bmap_set[thebit];
    }
    */

  } //end BMap namespace
} //end gft namespace

