
#ifndef _GFT_STRING_H_
#define _GFT_STRING_H_

#include "gft_common.h"

namespace gft{
  namespace String{

    //Removes the leading and trailing white space.
    void      Trim(char *str);
    void      SubString(char *str,
			int beginIndex,
			int endIndex);
    void      ReplaceCharacter(char *str,
			       char old_c,
			       char new_c);

  } //end String namespace
} //end gft namespace


#endif



