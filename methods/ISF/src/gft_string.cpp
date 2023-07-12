
#include "gft_string.h"

namespace gft{
  namespace String{

    void Trim(char *str){
      int i,j,s,space;
      
      s = strlen(str);
      if(s==0) return;
      
      i = 0;
      j = s-1;
      
      while(i<s){
	space = 0;
	if(str[i]==' '  ||
	   str[i]=='\t' ||
	   str[i]=='\n')
	  space = 1;
	if(!space)
	  break;
	i++;
      }
      
      while(j>=0){
	space = 0;
	if(str[j]==' '  ||
	   str[j]=='\t' ||
	   str[j]=='\n')
	  space = 1;
	if(!space)
	  break;
	j--;
      }
      
      SubString(str, i, j);
    }



    void      SubString(char *str,
			int beginIndex,
			int endIndex){
      int i,j;
      
      if(beginIndex>endIndex){
	str[0] = '\0';
	return;
      }
      
      j = 0;
      for(i=beginIndex; i<=endIndex; i++){
	str[j] = str[i];
	j++;
      }
      str[j] = '\0';
    }
    


    void      ReplaceCharacter(char *str,
			       char old_c,
			       char new_c){
      int i;
      
      if(old_c=='\0') return;
      i=0;
      while(str[i]!='\0'){
	if(str[i]==old_c)
	  str[i] = new_c;
	i++;
      }
    }
    

  } //end String namespace
} //end gft namespace

