
#include "gft_filelist.h"

namespace gft{
  namespace FileList{

    FileList *Create(int cap){
      FileList *L;
      
      L = (FileList *) calloc(1,sizeof(FileList));
      if(L == NULL)
	gft::Error((char *)MSG1,(char *)"FileList::Create");
      
      L->A = ArrayList::Create(cap);
      L->n = 0;
      
      return L;
    }


    void      Destroy(FileList **L){
      FileList *aux;
      
      aux = *L;
      if(aux != NULL){
	ArrayList::Destroy(&aux->A);
	free(aux);
	*L = NULL;
      }
    }


    FileList *Read(char *filename){
      char dir[512];
      char rel[512];
      char name[512];
      char msg[512];
      FileList *L;
      FILE *fp;
      int ret;

      fp = fopen(filename,"r");
      if(fp==NULL){
	sprintf(msg,"Cannot open %s",filename);
	gft::Error((char *)msg,(char *)"FileList::Read");
      }

      L = Create(100);

      dir[0]='\0';
      while(1){
	ret = fscanf(fp," %[^\n]", rel);
	if(ret==EOF)
	  break;
	
	String::Trim(rel);
	
	//Ignore a comment line:
	if(rel[0]=='#')
	  continue;
	
	//Change base dir:
	if(rel[0]=='/'){
	  if(rel[1]=='/'){
	    rel[0] = ' ';
	    rel[1] = ' ';
	    String::Trim(rel);
	    strcpy(dir, rel);
	    continue;
	  }
	}
	
	/*Add relative or absolute path:*/
	if(dir[0]!='\0'){
	  strcpy(name, dir);
	  MergeRelativePath(name, rel);
	}
	else
	  strcpy(name, rel);
	AddFile(L, name);
      }
      fclose(fp);
      return L;
    }


    void      Write(FileList *L,
		    char *filename){
      char msg[512];
      FILE *fp;
      int i;
      
      fp = fopen(filename,"w");
      if(fp == NULL){
	sprintf(msg,"Cannot open %s",filename);
	gft::Error((char *)msg,(char *)"FileList::Write");
      }
      
      for(i=0; i<L->n; i++)
	fprintf(fp,"%s\n",GetFile(L, i));
      
      fclose(fp);
    }


    void      AddFile(FileList *L, char *file){
      char *aux,*trim;
      int s;
      
      s = strlen(file);
      trim = (char *) calloc(s+1,sizeof(char));
      strcpy(trim, file);
      
      String::Trim(trim);
      s = strlen(trim);
      
      aux = (char *) calloc(s+1,sizeof(char));
      strcpy(aux, trim);
      free(trim);
      L->n++;
      
      ArrayList::AddElement(L->A, (void *)aux);
    }
    

    char     *GetFile(FileList *L, int index){
      return (char *)ArrayList::GetElement(L->A, index);
    }
    

    bool      HasFile(FileList *L, char *file){
      char *file_i;
      int i;
      
      for(i=0; i<L->n; i++){
	file_i = (char *)ArrayList::GetElement(L->A, i);
	if(strcmp(file_i,file)==0) return true;
      }
      return false;
    }


    void AddFilesInDir(FileList *L,
		       char *dir){
      glob_t data;
      char pattern[1024];
      int s,i;
      
      s = strlen(dir);
      if(s>0)
	if(dir[s-1]=='/')
	  dir[s-1] = '\0';
      
      sprintf(pattern,"%s/*",dir);
      
      switch( glob(pattern, GLOB_MARK, NULL, &data ) ){
      case 0:
	break;
      case GLOB_NOSPACE:
	printf( "Out of memory\n" );
	break;
      case GLOB_ABORTED:
	printf( "Reading error\n" );
	break;
      case GLOB_NOMATCH:
	printf( "No files found\n" );
	break;
      default:
	break;
      }
      
      for(i=0; i<(int)data.gl_pathc; i++){
	AddFile(L, data.gl_pathv[i]);
	//printf( "%s\n", data.gl_pathv[i] );
      }
      
      globfree( &data );
    }
    
    
    void AddFilesInDirRec(FileList *L,
			  char *dir){
      FileList *T;
      char *str;
      int i,s;
      T = Create(100);
      AddFilesInDir(T, dir);
      for(i=0; i<T->n; i++){
	str = GetFile(T, i);
	s = strlen(str);
	
	if(s>0 && str[s-1]=='/')
	  AddFilesInDirRec(L, str);
	else
	  AddFile(L, str);
      }
      Destroy(&T);
    }
    

    /*It shuffles the files at random.
      Should use "void srand(unsigned int seed);" before calling.*/
    void  Randomize(FileList *L){
      int i,j;
      void *tmp;
      
      for(i=0; i<L->n; i++){
	j = gft::RandomInteger(0, L->n-1);
	
	tmp = L->A->array[i];
	L->A->array[i] = L->A->array[j];
	L->A->array[j] = tmp;
      }
    }
    
    void  Resize(FileList *L, int n){
      L->n = MIN(n, L->n);
      ArrayList::Resize(L->A, n);
    }
    
    /*Trims the capacity of this FileList instance 
      to be the list's current size.*/
    void  Trim2Size(FileList *L){
      ArrayList::Trim2Size(L->A);
    }


    void      DeleteFilesInFileList(FileList *L){
      char command[5000];
      int i;
      
      for(i=0; i<L->n; i++){
	sprintf(command,"rm %s -f",GetFile(L, i));
	system(command);
      }
    }

    
    bool FileExists(char *file){
      FILE *fp=NULL;
      fp = fopen(file,"r");
      if(fp == NULL) 
	return false;
      else{
	fclose(fp);
	return true;
      }
    }


    void   RemoveFileExtension(char *file){
      int n = strlen(file);
      
      while(n>0){
	n--;
	if(file[n]=='/') break;
	if(file[n]=='.') file[n] = '\0';
      }
    }
    

    void   RemoveFileDirectory(char *file){
      int n,i;
      
      n = strlen(file);
      while(n>0){
	n--;
	if(file[n]=='/') break;
      }
      
      if(file[n]=='/') n++;
      
      i = 0;
      while(file[n]!='\0'){
	file[i] = file[n];
	i++; n++;
      }
      file[i] = '\0';
    }


    /*the "dir" string must have enough space for the result.*/
    void MergeRelativePath(char *dir, char *rel){
      int s;
      
      s = strlen(dir);
      if(s>0 && dir[s-1]=='/')
	dir[s-1] = '\0';
      
      while(rel[0]=='.' || rel[0]=='/'){
	if(rel[1]=='.'){
	  s = strlen(dir);
	  while(s>0){
	    s--;
	    if(dir[s]=='/'){
	      dir[s]='\0';
	      break;
	    }
	  }
	  rel++;
	}
	rel++;
      }
      strcat(dir,"/");
      strcat(dir,rel);
    }



  } /*end FileList namespace*/
} /*end gft namespace*/

