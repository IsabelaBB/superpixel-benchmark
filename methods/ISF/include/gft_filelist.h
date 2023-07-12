// version 00.00.04

#ifndef _GFT_FILELIST_H_
#define _GFT_FILELIST_H_

#include "gft_common.h"
#include "gft_arraylist.h"
#include "gft_string.h"

extern "C" {
#include <glob.h>
}


namespace gft{
  namespace FileList{

    typedef struct _fileList {
      ArrayList::ArrayList *A;
      int n;   //Number of files added.
    } FileList;

    FileList *Create(int cap);
    void      Destroy(FileList **L);
    
    FileList *Read(char *filename);
    void      Write(FileList *L,
		    char *filename);

    void      AddFile(FileList *L, char *file);
    char     *GetFile(FileList *L, int index);
    bool      HasFile(FileList *L, char *file);
    
    
    void AddFilesInDir(FileList *L,
		       char *dir);
    void AddFilesInDirRec(FileList *L,
			  char *dir);

    //It shuffles the files at random.
    //Should use "void srand(unsigned int seed);" before calling.
    void  Randomize(FileList *L);
    
    void  Resize(FileList *L, int n);
    
    //Trims the capacity of this FileList instance 
    //to be the list's current size.
    void  Trim2Size(FileList *L);

    void  DeleteFilesInFileList(FileList *L);
    
    bool      FileExists(char *file);
    void      RemoveFileDirectory(char *file);
    void      RemoveFileExtension(char *file);
    void      MergeRelativePath(char *dir, char *rel);

  } /*end FileList namespace*/
} /*end gft namespace*/

#endif

