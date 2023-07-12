
#include "gft_scene.h"

namespace gft{
  namespace Scene{

    Scene *Create(int xsize,int ysize,int zsize, int nbits, gft_SceneType type){
      Scene *scn;
      scn = (Scene *) calloc(1,sizeof(Scene));
      if(scn == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::Create");

      scn->nbits = nbits;
      scn->type = type;
      switch(type){
      case integer:
	switch(nbits){
	case  8: 
	  scn->ptr.scn8  = Scene8::Create(xsize, ysize, zsize);  
	  break;
	case 16: 
	  scn->ptr.scn16 = Scene16::Create(xsize, ysize, zsize); 
	  break;
	case 32: 
	  scn->ptr.scn32 = Scene32::Create(xsize, ysize, zsize); 
	  break;
	}
	break;
      case reals:
	switch(nbits){
	case 32: 
	  scn->ptr.scn32f = Scene32f::Create(xsize, ysize, zsize); 
	  break;
	case 64: 
	  scn->ptr.scn64f = Scene64f::Create(xsize, ysize, zsize); 
	  break;
	}
	break;
      }
      return scn;
    }

    Scene *Create(int xsize,int ysize,int zsize, int nbits){
      Scene *scn;
      scn = (Scene *) calloc(1,sizeof(Scene));
      if(scn == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::Create");

      scn->nbits = nbits;
      scn->type = integer;
      switch(nbits){
      case  8: 
	scn->ptr.scn8  = Scene8::Create(xsize, ysize, zsize);  
	break;
      case 16: 
	scn->ptr.scn16 = Scene16::Create(xsize, ysize, zsize); 
	break;
      case 32: 
	scn->ptr.scn32 = Scene32::Create(xsize, ysize, zsize); 
	break;
      }
      return scn;
    }

    Scene *Create(Scene *scn){
      Scene *new_scn;
      new_scn = (Scene *) calloc(1,sizeof(Scene));
      if(new_scn == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::Create");
      
      new_scn->nbits = scn->nbits;
      new_scn->type = scn->type;
      switch(scn->type){
      case integer:
	switch(scn->nbits){
	case  8: 
	  new_scn->ptr.scn8  = Scene8::Create(scn->ptr.scn8);
	  break;
	case 16: 
	  new_scn->ptr.scn16 = Scene16::Create(scn->ptr.scn16);
	  break;
	case 32: 
	  new_scn->ptr.scn32 = Scene32::Create(scn->ptr.scn32);
	  break;
	}
	break;
      case reals:
	switch(scn->nbits){
	case 32:
	  new_scn->ptr.scn32f  = Scene32f::Create(scn->ptr.scn32f);
	  break;
	case 64:
	  new_scn->ptr.scn64f  = Scene64f::Create(scn->ptr.scn64f);
	  break;
	}
	break;
      }
      return new_scn;
    }

    void   Destroy(Scene **scn){
      Scene *aux = *scn;
      if(aux==NULL) return;
      switch(aux->type){
      case integer:
	switch(aux->nbits){
	case  8: 
	  Scene8::Destroy(&aux->ptr.scn8);
	  break;
	case 16: 
	  Scene16::Destroy(&aux->ptr.scn16);
	  break;
	case 32: 
	  Scene32::Destroy(&aux->ptr.scn32);
	  break;
	}
	break;
      case reals:
	switch(aux->nbits){
	case 32: 
	  Scene32f::Destroy(&aux->ptr.scn32f);
	  break;
	case 64: 
	  Scene64f::Destroy(&aux->ptr.scn64f);
	  break;
	}
	break;
      }
      free(aux);
      *scn = NULL;
    }


    void   Copy(Scene *dest, Scene *src){
      if(dest->type==src->type) {
	switch( dest->type ) {
	case integer:
	  if(dest->nbits==src->nbits){
	    switch(dest->nbits){
	    case  8:
	      Scene8::Copy(dest->ptr.scn8, 
			   src->ptr.scn8);
	      break;
	    case 16:
	      Scene16::Copy(dest->ptr.scn16, 
			    src->ptr.scn16);
	      break;
	    case 32:
	      Scene32::Copy(dest->ptr.scn32, 
			    src->ptr.scn32);
	      break;
	    }
	  }
	  if(dest->nbits>src->nbits){
	    int p,n;
	    n = GetNumberOfVoxels(dest);
	    for(p=0; p<n; p++)
	      SetValue(dest, p, GetValue(src, p));
	  }
	  else
	    gft::Error((char *)"Incompatible number of bits",
		       (char *)"Scene::Copy");
	  break;
	case reals:
	  if(dest->nbits==src->nbits){
	    switch(dest->nbits){
	    case 32:
	      Scene32f::Copy(dest->ptr.scn32f, 
			     src->ptr.scn32f);
	      break;
	    case 64:
	      Scene64f::Copy(dest->ptr.scn64f, 
			     src->ptr.scn64f);
	      break;
	    }
	  }
	  if(dest->nbits>src->nbits){
	    int p,n;
	    n = GetNumberOfVoxels(dest);
	    for(p=0; p<n; p++)
	      SetValue(dest, p, GetValue(src, p));
	  }
	  else
	    gft::Error((char *)"Incompatible number of bits",
		       (char *)"Scene::Copy");
	  break;
	}
      }
    }


    void   Copy(Scene *dest, Scene *src, Voxel v){
      if(dest->nbits==src->nbits){
	switch(dest->nbits){
	case  8:
	  Scene8::Copy(dest->ptr.scn8, 
		       src->ptr.scn8, v);
	  break;
	case 16:
	  Scene16::Copy(dest->ptr.scn16, 
			src->ptr.scn16, v);
	  break;
	case 32:
	  Scene32::Copy(dest->ptr.scn32, 
			src->ptr.scn32, v);
	  break;
	}
      }
      else
	gft::Error((char *)"Incompatible number of bits",
		   (char *)"Scene::Copy");
    }


    Scene *Clone(Scene *scn){
      Scene *aux;

      aux = (Scene *) calloc(1,sizeof(Scene));
      if(aux == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::Clone");

      aux->nbits = scn->nbits;
      switch(scn->nbits){
      case  8: 
	aux->ptr.scn8  = Scene8::Clone(scn->ptr.scn8);
	break;
      case 16:
	aux->ptr.scn16 = Scene16::Clone(scn->ptr.scn16);
	break;
      case 32: 
	aux->ptr.scn32 = Scene32::Clone(scn->ptr.scn32);
	break;
      }
      return aux;
    }


    Scene *SubScene(Scene *scn, Voxel l, Voxel h){
      return SubScene(scn,
		      l.c.x, l.c.y, l.c.z,
		      h.c.x, h.c.y, h.c.z);
    }


    Scene *SubScene(Scene *scn,
		    int xl, int yl, int zl,
		    int xh, int yh, int zh){
      Scene *aux;
      aux = (Scene *) calloc(1,sizeof(Scene));
      if(aux == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::SubScene");

      aux->nbits = scn->nbits;
      switch(scn->nbits){
      case  8:
	aux->ptr.scn8  = Scene8::SubScene(scn->ptr.scn8, 
					  xl,yl,zl, 
					  xh,yh,zh);
	break;
      case 16:
	aux->ptr.scn16 = Scene16::SubScene(scn->ptr.scn16, 
					   xl,yl,zl, 
					   xh,yh,zh);
	break;
      case 32:
	aux->ptr.scn32 = Scene32::SubScene(scn->ptr.scn32, 
					   xl,yl,zl, 
					   xh,yh,zh);
	break;
      }
      return aux;
    }


    void   Fill(Scene *scn, int value){
      switch(scn->nbits){
      case  8: 
	Scene8::Fill(scn->ptr.scn8, (uchar)value);
	break;
      case 16:
	Scene16::Fill(scn->ptr.scn16, (ushort)value);
	break;
      case 32: 
	Scene32::Fill(scn->ptr.scn32, value);
	break;
      }
    }


    Scene *Read(char *filename){
      return NULL;
    }


    void   Write(Scene *scn, char *filename){
      switch(scn->nbits){
      case  8: 
	Scene8::Write(scn->ptr.scn8, filename);
	break;
      case 16:
	Scene16::Write(scn->ptr.scn16, filename);
	break;
      case 32: 
	Scene32::Write(scn->ptr.scn32, filename);
	break;
      }
    }


    void   SetValue(Scene *scn, int p, int value){
      switch(scn->nbits){
      case  8: 
	(scn->ptr.scn8)->data[p] = (uchar)value;
	break;
      case 16:
	(scn->ptr.scn16)->data[p] = (ushort)value;
	break;
      case 32: 
	(scn->ptr.scn32)->data[p] = value;
	break;
      }
    }

    void   SetValue(Scene *scn, int p, double value){
      switch(scn->nbits){
      case 32: 
	(scn->ptr.scn32f)->data[p] = value;
	break;
      case 64: 
	(scn->ptr.scn64f)->data[p] = value;
	break;
      }
    }

    int    GetValue(Scene *scn, Voxel v){
      switch(scn->nbits){
      case  8: 
	return (int)Scene8::GetValue(scn->ptr.scn8, v);
      case 16:
	return (int)Scene16::GetValue(scn->ptr.scn16, v);
      case 32:
	return Scene32::GetValue(scn->ptr.scn32, v);
      }
      return NIL;
    }

    double    GetValue64f(Scene *scn, Voxel v){
      switch(scn->nbits){
      case 32:
	return Scene32f::GetValue(scn->ptr.scn32f, v);
      case 64:
	return Scene64f::GetValue(scn->ptr.scn64f, v);
      }
      return NIL;
    }

    int    GetValue(Scene *scn, int p){
      switch(scn->nbits){
      case  8: 
	return (int)(scn->ptr.scn8)->data[p];
      case 16:
	return (int)(scn->ptr.scn16)->data[p];
      case 32:
	return (scn->ptr.scn32)->data[p];
      }
      return NIL;
    }

    double    GetValue64f(Scene *scn, int p){
      switch(scn->nbits){
      case 32:
	return (scn->ptr.scn32f)->data[p];
      case 64:
	return (scn->ptr.scn64f)->data[p];
      }
      return NIL;
    }

    int    GetValue(Scene *scn, int x, int y, int z){
      switch(scn->nbits){
      case  8: 
	return (int)(scn->ptr.scn8)->array[z][y][x];
      case 16:
	return (int)(scn->ptr.scn16)->array[z][y][x];
      case 32:
	return (scn->ptr.scn32)->array[z][y][x];
      }
      return NIL;
    }

    double    GetValue64f(Scene *scn, int x, int y, int z){
      switch(scn->nbits){
      case 32: 
	return (double)(scn->ptr.scn32f)->array[z][y][x];
      case 64:
	return (double)(scn->ptr.scn64f)->array[z][y][x];
      }
      return NIL;
    }

    int    GetValue_nn(Scene *scn, float x, float y, float z){
      switch(scn->nbits){
      case  8: 
	return (int)Scene8::GetValue_nn(scn->ptr.scn8, x,y,z);
      case 16:
	return (int)Scene16::GetValue_nn(scn->ptr.scn16, x,y,z);
      case 32:
	return Scene32::GetValue_nn(scn->ptr.scn32, x,y,z);
      }
      return NIL;
    }

    double   GetValue64f_nn(Scene *scn, float x, float y, float z){
      switch(scn->nbits){
      case 32:
	return (double)Scene32::GetValue_nn(scn->ptr.scn32, x,y,z);
      case 64:
	return Scene64f::GetValue_nn(scn->ptr.scn64f, x,y,z);
      }
      return NIL;
    }

    int    GetNumberOfVoxels(Scene *scn){
      switch(scn->type) {
      case integer:
	switch(scn->nbits){
	case  8:
	  return (scn->ptr.scn8)->n;
	case 16:
	  return (scn->ptr.scn16)->n;
	case 32:
	  return (scn->ptr.scn32)->n;
	}
	break;
      case reals:
	switch(scn->nbits){
	case 32:
	  return (scn->ptr.scn32f)->n;
	case 64:
	  return (scn->ptr.scn64f)->n;
	}
	break;
      }
      return NIL;
    }
  // Parei aqui.
    int    GetAddressX(Scene *scn, int p){
      switch(scn->nbits){
      case  8:
	return Scene8::GetAddressX(scn->ptr.scn8, p);
      case 16:
	return Scene16::GetAddressX(scn->ptr.scn16, p);
      case 32:
	return Scene32::GetAddressX(scn->ptr.scn32, p);
      }
      return NIL;
    }

    int    GetAddressY(Scene *scn, int p){
      switch(scn->nbits){
      case  8:
	return Scene8::GetAddressY(scn->ptr.scn8, p);
      case 16:
	return Scene16::GetAddressY(scn->ptr.scn16, p);
      case 32:
	return Scene32::GetAddressY(scn->ptr.scn32, p);
      }
      return NIL;
    }

    int    GetAddressZ(Scene *scn, int p){
      switch(scn->nbits){
      case  8:
	return Scene8::GetAddressZ(scn->ptr.scn8, p);
      case 16:
	return Scene16::GetAddressZ(scn->ptr.scn16, p);
      case 32:
	return Scene32::GetAddressZ(scn->ptr.scn32, p);
      }
      return NIL;
    }

    int    GetVoxelAddress(Scene *scn, Voxel v){
      switch(scn->nbits){
      case  8:
	return Scene8::GetVoxelAddress(scn->ptr.scn8, v);
      case 16:
	return Scene16::GetVoxelAddress(scn->ptr.scn16, v);
      case 32:
	return Scene32::GetVoxelAddress(scn->ptr.scn32, v);
      }
      return NIL;
    }

    int    GetVoxelAddress(Scene *scn, int x, int y, int z){
      switch(scn->nbits){
      case  8:
	return Scene8::GetVoxelAddress(scn->ptr.scn8, x,y,z);
      case 16:
	return Scene16::GetVoxelAddress(scn->ptr.scn16, x,y,z);
      case 32:
	return Scene32::GetVoxelAddress(scn->ptr.scn32, x,y,z);
      }
      return NIL;
    }

    bool   IsValidVoxel(Scene *scn, int x, int y, int z){
      switch(scn->nbits){
      case  8:
	return Scene8::IsValidVoxel(scn->ptr.scn8, x,y,z);
      case 16:
	return Scene16::IsValidVoxel(scn->ptr.scn16, x,y,z);
      case 32:
	return Scene32::IsValidVoxel(scn->ptr.scn32, x,y,z);
      }
      return false;
    }

    bool   IsValidVoxel(Scene *scn, Voxel v){
      switch(scn->nbits){
      case  8:
	return Scene8::IsValidVoxel(scn->ptr.scn8, v);
      case 16:
	return Scene16::IsValidVoxel(scn->ptr.scn16, v);
      case 32:
	return Scene32::IsValidVoxel(scn->ptr.scn32, v);
      }
      return false;
    }
  
    int    GetMaximumValue(Scene *scn){
      switch(scn->nbits){
      case  8:
	return Scene8::GetMaximumValue(scn->ptr.scn8);
      case 16:
	return Scene16::GetMaximumValue(scn->ptr.scn16);
      case 32:
	return Scene32::GetMaximumValue(scn->ptr.scn32);
      }
      return NIL;
    }

    int    GetMinimumValue(Scene *scn){
      switch(scn->nbits){
      case  8:
	return Scene8::GetMinimumValue(scn->ptr.scn8);
      case 16:
	return Scene16::GetMinimumValue(scn->ptr.scn16);
      case 32:
	return Scene32::GetMinimumValue(scn->ptr.scn32);
      }
      return NIL;
    }


    Scene *MBB(Scene *scn){
      Scene *aux;
      aux = (Scene *) calloc(1,sizeof(Scene));
      if(aux == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::MBB");

      aux->nbits = scn->nbits;
      switch(scn->nbits){
      case  8:
	aux->ptr.scn8  = Scene8::MBB(scn->ptr.scn8);
	break;
      case 16:
	aux->ptr.scn16 = Scene16::MBB(scn->ptr.scn16);
	break;
      case 32:
	aux->ptr.scn32 = Scene32::MBB(scn->ptr.scn32);
	break;
      }
      return aux;
    }


    void   MBB(Scene *scn, Voxel *l, Voxel *h){
      switch(scn->nbits){
      case  8:
	Scene8::MBB(scn->ptr.scn8, l, h);
	break;
      case 16:
	Scene16::MBB(scn->ptr.scn16, l, h);
	break;
      case 32:
	Scene32::MBB(scn->ptr.scn32, l, h);
	break;
      }
    }

    
    Scene *AddFrame(Scene *scn,  int sz, int value){
      Scene *aux;
      aux = (Scene *) calloc(1,sizeof(Scene));
      if(aux == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::AddFrame");

      aux->nbits = scn->nbits;
      switch(scn->nbits){
      case  8:
	aux->ptr.scn8  = Scene8::AddFrame(scn->ptr.scn8, sz, (uchar)value);
	break;
      case 16:
	aux->ptr.scn16 = Scene16::AddFrame(scn->ptr.scn16, sz, (ushort)value);
	break;
      case 32:
	aux->ptr.scn32 = Scene32::AddFrame(scn->ptr.scn32, sz, value);
	break;
      }
      return aux;
    }


    Scene *RemFrame(Scene *fscn, int sz){
      Scene *aux;
      aux = (Scene *) calloc(1,sizeof(Scene));
      if(aux == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::RemFrame");

      aux->nbits = fscn->nbits;
      switch(fscn->nbits){
      case  8:
	aux->ptr.scn8  = Scene8::RemFrame(fscn->ptr.scn8, sz);
	break;
      case 16:
	aux->ptr.scn16 = Scene16::RemFrame(fscn->ptr.scn16, sz);
	break;
      case 32:
	aux->ptr.scn32 = Scene32::RemFrame(fscn->ptr.scn32, sz);
	break;
      }
      return aux;
    }



  } //end Scene namespace
} //end gft namespace

