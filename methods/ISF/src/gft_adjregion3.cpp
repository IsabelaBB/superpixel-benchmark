
#include "gft_adjregion3.h"

namespace gft{
  namespace AdjRegion3{

    AdjRegion3 *Create(int n){
      AdjRegion3 *C=NULL;

      C = (AdjRegion3 *) calloc(1,sizeof(AdjRegion3));
      if(C != NULL){
	if(n>0){
	  C->d = (Displacement3 *)_mm_malloc(sizeof(Displacement3)*n, 16);
	  if(C->d==NULL)
	    gft::Error((char *)MSG1,
		       (char *)"AdjRegion3::Create");
	  memset((void *)C->d, 0, sizeof(Displacement3)*n);
	  C->n = n;
	}
	else{
	  C->d = NULL; 
	  C->n  = 0;
	}
	C->dp = NULL;
	C->xsize  = -1;
	C->xysize = -1;
	C->max[0] = INT_MIN;
	C->max[1] = INT_MIN;
	C->max[2] = INT_MIN;
	C->min[0] = INT_MAX;
	C->min[1] = INT_MAX;
	C->min[2] = INT_MAX;
      }
      else
	gft::Error((char *)MSG1,
		   (char *)"AdjRegion3::Create");
      return(C);
    }


    AdjRegion3 *Create(AdjRel3::AdjRel3 *A){
      AdjRegion3 *C;
      int i;

      C = Create(A->n);
      for(i=0; i < A->n; i++){
	C->d[i].v = A->d[i].v;
      }
      return C;
    }


    AdjRegion3 *Create(Scene8::Scene8 *mask,
		       Voxel Ref){
      AdjRegion3 *C;
      Voxel v;
      int i,p,n,nelems=0;

      n = mask->n;
      for(p=0; p<n; p++)
	if(mask->data[p]>0)
	  nelems++;

      C = Create(nelems);
      i = 0;
      for(p=0; p<n; p++){
	if(mask->data[p]>0){
	  v.c.x = Scene8::GetAddressX(mask, p);
	  v.c.y = Scene8::GetAddressY(mask, p);
	  v.c.z = Scene8::GetAddressZ(mask, p);
	  C->d[i].v = v.v - Ref.v;
	  i++;
	}
      }
      return C;
    }


    void Destroy(AdjRegion3 **adjreg){
      AdjRegion3 *aux;
      aux = *adjreg;
      if(aux != NULL){
	if(aux->d  != NULL) _mm_free(aux->d);
	if(aux->dp != NULL) gft::FreeIntArray(&aux->dp);
	free(aux);
	*adjreg = NULL;
      }   
    }


    AdjRegion3 *Clone(AdjRegion3 *adjreg){
      AdjRegion3 *C;
      int i;

      C = Create(adjreg->n);
      for(i=0; i < adjreg->n; i++){
	C->d[i] = adjreg->d[i];
      }
      return C;
    }


    AdjRegion3 *Merge(AdjRegion3 *r1, AdjRegion3 *r2){
      int dx_max,dy_max,dz_max;
      int dx_min,dy_min,dz_min;
      Scene8::Scene8 *mask=NULL;
      AdjRegion3 *r3=NULL;
      Voxel v;
      
      RefreshLimits(r1);
      RefreshLimits(r2);
      dx_max = MAX( r1->max[0], r2->max[0] );
      dy_max = MAX( r1->max[1], r2->max[1] );
      dz_max = MAX( r1->max[2], r2->max[2] );
      dx_min = MIN( r1->min[0], r2->min[0] );
      dy_min = MIN( r1->min[1], r2->min[1] );
      dz_min = MIN( r1->min[2], r2->min[2] );

      mask = Scene8::Create(abs(dx_min)+dx_max+1,
			    abs(dy_min)+dy_max+1,
			    abs(dz_min)+dz_max+1);
      v.c.x = abs(dx_min);
      v.c.y = abs(dy_min);
      v.c.z = abs(dz_min);
      Draw(r1, mask, v, 1);
      Draw(r2, mask, v, 1);
      r3 = Create(mask, v);
      Scene8::Destroy(&mask);
      return r3;
    }



    Scene8::Scene8 *Export2Mask(AdjRegion3 *adjreg){
      int i,dx,dy,dz,dxmax,dymax,dzmax;
      Scene8::Scene8 *mask;
      Voxel v;
      dxmax = dymax = dzmax = 0;
      for(i=0; i<adjreg->n; i++){
	dx = abs(adjreg->d[i].axis.x);
	dy = abs(adjreg->d[i].axis.y);
	dz = abs(adjreg->d[i].axis.z);
	if(dx>dxmax) dxmax = dx;
	if(dy>dymax) dymax = dy;
	if(dz>dzmax) dzmax = dz;
      }
      mask = Scene8::Create(dxmax*2+1,
			    dymax*2+1,
			    dzmax*2+1);
      v.c.x = dxmax;
      v.c.y = dymax;
      v.c.z = dzmax;
      Draw(adjreg, mask, v, 1);
      return mask;
    }


    void   Draw(AdjRegion3 *adjreg,
		Scene8::Scene8 *scn,
		Voxel u,
		uchar val){
      Voxel v;
      int i;
      for(i=0; i<adjreg->n; i++){
	v.v = u.v + adjreg->d[i].v;
	if(Scene8::IsValidVoxel(scn,v))
	  scn->array[v.c.z][v.c.y][v.c.x] = val;
      }
    }


    void    Draw(AdjRegion3 *adjreg,
		 Scene16::Scene16 *scn,
		 Voxel u,
		 ushort val){
      Voxel v;
      int i;
      for(i=0; i<adjreg->n; i++){
	v.v = u.v + adjreg->d[i].v;
	if(Scene16::IsValidVoxel(scn,v))
	  scn->array[v.c.z][v.c.y][v.c.x] = val;
      }
    }


    void    Draw(AdjRegion3 *adjreg,
		 Scene32::Scene32 *scn,
		 Voxel u,
		 int val){
      Voxel v;
      int i;
      for(i=0; i<adjreg->n; i++){
	v.v = u.v + adjreg->d[i].v;
	if(Scene32::IsValidVoxel(scn,v))
	  scn->array[v.c.z][v.c.y][v.c.x] = val;
      }
    }


    void   DrawOpt(AdjRegion3 *adjreg,
		   Scene8::Scene8 *scn,
		   int p, uchar val){
      int i,q;
      Optimize(adjreg, scn);
      for(i=0; i<adjreg->n; i++){
	q = p + adjreg->dp[i];
	scn->data[q] = val;
      }
    }


    void    DrawOpt(AdjRegion3 *adjreg,
		    Scene16::Scene16 *scn,
		    int p, ushort val){
      int i,q;
      Optimize(adjreg, scn);
      for(i=0; i<adjreg->n; i++){
	q = p + adjreg->dp[i];
	scn->data[q] = val;
      }
    }


    void    DrawOpt(AdjRegion3 *adjreg,
		    Scene32::Scene32 *scn,
		    int p, int val){
      int i,q;
      Optimize(adjreg, scn);
      for(i=0; i<adjreg->n; i++){
	q = p + adjreg->dp[i];
	scn->data[q] = val;
      }
    }


    void    Optimize(AdjRegion3 *adjreg,
		     int xsize, int ysize){
      AdjRegion3 *C = adjreg;
      int xysize,i;
      xysize = xsize*ysize;

      if(C->xsize==xsize && C->xysize==xysize && C->dp!=NULL)
	return;
  
      if(C->dp!=NULL) gft::FreeIntArray(&C->dp);
      C->dp = gft::AllocIntArray(C->n);
      
      for(i=0; i<C->n; i++)
	C->dp[i] = (C->d[i].axis.x + 
		    C->d[i].axis.y*xsize + 
		    C->d[i].axis.z*xysize);
      C->xsize  = xsize;
      C->xysize = xysize;
    }


    void    Optimize(AdjRegion3 *adjreg,
		     Scene8::Scene8 *scn){
      Optimize(adjreg, scn->xsize, scn->ysize);
    }

    void    Optimize(AdjRegion3 *adjreg,
		     Scene16::Scene16 *scn){
      Optimize(adjreg, scn->xsize, scn->ysize);
    }

    void    Optimize(AdjRegion3 *adjreg,
		     Scene32::Scene32 *scn){
      Optimize(adjreg, scn->xsize, scn->ysize);
    }


    void    RefreshLimits(AdjRegion3 *adjreg){
      int dx_min, dy_min, dz_min;
      int dx_max, dy_max, dz_max;
      int i,dx,dy,dz;

      dx_min = dy_min = dz_min = INT_MAX;
      dx_max = dy_max = dz_max = INT_MIN;
      for(i=0; i<adjreg->n; i++){
	dx = adjreg->d[i].axis.x;
	dy = adjreg->d[i].axis.y;
	dz = adjreg->d[i].axis.z;
	if(dx > dx_max) dx_max = dx;
	if(dy > dy_max) dy_max = dy;
	if(dz > dz_max) dz_max = dz;
	if(dx < dx_min) dx_min = dx;
	if(dy < dy_min) dy_min = dy;
	if(dz < dz_min) dz_min = dz;
      }
      adjreg->max[0] = dx_max;
      adjreg->max[1] = dy_max;
      adjreg->max[2] = dz_max;
      adjreg->min[0] = dx_min;
      adjreg->min[1] = dy_min;
      adjreg->min[2] = dz_min;
    }


    void    GetLimits(AdjRegion3 *adjreg,
		      int *dx_min, int *dy_min, int *dz_min,
		      int *dx_max, int *dy_max, int *dz_max){
      if(adjreg->max[0]<adjreg->min[0])
	RefreshLimits(adjreg);
  
      *dx_max = adjreg->max[0];
      *dy_max = adjreg->max[1];
      *dz_max = adjreg->max[2];
      *dx_min = adjreg->min[0];
      *dy_min = adjreg->min[1];
      *dz_min = adjreg->min[2];
    }


    bool    FitInside(AdjRegion3 *adjreg, Voxel vx,
		      int xsize, int ysize, int zsize,
		      int sz){
      int max[3],min[3];
      GetLimits(adjreg,
		&min[0], &min[1], &min[2],
		&max[0], &max[1], &max[2]);

      if(vx.c.x < sz-min[0] ||
	 vx.c.y < sz-min[1] ||
	 vx.c.z < sz-min[2] ||
	 vx.c.x+max[0]>=xsize-sz ||
	 vx.c.y+max[1]>=ysize-sz ||
	 vx.c.z+max[2]>=zsize-sz)
	return false;
      else
	return true;
    }


    bool    FitInside(AdjRegion3 *adjreg, Voxel vx,
		      Scene8::Scene8 *scn, int sz){
      return FitInside(adjreg, vx, 
		       scn->xsize, scn->ysize, scn->zsize, 
		       sz);
    }

    bool    FitInside(AdjRegion3 *adjreg, Voxel vx,
		      Scene16::Scene16 *scn, int sz){
      return FitInside(adjreg, vx, 
		       scn->xsize, scn->ysize, scn->zsize, 
		       sz);
    }

    bool    FitInside(AdjRegion3 *adjreg, Voxel vx,
		      Scene32::Scene32 *scn, int sz){
      return FitInside(adjreg, vx, 
		       scn->xsize, scn->ysize, scn->zsize, 
		       sz);
    }


    float   InnerMean(AdjRegion3 *adjreg,
		      Voxel vx,
		      Scene16::Scene16 *scn){
      Voxel v;
      float sum=0.0;
      int i,p;
      for(i=0; i<adjreg->n; i++){
	v.v = vx.v + adjreg->d[i].v;
	if(Scene16::IsValidVoxel(scn,v)){
	  p = Scene16::GetVoxelAddress(scn,v);
	  sum += scn->data[p];
	}
      }
      return (sum/adjreg->n);
    }



    float   InnerSum(AdjRegion3 *adjreg,
		     Voxel vx,
		     Scene16::Scene16 *scn){
      Voxel v;
      float sum=0.0;
      int i,p;
      for(i=0; i<adjreg->n; i++){
	v.v = vx.v + adjreg->d[i].v;
	if(Scene16::IsValidVoxel(scn,v)){
	  p = Scene16::GetVoxelAddress(scn,v);
	  sum += scn->data[p];
	}
      }
      return sum;
    }


    //-----------multithreading-----------------
    /*
    extern "C"{
    #include <pthread.h>
    #include "processors.h"
    }

    /// \cond
    typedef struct _argdraw {
      AdjRegion3 *adjreg;
      Scene8::Scene8   *scn8;
      Scene16::Scene16 *scn16;
      Scene32::Scene32 *scn32;
      Voxel u;
      int p;
      int val;
      int i; 
      int j;
    } ArgDraw;
    /// \endcond

    
    /// \cond
    void *ThreadDraw8(void *arg){
      ArgDraw *p_arg;
      AdjRegion3 *adjreg;
      Scene8::Scene8 *scn;
      Voxel v,u;
      int i;
      uchar val;
      
      p_arg  = (ArgDraw *)arg;
      adjreg = p_arg->adjreg;
      scn    = p_arg->scn8;
      u      = p_arg->u;
      val    = p_arg->val;

      for(i=p_arg->i; i<=p_arg->j; i++){
	v.v = u.v + adjreg->d[i].v;
	if(Scene8::IsValidVoxel(scn,v))
	  scn->array[v.c.z][v.c.y][v.c.x] = val;
      }
      return NULL;
    }
    /// \endcond

    
    void  MT_Draw(AdjRegion3 *adjreg,
		  Scene8::Scene8 *scn,
		  Voxel u,
		  uchar val){
      ArgDraw arg;
      ArgDraw args[8];
      int i,nprocs;
      int first,last,nelems,de;
      pthread_t thread_id[8];
      int iret[8];
      
      nprocs = GetNumberOfProcessors();
      //printf("nprocs: %d\n",nprocs);
      
      if(nprocs<=1){ 
	Draw(adjreg,scn,u,val);
	return;
      }
      if(nprocs>=8) nprocs = 8;
      
      arg.adjreg = adjreg;
      arg.scn8   = scn;
      arg.u      = u;
      arg.p      = NIL;
      arg.val    = val;
      
      first  = 0;
      last   = adjreg->n-1;
      nelems = last-first+1;
      de     = nelems/nprocs;
      arg.i = NIL;
      arg.j  = first-1;
      for(i=0; i<nprocs; i++){
	args[i] = arg;
	
	args[i].i = arg.j+1;
	if(i<nprocs-1) args[i].j = args[i].i+(de-1);
	else           args[i].j = last;
	
	//Create independent threads each of which will execute function
	iret[i] = pthread_create(&thread_id[i], NULL, 
				 ThreadDraw8,
				 (void*)&args[i]);
	arg = args[i];
      }
 
      //Wait till threads are complete before main continues.
      for(i=0; i<nprocs; i++)
	pthread_join(thread_id[i], NULL);
    }
    

    //----------------------

   /// \cond
    void *ThreadDrawOpt8(void *arg){
      ArgDraw *p_arg;
      AdjRegion3 *adjreg;
      Scene8::Scene8 *scn;
      int i,p,q;
      uchar val;

      p_arg  = (ArgDraw *)arg;
      adjreg = p_arg->adjreg;
      scn    = p_arg->scn8;
      p      = p_arg->p;
      val    = p_arg->val;

      for(i=p_arg->i; i<=p_arg->j; i++){
	q = p + adjreg->dp[i];
	scn->data[q] = val;
      }
      return NULL;
    }
    /// \endcond


    void  MT_DrawOpt(AdjRegion3 *adjreg,
		     Scene8::Scene8 *scn,
		     int p, uchar val){
      ArgDraw arg;
      ArgDraw args[8];
      int i,nprocs;
      int first,last,nelems,de;
      pthread_t thread_id[8];
      int iret[8];

      Optimize(adjreg, scn);

      nprocs = GetNumberOfProcessors();
      //printf("nprocs: %d\n",nprocs);

      if(nprocs<=1){
	DrawOpt(adjreg,scn,p,val);
	return;
      }
      if(nprocs>=8) nprocs = 8;

      arg.adjreg = adjreg;
      arg.scn8  = scn;
      arg.p     = p;
      arg.val   = val;
      
      first  = 0;
      last   = adjreg->n-1;
      nelems = last-first+1;
      de     = nelems/nprocs;
      arg.i  = NIL;
      arg.j  = first-1;
      for(i=0; i<nprocs; i++){
	args[i] = arg;

	args[i].i = arg.j+1;
	if(i<nprocs-1) args[i].j = args[i].i+(de-1);
	else           args[i].j = last;

	//Create independent threads each of which will execute function
	iret[i] = pthread_create(&thread_id[i], NULL, 
				 ThreadDrawOpt8,
				 (void*)&args[i]);
	arg = args[i];
      }

      //Wait till threads are complete before main continues.
      for(i=0; i<nprocs; i++)
	pthread_join(thread_id[i], NULL);
    }
    */

  } //end AdjRegion3 namespace
} //end gft namespace



