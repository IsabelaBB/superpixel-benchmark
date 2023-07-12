
#include "gft_features.h"


namespace gft{

    namespace Features{


      Features *Create(int n, int nfeats, int style){
	Features *f = NULL;
	int i;
	f = (Features *) calloc(1,sizeof(Features));
	if(f == NULL)
	  gft::Error((char *)MSG1,(char *)"Features::Create");

	f->n = n;
	f->nfeats = nfeats;

	f->label = NULL;
	if(style & gftFeature_ALLOCLABEL)
	  f->label = gft::AllocIntArray(n);

	f->index = NULL;
	if(style & gftFeature_ALLOCINDEX)
	  f->index = gft::AllocIntArray(n);

	f->fv = NULL;
	if(style & gftFeature_ALLOCFV){
	  f->fv = (float **)calloc(n, sizeof(float *));
	  if(f->fv == NULL)
	    gft::Error((char *)MSG1,(char *)"Features::Create");
	  for(i = 0; i < n; i++)
	    f->fv[i] = gft::AllocFloatArray(nfeats);
	}

	f->dist = NULL;
	if(style & gftFeature_ALLOCDIST){
	  f->dist = (float **)calloc(n, sizeof(float *));
	  if(f->dist == NULL)
	    gft::Error((char *)MSG1,(char *)"Features::Create");
	  for(i = 0; i < n; i++)
	    f->dist[i] = gft::AllocFloatArray(n);
	}
	
	return f;
      }

      
      Features *Clone(Features *f){
	Features *c;
	int i, style = 0;
	if(f->label != NULL) style = style | gftFeature_ALLOCLABEL;
	if(f->index != NULL) style = style | gftFeature_ALLOCINDEX;
	if(f->fv != NULL)    style = style | gftFeature_ALLOCFV;
	if(f->dist != NULL)  style = style | gftFeature_ALLOCDIST;
	c = Create(f->n, f->nfeats, style);

	if(style & gftFeature_ALLOCLABEL)
	  memcpy(c->label, f->label, f->n*sizeof(int));
	if(style & gftFeature_ALLOCINDEX)	
	  memcpy(c->index, f->index, f->n*sizeof(int));

	for(i = 0; i < f->n; i++){
	  if(style & gftFeature_ALLOCFV)
	    memcpy(c->fv[i], f->fv[i], f->nfeats*sizeof(float));
	  if(style & gftFeature_ALLOCDIST)
	    memcpy(c->dist[i], f->dist[i], f->n*sizeof(float));
	}
	return c;
      }


      void     *Destroy(Features **f){
	Features *tmp;
	int i;
	tmp = *f;
	if(tmp != NULL){
	  if(tmp->label != NULL)
	    gft::FreeIntArray(&tmp->label);
	  
	  if(tmp->index != NULL)
	    gft::FreeIntArray(&tmp->index);
	  
	  if(tmp->fv != NULL){
	    for(i = 0; i < tmp->n; i++)
	      gft::FreeFloatArray(&tmp->fv[i]);
	    free(tmp->fv);
	  }
	  
	  if(tmp->dist != NULL){
	    for(i = 0; i < tmp->n; i++)
	      gft::FreeFloatArray(&tmp->dist[i]);
	    free(tmp->dist);
	  }
	  
	  free(tmp);
	}
	*f = NULL;
      }



      void Randomize(Features *f){
	int i,j,k,tmp;
	float ftmp;
	float *ptmp;
	
	srand((int)time(NULL));
	for(i = 0; i < f->n; i++){
	  j = gft::RandomInteger(0, f->n-1);

	  if(f->label != NULL){
	    tmp = f->label[i];
	    f->label[i] = f->label[j];
	    f->label[j] = tmp;
	  }

	  if(f->index != NULL){
	    tmp = f->index[i];
	    f->index[i] = f->index[j];
	    f->index[j] = tmp;
	  }

	  if(f->fv != NULL){
	    ptmp = f->fv[i];
	    f->fv[i] = f->fv[j];
	    f->fv[j] = ptmp;
	  }

	  if(f->dist != NULL){
	    ptmp = f->dist[i];
	    f->dist[i] = f->dist[j];
	    f->dist[j] = ptmp;

	    for(k = 0; k < f->n; k++){
	      ftmp = f->dist[k][i];
	      f->dist[k][i] = f->dist[k][j];
	      f->dist[k][j] = ftmp;
	    }
	  }
	  
	}
      }
      
      
      
      Features *RemoveSamples(Features **f,
			      float rate){
	Features *f0 = NULL,*f1 = NULL,*f2 = NULL;
	int *hist,*hist1,*hist2,*index,*pos;
	int Lmax,i,j,l,nsamples,nfeats,p1,p2,style = 0;

	f0 = *f;
	
	if(f0->label != NULL) style = style | gftFeature_ALLOCLABEL;
	if(f0->index != NULL) style = style | gftFeature_ALLOCINDEX;
	if(f0->fv != NULL   ) style = style | gftFeature_ALLOCFV;
	if(f0->dist != NULL ) style = style | gftFeature_ALLOCDIST;

	Lmax = 0;
	for(i = 0; i < f0->n; i++)
	  if(f0->label[i] > Lmax)
	    Lmax = f0->label[i];
	hist  = gft::AllocIntArray(Lmax+1);
	hist1 = gft::AllocIntArray(Lmax+1);
	hist2 = gft::AllocIntArray(Lmax+1);
	pos   = gft::AllocIntArray(Lmax+1);
	
	for(i = 0; i < f0->n; i++)
	  hist[f0->label[i]]++;

	nsamples = 0;
	for(i = 0; i <= Lmax; i++){
	  hist1[i] = ceil(hist[i]*rate);
	  hist2[i] = hist[i] - hist1[i];
	  nsamples += hist1[i];
	}

	pos[0] = 0;
	for(i = 1; i <= Lmax; i++)
	  pos[i] = pos[i-1]+hist[i-1];

	nfeats = f0->nfeats;
	f1 = Create(nsamples, nfeats, style);
       	f2 = Create(f0->n-nsamples, nfeats, style);

	//Index to positions ordered by label.
	index = gft::AllocIntArray(f0->n);
	for(i = 0; i < f0->n; i++){
	  l = f0->label[i];
	  j = pos[l];
	  index[j] = i;
	  pos[l]++;
	}
	
	p1 = p2 = 0;
	for(l = 0; l <= Lmax; l++){

	  for(i = 0; i < hist1[l]; i++){
	    if(f0->fv != NULL)
	      memcpy(f1->fv[p1],
		     f0->fv[index[p1+p2]],
		     nfeats*sizeof(float));
	    if(f0->label != NULL)
	      f1->label[p1] = f0->label[index[p1+p2]];
	    if(f0->index != NULL)
	      f1->index[p1] = f0->index[index[p1+p2]];
	    p1++;
	  }

	  for(j = 0; j < hist2[l]; j++){
	    if(f0->fv != NULL)
	      memcpy(f2->fv[p2],
		     f0->fv[index[p1+p2]],
		     nfeats*sizeof(float));
	    if(f0->label != NULL)
	      f2->label[p2] = f0->label[index[p1+p2]];
	    if(f0->index != NULL)
	      f2->index[p2] = f0->index[index[p1+p2]];
	    p2++;
	  }
	}
	
	//Changes lf e lf2.
	*f = f2;
	Destroy(&f0);
	
	gft::FreeIntArray(&index);
	gft::FreeIntArray(&hist);
	gft::FreeIntArray(&hist1);
	gft::FreeIntArray(&hist2);
	gft::FreeIntArray(&pos);
  
	return f1;
      }

      
      
      Features *Merge(Features *f1, 
		      Features *f2){
	Features *f;
	int i,n,nfeats,style = 0;
	
	if(f1==NULL && f2==NULL) return NULL;
	if(f1==NULL)
	  return Clone(f2);
	if(f2==NULL)
	  return Clone(f1);
	
	n = f1->n + f2->n;
	nfeats = f1->nfeats;
	if(nfeats != f2->nfeats) return NULL;

	if(f1->label != NULL || f2->label != NULL) style = style | gftFeature_ALLOCLABEL;
	if(f1->index != NULL || f2->index != NULL) style = style | gftFeature_ALLOCINDEX;
	if(f1->fv != NULL    || f2->fv != NULL)    style = style | gftFeature_ALLOCFV;
	if(f1->dist != NULL  || f2->dist != NULL)  style = style | gftFeature_ALLOCDIST;
	
	f = Create(n, nfeats, style);

	if(f1->label != NULL)
	  memcpy(f->label,           f1->label, f1->n*sizeof(int));
	if(f2->label != NULL)
	  memcpy(&(f->label[f1->n]), f2->label, f2->n*sizeof(int));

	if(f1->index != NULL)
	  memcpy(f->index,           f1->index, f1->n*sizeof(int));
	if(f2->index != NULL)
	  memcpy(&(f->index[f1->n]), f2->index, f2->n*sizeof(int));

	if(f1->fv != NULL)
	  for(i = 0; i < f1->n; i++)
	    memcpy(f->fv[i], f1->fv[i], nfeats*sizeof(float));

	if(f2->fv != NULL)
	  for(i = 0; i < f2->n; i++)
	    memcpy(f->fv[i+f1->n], f2->fv[i], nfeats*sizeof(float));

	if(f1->dist != NULL)
	  for(i = 0; i < f1->n; i++)
	    memcpy(f->dist[i], f1->dist[i], f1->n*sizeof(float));

	if(f2->dist != NULL)
	  for(i = 0; i < f2->n; i++)
	    memcpy(&(f->dist[i+f1->n][f1->n]), f2->dist[i], f2->n*sizeof(float));
	
	return f;
      }
      

      
    } //end Features namespace
} //end gft namespace


