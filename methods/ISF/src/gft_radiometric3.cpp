
#include "gft_radiometric3.h"


namespace gft{

  namespace Scene32{

    Curve::Curve *Histogram(Scene32 *scn, int binwidth){
      int i,nbins,maxbins,b;
      Curve::Curve *hist = NULL;

      maxbins = GetMaximumValue(scn)+1;

      nbins = maxbins/binwidth;
      if(maxbins%binwidth!=0) nbins+=1;

      hist  = Curve::Create(nbins);

      for(i=0; i<scn->n; i++){
	b = scn->data[i]/binwidth;
	hist->Y[b]++;
      }
      for(i=0; i<nbins; i++)
	hist->X[i] = i*binwidth + (float)(binwidth-1)/2.0;
      return (hist);
    }


    void LinearStretchinplace(Scene32 *scn, 
			      int omin,int omax,
			      int nmin,int nmax){
      float tmp,d,dn;
      int   p;

      d  = (float)(omax - omin);
      dn = (float)(nmax - nmin);
      for(p=0; p<scn->n; p++){
	if(scn->data[p] < omin)
	  scn->data[p] = nmin;
	else if(scn->data[p] > omax)
	  scn->data[p] = nmax;
	else if( d <= 0.0 ) scn->data[p] = nmin;
	else{
	  tmp = ((float)(scn->data[p] - omin))/d;
	  scn->data[p] = nmin + ROUND(tmp*dn);
	}
      }
      scn->maxval = nmax;
    }


    Scene32 *LinearStretch(Scene32 *scn, 
			   int omin,int omax,
			   int nmin,int nmax){
      Scene32 *sscn;
      sscn = Clone(scn);
      LinearStretchinplace(sscn, omin,omax, nmin,nmax);
      return sscn;
    }


  } //end Scene32 namespace


  namespace Scene16{

    Curve::Curve *Histogram(Scene16 *scn, int binwidth){
      int i,nbins,maxbins,b;
      Curve::Curve *hist = NULL;
      
      maxbins = GetMaximumValue(scn)+1;
      
      nbins = maxbins/binwidth;
      if(maxbins%binwidth!=0) nbins+=1;
      
      hist  = Curve::Create(nbins);
      
      for(i=0; i<scn->n; i++){
	b = scn->data[i]/binwidth;
	hist->Y[b]++;
      }
      for(i=0; i<nbins; i++)
	hist->X[i] = i*binwidth + (float)(binwidth-1)/2.0;
      return (hist);
    }


    void LinearStretchinplace(Scene16 *scn, 
			      int omin,int omax,
			      int nmin,int nmax){
      float tmp,d,dn;
      int   p;
      
      d  = (float)(omax - omin);
      dn = (float)(nmax - nmin);
      for(p=0; p<scn->n; p++){
	if(scn->data[p] < omin)
	  scn->data[p] = nmin;
	else if(scn->data[p] > omax)
	  scn->data[p] = nmax;
	else if( d <= 0.0 ) scn->data[p] = nmin;
	else{
	  tmp = ((float)(scn->data[p] - omin))/d;
	  scn->data[p] = nmin + ROUND(tmp*dn);
	}
      }
      scn->maxval = nmax;
    }


    Scene16 *LinearStretch(Scene16 *scn, 
			   int omin,int omax,
			   int nmin,int nmax){
      Scene16 *sscn;
      sscn = Clone(scn);
      LinearStretchinplace(sscn, omin,omax, nmin,nmax);
      return sscn;
    }

  } //end Scene16 namespace


  namespace Scene8{

    Curve::Curve *Histogram(Scene8 *scn, int binwidth){
      int i,nbins,maxbins,b;
      Curve::Curve *hist = NULL;

      maxbins = GetMaximumValue(scn)+1;

      nbins = maxbins/binwidth;
      if(maxbins%binwidth!=0) nbins+=1;

      hist  = Curve::Create(nbins);

      for(i=0; i<scn->n; i++){
	b = scn->data[i]/binwidth;
	hist->Y[b]++;
      }
      for(i=0; i<nbins; i++)
	hist->X[i] = i*binwidth + (float)(binwidth-1)/2.0;
      return (hist);
    }

    void LinearStretchinplace(Scene8 *scn, 
			      int omin,int omax,
			      int nmin,int nmax){
      float tmp,d,dn;
      int   p;
      
      d  = (float)(omax - omin);
      dn = (float)(nmax - nmin);
      for(p=0; p<scn->n; p++){
	if(scn->data[p] < omin)
	  scn->data[p] = nmin;
	else if(scn->data[p] > omax)
	  scn->data[p] = nmax;
	else if( d <= 0.0 ) scn->data[p] = nmin;
	else{
	  tmp = ((float)(scn->data[p] - omin))/d;
	  scn->data[p] = nmin + ROUND(tmp*dn);
	}
      }
      scn->maxval = nmax;
    }


    Scene8 *LinearStretch(Scene8 *scn, 
			  int omin,int omax,
			  int nmin,int nmax){
      Scene8 *sscn;
      sscn = Clone(scn);
      LinearStretchinplace(sscn, omin,omax, nmin,nmax);
      return sscn;
    }

  } //end Scene8 namespace


  namespace Scene{

    Curve::Curve *Histogram(Scene *scn, int binwidth){
      switch(scn->nbits){
      case  8:
	return Scene8::Histogram(scn->ptr.scn8, binwidth);
      case 16:
	return Scene16::Histogram(scn->ptr.scn16, binwidth);
      case 32:
	return Scene32::Histogram(scn->ptr.scn32, binwidth);
      }
      return NULL;
    }


    void LinearStretchinplace(Scene *scn, 
			      int omin,int omax,
			      int nmin,int nmax){
      switch(scn->nbits){
      case  8:
	Scene8::LinearStretchinplace(scn->ptr.scn8,omin,omax,nmin,nmax);
	break;
      case 16:
	Scene16::LinearStretchinplace(scn->ptr.scn16,omin,omax,nmin,nmax);
	break;
      case 32:
	Scene32::LinearStretchinplace(scn->ptr.scn32,omin,omax,nmin,nmax);
	break;
      }
    }


    Scene *LinearStretch(Scene *scn, 
			 int omin,int omax,
			 int nmin,int nmax){
      Scene *sscn;
      sscn = Clone(scn);
      LinearStretchinplace(sscn, omin,omax, nmin,nmax);
      return sscn;
    }


  } //end Scene namespace


} //end gft namespace




