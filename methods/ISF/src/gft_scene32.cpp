
#include "gft_scene32.h"

extern "C" {
  #include "gft_bzlib.h"
  #include "nifti1_io.h"
}


namespace gft{
  namespace Scene32{

  Scene32 *Create(int xsize,int ysize,int zsize){
    Scene32 *scn=NULL;
    int **tmp=NULL;
    int xysize;
    int i,j,p,N;

    scn = (Scene32 *) calloc(1,sizeof(Scene32));
    if(scn == NULL) 
      gft::Error((char *)MSG1,(char *)"Scene32::Create");
 
    scn->xsize  = xsize;
    scn->ysize  = ysize;
    scn->zsize  = zsize;
    scn->dx     = 1.0;
    scn->dy     = 1.0;
    scn->dz     = 1.0;
    scn->maxval = 0;
    scn->n      = xsize*ysize*zsize;
    scn->nii_hdr = NULL;

    //For SSE optimization we allocate a multiple of 4.
    N = scn->n;
    if(N%4!=0) N += (4-N%4);

    scn->data   = gft::AllocIntArray(N);
    if(scn->data==NULL) 
      gft::Error((char *)MSG1,(char *)"Scene32::Create");

    scn->array = (int ***) calloc(zsize, sizeof(int **));
    if(scn->array==NULL) 
      gft::Error((char *)MSG1,(char *)"Scene32::Create");

    tmp = (int **) calloc(zsize*ysize, sizeof(int *));
    if(tmp==NULL) 
      gft::Error((char *)MSG1,(char *)"Scene32::Create");

    scn->array[0] = tmp;
    for(i=1; i<zsize; i++)
      scn->array[i] = scn->array[i-1] + ysize;

    xysize = xsize*ysize;
    for(i=0; i<zsize; i++){
      for(j=0; j<ysize; j++){
	p = j + i*ysize;
	tmp[p] = scn->data + xysize*i + xsize*j;
      }}
    return(scn);
  }


  Scene32 *Create(Scene32 *scn){
    Scene32 *new_scn=NULL;
    new_scn = Create(scn->xsize,
		     scn->ysize,
		     scn->zsize);
    new_scn->dx = scn->dx;
    new_scn->dy = scn->dy;
    new_scn->dz = scn->dz;
    return new_scn;
  }


  void     Destroy(Scene32 **scn){
    Scene32 *aux;
    
    aux = *scn;
    if(aux != NULL){
      if(aux->data     != NULL) gft::FreeIntArray(&aux->data);
      if(aux->array[0] != NULL) free(aux->array[0]);
      if(aux->array    != NULL) free(aux->array);
      free(aux);    
      *scn = NULL;
    }
  }

  Scene32 *Clone(Scene32 *scn){
    Scene32 *aux = Create(scn->xsize,scn->ysize,scn->zsize);
    aux->dx     = scn->dx;
    aux->dy     = scn->dy;
    aux->dz     = scn->dz;
    aux->maxval = scn->maxval;
    memcpy(aux->data,scn->data,sizeof(int)*aux->n);
    return aux;
  }

  Scene32 *SubScene(Scene32 *scn, Voxel l, Voxel h){
    return SubScene(scn,
		    l.c.x, l.c.y, l.c.z,
		    h.c.x, h.c.y, h.c.z);
  }

  Scene32 *SubScene(Scene32 *scn,
		    int xl, int yl, int zl,
		    int xh, int yh, int zh){
    Scene32 *sub=NULL;
    Voxel v;
    int i,j;
    //v4si *p1,*p2;
  
    if(!IsValidVoxel(scn,xl,yl,zl)||
       !IsValidVoxel(scn,xh,yh,zh)||
       (xl > xh)||(yl>yh)||(zl>zh))
      return NULL;

    sub = Create(xh-xl+1,yh-yl+1,zh-zl+1);
    sub->dx = scn->dx;
    sub->dy = scn->dy;
    sub->dz = scn->dz;
    j = 0;
    for(v.c.z=zl; v.c.z<=zh; v.c.z++){
      for(v.c.y=yl; v.c.y<=yh; v.c.y++){
	/*
	for(v.c.x=xl; v.c.x<=xh-3; v.c.x+=4){
	  i = GetVoxelAddress(scn, v);
	  p1 = (v4si *)(&sub->data[j]);
	  p2 = (v4si *)(&scn->data[i]);
	  *p1 = *p2;
	  j+=4;
	  }
	for(; v.c.x <= xh; v.c.x++){
	  i = GetVoxelAddress(scn, v);
	  sub->data[j] = scn->data[i];
	  j++;
	}
	*/
	v.c.x=xl;
	i = GetVoxelAddress(scn, v);
	memcpy(&sub->data[j],&scn->data[i],sizeof(int)*sub->xsize);
	j += sub->xsize;
      }
    }
    return(sub);
  }


  void     Copy(Scene32 *dest, Scene32 *src){
    if(dest->xsize!=src->xsize ||
       dest->ysize!=src->ysize ||
       dest->zsize!=src->zsize)
      gft::Error((char *)"Incompatible data",(char *)"Scene32::Copy");
    dest->dx     = src->dx;
    dest->dy     = src->dy;
    dest->dz     = src->dz;
    dest->maxval = src->maxval;
    memcpy(dest->data,src->data,sizeof(int)*src->n);
  }


  void     Copy(Scene32 *dest, Scene32 *src, Voxel v){
    int x1,y1,z1,x2,y2,z2,dx;
    Voxel t,u;
    int p,q;
    x1 = MAX(0, -v.c.x);
    x2 = MIN(src->xsize-1, dest->xsize-1 - v.c.x);
    dx = x2 - x1 + 1;
    if(dx<=0) return;

    y1 = MAX(0, -v.c.y);
    y2 = MIN(src->ysize-1, dest->ysize-1 - v.c.y);

    z1 = MAX(0, -v.c.z);
    z2 = MIN(src->zsize-1, dest->zsize-1 - v.c.z);
    
    for(t.c.z=z1; t.c.z<=z2; t.c.z++){
      for(t.c.y=y1; t.c.y<=y2; t.c.y++){
	t.c.x=x1;
	p = GetVoxelAddress(src, t);

	u.v = v.v + t.v;
	q = GetVoxelAddress(dest, u);

	memcpy(&dest->data[q],&src->data[p],sizeof(int)*dx);
      }
    }
  }


  void     Fill(Scene32 *scn, int value){
    if(value==0)
      memset((void *)scn->data, 0, sizeof(int)*scn->n);
    else{
      int p;
      v4si v;
      v4si *ptr;
      ((int *)(&v))[0] = value;
      ((int *)(&v))[1] = value;
      ((int *)(&v))[2] = value;
      ((int *)(&v))[3] = value;
      for(p=0; p<scn->n; p+=4){
	ptr = (v4si *)(&scn->data[p]);
	*ptr = v;
      }
    }
    scn->maxval = value;
  }


  Scene32 *Read(char *filename){
    Scene32  *scn=NULL;	
    FILE   *fp=NULL;
    uchar  *data8=NULL;
#if _WIN32 || __BYTE_ORDER==__LITTLE_ENDIAN 
    ushort *data16=NULL;
#endif
    char    type[10];
    int     i,n,v,xsize,ysize,zsize;
    long    pos;
    int     null;
    
    // Checking file type
    int len = strlen(filename);
    if ( (len>=4) && ((strcasecmp(filename + len - 4, ".hdr")==0) || (strcasecmp(filename + len - 4, ".img")==0))) 
      /*return ReadScene_Analyze(filename); */
      return ReadNifti1(filename);
    if ( (len>=8) && (strcasecmp(filename + len - 8, ".scn.bz2")==0))
      return ReadCompressed(filename);
    if ( (len>=4) && (strcasecmp(filename + len - 4, ".nii")==0))
      return ReadNifti1(filename);
    if ( (len>=7) && (strcasecmp(filename + len - 7, ".nii.gz")==0))
      return ReadNifti1(filename);
    if ( (len<=4) || (strcasecmp(filename + len - 4, ".scn")!=0)) {
      gft::Error((char *)MSG2,(char *)"Scene32::Read: Invalid file name or extension.");
      return NULL;
    }
    
    // Read the scn file
    fp = fopen(filename,"rb");
    if (fp == NULL){
      gft::Error((char *)MSG2,(char *)"Scene32::Read");
    }
    null = fscanf(fp,"%s\n",type);
    if((strcmp(type,"SCN")==0)){
      null = fscanf(fp,"%d %d %d\n",&xsize,&ysize,&zsize);
      scn = Create(xsize,ysize,zsize);
      n = xsize*ysize*zsize;
      null = fscanf(fp,"%f %f %f\n",&scn->dx,&scn->dy,&scn->dz);
      null = fscanf(fp,"%d",&v);
      pos = ftell(fp);
      //printf(" current relative position in file: %ld\n",pos);///	
      fseek(fp,(pos+1)*sizeof(char),SEEK_SET); // +1 for the EOL \n character not included in last fscanf() call
      if (v==8){
	data8  = gft::AllocUCharArray(n);			
	null = fread(data8,sizeof(uchar),n,fp);
	for (i=0; i < n; i++) 
	  scn->data[i] = (int) data8[i];
	gft::FreeUCharArray(&data8);
      } else if (v==16) {
	// assuming that data was written with LITTLE-ENDIAN Byte ordering, i.e. LSB...MSB
#if _WIN32 || __BYTE_ORDER==__LITTLE_ENDIAN
	// for PCs, Intel microprocessors (LITTLE-ENDIAN too) -> no change
	data16 = gft::AllocUShortArray(n);
	null = fread(data16,sizeof(ushort),n,fp);
	for (i=0; i < n; i++)
	  scn->data[i] = (int) data16[i];
	gft::FreeUShortArray(&data16);
#else
	// for Motorola, IBM, SUN (BIG-ENDIAN) -> SWAp Bytes!
	gft::Warning("Data is converted from LITTLE to BIG-ENDIAN","Scene32::Read");
	data8 = gft::AllocUCharArray(2*n);
	null = fread(data8,sizeof(uchar),2*n,fp);
	j=0;
	for (i=0; i < 2*n; i+=2) {
	  scn->data[j] = (int) data8[i] + 256 * (int) data8[i+1];
	  j++;
	}
	gft::FreeUCharArray(&data8);
#endif
      } else { /* n = 32 */
	//Warning("32-bit data. Values may be wrong (little/big endian byte ordering)","ReadScene");
	n = xsize*ysize*zsize;
	null = fread(scn->data,sizeof(int),n,fp);      
      }
      fclose(fp);
    } else {
      fprintf(stderr,"Input scene must be SCN\n");
      exit(-1);
    }
    scn->maxval = GetMaximumValue(scn);
    return(scn);
  }


  Scene32 *ReadCompressed(char *filename) {
    FILE   *f;
    gft_BZFILE *b;
    int     nBuf;
    char   *buf;
    int     bzerror;
    Scene32 *scn = NULL;
    int     xsize, ysize, zsize, pos;
    int     nbits, hdr_size, i, p;
    float   dx, dy, dz;
    char    s[ 2000 ];
    
    f = fopen ( filename, "r" );
    buf = ( char* ) calloc( 1025, sizeof( char ) );
    if ( !f ) {
      printf("Erro ao abrir arquivos!\n");
      return NULL;
    }
    b = gft_BZ2_bzReadOpen ( &bzerror, f, 0, 0, NULL, 0 );
    if ( bzerror != gft_BZ_OK ) {
      gft_BZ2_bzReadClose ( &bzerror, b );
      fclose(f);
      printf("Erro ao abrir cena compactada!\n");
      return NULL;
    }
    
    nBuf = gft_BZ2_bzRead ( &bzerror, b, buf, 1024 );
    hdr_size = 0;
    if ( ( ( bzerror == gft_BZ_OK ) || ( bzerror == gft_BZ_STREAM_END ) ) && (nBuf > 0) ) {
      sscanf( buf, "%[^\n]\n%d %d %d\n%f %f %f\n%d\n", s, &xsize, &ysize, &zsize, &dx, &dy, &dz, &nbits );
      //printf( "s:%s, xsize:%d, ysize:%d, zsize:%d, dx:%f, dy:%f, dz:%f, bits:%d\n", s, xsize, ysize, zsize, dx, dy, dz, nbits );
      if ( strcmp( s, "SCN" ) != 0 ) {
	printf( "Format must by a compressed scene.\nFound %s\n", s );
	return 0;
      }
      scn = Create( xsize, ysize, zsize );
      scn->dx = dx;
      scn->dy = dy;
      scn->dz = dz;
      
      hdr_size = strlen( s ) + strlen( "\n" );
      for( i = 0; i < 3; i++ ) {
	sscanf( &buf[ hdr_size ], "%[^\n]\n", s );
	hdr_size += strlen( s ) + strlen( "\n" );
      }
    }
    else {
      printf("Erro ao ler cena compactada!\n");
      return 0;
    }
    
    gft_BZ2_bzReadClose ( &bzerror, b );
    rewind( f );
    b = gft_BZ2_bzReadOpen ( &bzerror, f, 0, 0, NULL, 0 );
    if ( bzerror != gft_BZ_OK ) {
      gft_BZ2_bzReadClose ( &bzerror, b );
      fclose(f);
      printf("Erro ao abrir cena compactada!\n");
      return NULL;
    }
    
    nBuf = gft_BZ2_bzRead ( &bzerror, b, buf, hdr_size );
    p = 0;
    do {
      nBuf = gft_BZ2_bzRead ( &bzerror, b, buf, 1024 );
      if ( ( ( bzerror == gft_BZ_OK ) || ( bzerror == gft_BZ_STREAM_END ) ) && ( nBuf > 0 ) ) {
	pos = 0;
	while( pos < nBuf - ( nbits / 8 ) + 1 ) {
	  if( nbits == 8 ) {
	    scn->data[ p ] = ( int ) ( ( uchar ) buf[ pos ] );
	    pos++;
	  }
	  else if( nbits == 16 ) {
	    scn->data[ p ] = ( int ) ( ( ( uchar ) buf[ pos + 1 ] ) * 256 + ( ( uchar ) buf[ pos ] ) );
	    if ( scn->data[ p ] < 0 ) printf("dado negativo.\n");
	    pos += 2;
	  }
	  else {
	    scn->data[ p ] = ( int ) ( ( uchar ) buf[ pos + 3 ] );
	    for( i = 2; i >= 0; i-- ) {
	      scn->data[ p ] = scn->data[ p ] * 256 + ( ( uchar ) buf[ pos + i ] );
	    }
	    pos +=4;
	  }
	  p++;
	}
      }
    } while( bzerror == gft_BZ_OK );
    //printf( "scn->n:%d, p:%d\n", scn->n, p );
    
    GetMaximumValue( scn );
    //printf( "max=%d\n", scn->maxval );
    gft_BZ2_bzReadClose ( &bzerror, b );
    fclose(f);
    free(buf);
    
    //printf( "fim!\n" );
    return( scn );
  }


  Scene32 *ReadNifti1(char *filename){
    nifti_image *nii;
    Scene32 *scn;
    float max;
    int p;
  
    nii = nifti_image_read( filename , 1 );
    scn = Create( nii->nx, nii->ny, nii->nz );
    scn->dx = nii->dx;
    scn->dy = nii->dy;
    scn->dz = nii->dz;
    scn->nii_hdr = nii;
    
    if( ( nii->datatype == NIFTI_TYPE_COMPLEX64 ) ||
	( nii->datatype == NIFTI_TYPE_FLOAT64 ) ||
	( nii->datatype == NIFTI_TYPE_RGB24 ) ||
	( nii->datatype == NIFTI_TYPE_RGB24 ) ||
	( nii->datatype >= NIFTI_TYPE_UINT32 ) ||
	( nii->dim[ 0 ] < 3 ) || ( nii->dim[ 0 ] > 4 ) ) {
      printf( "Error: Data format not supported, or header is corrupted.\n" );
      printf( "Data that is NOT supported: complex, double, RGB, unsigned integer of 32 bits, temporal series and statistical images.\n" );
      exit( -1 );
    }
  
    if( nii->datatype == NIFTI_TYPE_INT32 ) {
      //printf( "Integer src image.\n" );
      memcpy( scn->data, nii->data, nii->nvox * nii->nbyper );
    }
    if( ( nii->datatype == NIFTI_TYPE_INT16 ) || ( nii->datatype == NIFTI_TYPE_UINT16 ) ) {
      //printf( "Short src image.\n" );
      for( p = 0; p < scn->n; p++ ) {
	scn->data[ p ] = ( ( unsigned short* ) nii->data )[ p ];
      }
    }
    if( ( nii->datatype == NIFTI_TYPE_INT8 ) || ( nii->datatype == NIFTI_TYPE_UINT8 ) ) {
      //printf( "Char src image.\n" );
      for( p = 0; p < scn->n; p++ ) {
	scn->data[ p ] = ( ( unsigned char* ) nii->data )[ p ];
      }
    }
    if( nii->datatype == NIFTI_TYPE_FLOAT32 ) {
      //printf( "Float src image.\n" );
      // max used to set data to integer range without loosing its precision.
      max = 0.0;
      for( p = 0; p < scn->n; p++ ) {
	if( max < ( ( float* ) nii->data )[ p ] ) {
	  max = ( ( float* ) nii->data )[ p ];
	}
      }
      nii->flt_conv_factor = 10000.0 / max;
      for( p = 0; p < scn->n; p++ ) {
	scn->data[ p ] = ROUND( ( ( float* ) nii->data )[ p ] * nii->flt_conv_factor );
      }
    }
    
    free( nii->data );
    nii->data = NULL;
    return( scn );
  }


  void     Write(Scene32 *scn, char *filename){
    FILE *fp=NULL;
    int Imax;
    int i,n;
    uchar  *data8 =NULL;
    ushort *data16=NULL;

    // Checking file type
    int len = strlen(filename);
    if ( (len>=4) && ((strcasecmp(filename + len - 4, ".hdr")==0) || (strcasecmp(filename + len - 4, ".img")==0))) {
      //WriteScene_Analyze(scn, filename);
      WriteNifti1(scn, filename);
      return;
    }
    if ( (len>=8) && (strcasecmp(filename + len - 8, ".scn.bz2")==0)) {
      WriteCompressed(scn, filename);
      return;
    }
    if ( (len>=7) && (strcasecmp(filename + len - 7, ".nii.gz")==0)) {
      WriteNifti1( scn, filename );
      return;
    }
    if ( (len>=4) && (strcasecmp(filename + len - 4, ".nii")==0)) {
      WriteNifti1( scn, filename );
      return;
    }
    if ( (len<=4) || (strcasecmp(filename + len - 4, ".scn")!=0)) {
      gft::Error((char *)MSG2,(char *)"Scene32::Write: Invalid file name or extension.");
    }


    // Writing the scn file
    fp = fopen(filename,"wb"); 
    if(fp == NULL) 
      gft::Error((char *)MSG2,(char *)"Scene32::Write");

    fprintf(fp,"SCN\n");
    fprintf(fp,"%d %d %d\n",scn->xsize,scn->ysize,scn->zsize);
    fprintf(fp,"%f %f %f\n",scn->dx,scn->dy,scn->dz);
  
    Imax = GetMaximumValue(scn);
    
    n = scn->n;
    if(Imax < 256) {
      fprintf(fp,"%d\n",8);
      data8 = gft::AllocUCharArray(n);
      for(i=0; i<n; i++) 
	data8[i] = (uchar) scn->data[i];
      fwrite(data8,sizeof(uchar),n,fp);
      gft::FreeUCharArray(&data8);
    } else if(Imax < 65536) {
      fprintf(fp,"%d\n",16);
      data16 = gft::AllocUShortArray(n);
      for(i=0; i<n; i++)
	data16[i] = (ushort) scn->data[i];
      fwrite(data16,sizeof(ushort),n,fp);
      gft::FreeUShortArray(&data16);
    } else {
      fprintf(fp,"%d\n",32);
      fwrite(scn->data,sizeof(int),n,fp);
    }
    fclose(fp);
  }


  void WriteCompressed(Scene32 *scn, char *filename) {
    FILE   *f;
    //FILE   *in;
    gft_BZFILE *b;
    //int     nBuf = 1024;
    //char    buf[1024];
    int    Imax, n, i;
    int    bzerror;
    char   *data;
    uchar  *data8;
    ushort *data16;
    
    f = fopen( filename, "wb" );
    
    //if ( ( !f ) || ( !in ) ) {
    if ( !f ) {
      printf("Error on opening the compressed file %s\n", filename);
      return;
    }
    
    b = gft_BZ2_bzWriteOpen( &bzerror, f, 9, 0, 30 );
    if (bzerror != gft_BZ_OK) {
      gft_BZ2_bzWriteClose ( &bzerror, b, 0, 0, 0 );
      fclose(f);
      printf("Error on opening the file %s\n", filename);
    }
    
    //while ( (bzerror == gft_BZ_OK ) && (nBuf > 0) ) {
    if (bzerror == gft_BZ_OK) {
      /* get data to write into buf, and set nBuf appropriately */
      data = (char*) calloc(200, sizeof(int));
      n = scn->xsize*scn->ysize*scn->zsize;
      sprintf(data,"SCN\n");
      sprintf(data,"%s%d %d %d\n",data,scn->xsize,scn->ysize,scn->zsize);
      sprintf(data,"%s%f %f %f\n",data,scn->dx,scn->dy,scn->dz);
      Imax = GetMaximumValue(scn);
      if (Imax < 256) {
	//printf("8bits\n");
	sprintf(data,"%s%d\n",data,8);
	gft_BZ2_bzWrite ( &bzerror, b, data, strlen( data ) );
	data8 = gft::AllocUCharArray(n);
	for (i=0; i < n; i++) 
	  data8[i] = (uchar) scn->data[i];
	gft_BZ2_bzWrite ( &bzerror, b, data8, n );
	gft::FreeUCharArray(&data8);
      } else if (Imax < 65536) {
	//printf("16bits\n");
	sprintf(data,"%s%d\n",data,16);
	gft_BZ2_bzWrite ( &bzerror, b, data, strlen( data ) );
	data16 = gft::AllocUShortArray(n);
	for (i=0; i < n; i++)
	  data16[i] = (ushort) scn->data[i];
	gft_BZ2_bzWrite ( &bzerror, b, data16, 2 * n );
	gft::FreeUShortArray(&data16);
      } else {
	//printf("32bits\n");
	sprintf(data,"%s%d\n",data,32);
	gft_BZ2_bzWrite ( &bzerror, b, data, strlen( data ) );
	gft_BZ2_bzWrite ( &bzerror, b, scn->data, 4 * n );
      }
      free(data);
	/*
	  nBuf = fread(buf, sizeof(char), 1024, in);
	  if (nBuf > 0)
	  gft_BZ2_bzWrite ( &bzerror, b, buf, nBuf );
	  
	  if (bzerror == gft_BZ_IO_ERROR) { 
	  break;
	  }
	*/
    }
    
    gft_BZ2_bzWriteClose ( &bzerror, b, 0, 0, 0 );
    fclose(f);
    
    if (bzerror == gft_BZ_IO_ERROR) {
      printf("Error on writing to %s\n", filename);
    }
  }


  void WriteNifti1(Scene32 *scn, char *filename){
    nifti_image *nii;
    int len;
    int p;
    
    nii = scn->nii_hdr;
    if( nii == NULL )
      gft::Error((char *)MSG4, (char *)"Scene32::WriteNifti1" );
    
    if( nii->data != NULL )
      free( nii->data );
    nii->data = calloc( nii->nbyper, nii->nvox );
    
    if( nii->datatype == NIFTI_TYPE_INT32 ) {
      //printf( "Integer src image.\n" );
      memcpy( nii->data, scn->data, nii->nvox * nii->nbyper );
    }
    else if( ( nii->datatype == NIFTI_TYPE_INT16 ) || ( nii->datatype == NIFTI_TYPE_UINT16 ) ) {
      //printf( "Short src image.\n" );
      for( p = 0; p < scn->n; p++ ) {
	( ( unsigned short* ) nii->data )[ p ] = scn->data[ p ];
      }
    }
    else if( ( nii->datatype == NIFTI_TYPE_INT8 ) || ( nii->datatype == NIFTI_TYPE_INT8 ) || ( nii->datatype == NIFTI_TYPE_UINT8 ) ) {
      //printf( "Char src image.\n" );
      for( p = 0; p < scn->n; p++ ) {
	( ( unsigned char* ) nii->data )[ p ] = scn->data[ p ];
      }
    }
    else if( nii->datatype == NIFTI_TYPE_FLOAT32 ) {
      //printf( "Float src image.\n" );
      // max used to set data to integer range without loosing its precision.
      for( p = 0; p < scn->n; p++ ) {
	( ( float* ) nii->data )[ p ] = scn->data[ p ] / nii->flt_conv_factor;
      }
    }
    else 
      gft::Error((char *)MSG5, (char *)"Scene32::WriteNifti1" );
    
    len = strlen( filename );
    if( ( strcasecmp( filename + len - 4, ".hdr" ) == 0 ) || ( strcasecmp( filename + len - 4, ".img" ) == 0 ) )
      nii->nifti_type = NIFTI_FTYPE_NIFTI1_2;
    else
      nii->nifti_type = NIFTI_FTYPE_NIFTI1_1;
    nifti_set_filenames( nii, filename, 0, nii->byteorder );
    
    nifti_image_write( nii );
    free( nii->data );
    nii->data = NULL;
  }



  float GetValue_trilinear(Scene32 *scn, float x,float y, float z){
    int px,py,pz,i1,i2;
    float dx,dy,dz;
    int nbx=0,nby=0,nbz=0;
    float p1,p2,p3,p4,p5,p6,res;

    if(x>scn->xsize-1||y>scn->ysize-1||z>scn->zsize-1||x<0||y<0||z<0)
      gft::Error((char *)"Out-of-bounds",
		 (char *)"Scene32::GetValue_trilinear");

    px = (int)x;  dx=x-px;
    py = (int)y;  dy=y-py;
    pz = (int)z;  dz=z-pz;

    // If it's not on the border, it has a neighour ahead.
    if(px<scn->xsize-1) nbx=1; 
    if(py<scn->ysize-1) nby=1;
    if(pz<scn->zsize-1) nbz=1;

    // 1st: Interpolate in Z
    i1 = scn->data[GetVoxelAddress(scn,px,py,pz)];
    i2 = scn->data[GetVoxelAddress(scn,px,py,pz+nbz)];
    p1 = i1 + dz*(i2-i1);
    i1 = scn->data[GetVoxelAddress(scn,px+nbx,py,pz)];
    i2 = scn->data[GetVoxelAddress(scn,px+nbx,py,pz+nbz)];
    p2 = i1 + dz*(i2-i1);
    i1 = scn->data[GetVoxelAddress(scn,px,py+nby,pz)];
    i2 = scn->data[GetVoxelAddress(scn,px,py+nby,pz+nbz)];
    p3 = i1 + dz*(i2-i1);
    i1 = scn->data[GetVoxelAddress(scn,px+nbx,py+nby,pz)];
    i2 = scn->data[GetVoxelAddress(scn,px+nbx,py+nby,pz+nbz)];
    p4 = i1 + dz*(i2-i1);
    // 2nd: Interpolate in X
    p5 = p1 + dx*(p2-p1);
    p6 = p3 + dx*(p4-p3);
    // 3rd: Interpolate in Y
    res = p5 + dy*(p6-p5);
    return res;
  }


  // return the nearest voxel.
  int   GetValue_nn(Scene32 *scn, float x, float y, float z){
    Voxel v;
    if(x>scn->xsize-1||y>scn->ysize-1||z>scn->zsize-1||
       x<(float)0||y<(float)0||z<(float)0)
      gft::Error((char *)"Out-of-bounds",
		 (char *)"Scene32::GetValue_nn");

    v.c.x=(int)(x+0.5);
    v.c.y=(int)(y+0.5);
    v.c.z=(int)(z+0.5);
    return scn->array[v.c.z][v.c.y][v.c.x];
  }


  int          GetMaximumValue(Scene32 *scn){
    int p,Imax;
    
    Imax = INT_MIN;
    for(p=0; p<scn->n; p++) {
      if(scn->data[p] > Imax)
	Imax = scn->data[p];
    }
    scn->maxval = Imax;
    return(Imax); 
  }

  int          GetMinimumValue(Scene32 *scn){
    int p,Imin;
  
    Imin = INT_MAX; 
    for(p=0; p<scn->n; p++){
      if(scn->data[p] < Imin)
	Imin = scn->data[p];
    }    
    return(Imin); 
  }


  gft::Image32::Image32 *GetSliceX(Scene32 *scn, int x){
    gft::Image32::Image32 *img;
    int z,y;
    img = gft::Image32::Create(scn->ysize, scn->zsize);
    for(z = 0; z < scn->zsize; z++){
      for(y = 0; y < scn->ysize; y++){
	img->array[z][y] = scn->array[z][y][x];
      }
    }
    return img;
  }


  gft::Image32::Image32 *GetSliceY(Scene32 *scn, int y){
    gft::Image32::Image32 *img;
    int z,x;
    img = gft::Image32::Create(scn->xsize, scn->zsize);
    for(z = 0; z < scn->zsize; z++){
      for(x = 0; x < scn->xsize; x++){
	img->array[z][x] = scn->array[z][y][x];
      }
    }
    return img;
  }

    
  gft::Image32::Image32 *GetSliceZ(Scene32 *scn, int z){
    gft::Image32::Image32 *img;
    int *data;
    img = gft::Image32::Create(scn->xsize, scn->ysize);
    data = scn->data + z*img->n;
    memcpy(img->data, data, sizeof(int)*img->n);
    return img;
  }

  void PutSliceX(Scene32 *scn, Image32::Image32 *img, int x){
    int z,y;
    for(z = 0; z < scn->zsize; z++){
      for(y = 0; y < scn->ysize; y++){
	scn->array[z][y][x] = img->array[z][y];
      }
    }
  }
    

  void PutSliceY(Scene32 *scn, Image32::Image32 *img, int y){
    int z,x;
    for(z = 0; z < scn->zsize; z++){
      for(x = 0; x < scn->xsize; x++){
	scn->array[z][y][x] = img->array[z][x];
      }
    }
  }

  void PutSliceZ(Scene32 *scn, Image32::Image32 *img, int z){
    int *data;
    data = scn->data + z*img->n;
    memcpy(data, img->data, sizeof(int)*img->n);
  }
    

  void     MBB(Scene32 *scn, Voxel *l, Voxel *h){
    Voxel v;
    
    l->c.x  = scn->xsize-1;
    l->c.y  = scn->ysize-1;
    l->c.z  = scn->zsize-1;
    h->c.x = 0;
    h->c.y = 0;
    h->c.z = 0;
  
    for(v.c.z=0; v.c.z<scn->zsize; v.c.z++)
      for(v.c.y=0; v.c.y<scn->ysize; v.c.y++)
	for(v.c.x=0; v.c.x<scn->xsize; v.c.x++)    
	  if(scn->data[GetVoxelAddress(scn,v)] > 0){
	    if(v.c.x < l->c.x)
	      l->c.x = v.c.x;
	    if(v.c.y < l->c.y)
	      l->c.y = v.c.y;
	    if(v.c.z < l->c.z)
	      l->c.z = v.c.z;
	    if(v.c.x > h->c.x)
	      h->c.x = v.c.x;
	    if(v.c.y > h->c.y)
	      h->c.y = v.c.y;	
	    if(v.c.z > h->c.z)
	      h->c.z = v.c.z;	
	  }
  }


  Scene32 *MBB(Scene32 *scn){
    Voxel lower,higher;
    Scene32 *mbb=NULL;
    MBB(scn, &lower, &higher);
    mbb = SubScene(scn,
		   lower.c.x,  lower.c.y,  lower.c.z,
		   higher.c.x, higher.c.y, higher.c.z);
    return(mbb);
  }


  Scene32 *AddFrame(Scene32 *scn,  int sz, int value){
    Scene32 *fscn;
    int y, z,*dst,*src,nbytes,offset1, offset2;
    
    fscn = Create(scn->xsize+(2*sz),
		  scn->ysize+(2*sz), 
		  scn->zsize+(2*sz));
    fscn->dx = scn->dx;
    fscn->dy = scn->dy;
    fscn->dz = scn->dz;

    Fill(fscn,value);
    nbytes = sizeof(int)*scn->xsize;
  
    offset1 = 0;
    offset2 = GetVoxelAddress(fscn, sz, sz, sz);
    
    for(z=0; z<scn->zsize; z++){
      src = scn->data+offset1;
      dst = fscn->data+offset2;
      for(y=0; y<scn->ysize; y++){
	memcpy(dst,src,nbytes);
	src += scn->xsize;
	dst += fscn->xsize;
      }
      offset1 += scn->xsize*scn->ysize;
      offset2 += fscn->xsize*fscn->ysize;
    }
    return(fscn);
  }


  Scene32 *RemFrame(Scene32 *fscn, int sz){
    Scene32 *scn;
    int y,z,*dst,*src,nbytes,offset;
    
    scn = Create(fscn->xsize-(2*sz),
		 fscn->ysize-(2*sz),
		 fscn->zsize-(2*sz));
    scn->dx = fscn->dx;
    scn->dy = fscn->dy;
    scn->dz = fscn->dz;

    nbytes = sizeof(int)*scn->xsize;  
    offset = GetVoxelAddress(fscn, sz, sz, sz);

    src = fscn->data+offset;
    dst = scn->data;
    for(z=0; z<scn->zsize; z++,src+=2*sz*fscn->xsize) {
      for(y=0; y<scn->ysize; y++,src+=fscn->xsize,dst+=scn->xsize){
	memcpy(dst,src,nbytes);
      }
    }
    return(scn);
  }


  Scene16::Scene16 *ConvertTo16(Scene32 *scn){
    Scene16::Scene16 *scn16;
    int p;

    scn16 = Scene16::Create(scn->xsize,
			    scn->ysize,
			    scn->zsize);
    scn16->dx = scn->dx;
    scn16->dy = scn->dy;
    scn16->dz = scn->dz;
    for(p=0; p<scn->n; p++){
      scn16->data[p] = (ushort)scn->data[p];
    }
    return scn16;
  }

  Scene8::Scene8  *ConvertTo8(Scene32 *scn){
    Scene8::Scene8 *scn8;
    int p;

    scn8 = Scene8::Create(scn->xsize,
			  scn->ysize,
			  scn->zsize);
    scn8->dx = scn->dx;
    scn8->dy = scn->dy;
    scn8->dz = scn->dz;
    for(p=0; p<scn->n; p++){
      scn8->data[p] = (uchar)scn->data[p];
    }
    return scn8;
  }


  Scene64::Scene64 *ComputeIntegralScene(Scene32 *scn){
    Scene64::Scene64 *Iscn=NULL;
    int p,q,i,j,k;
    int xysize = scn->xsize*scn->ysize;
    long long sum;
    
    //iscn = CreateScene(scn->xsize, scn->ysize, scn->zsize);
    //SetVoxelSize(iscn, scn->dx, scn->dy, scn->dz);
    Iscn = Scene64::Create(scn->xsize, scn->ysize, scn->zsize);
    
    for(k = 0; k < scn->zsize; k++){
      for(i = 0; i < scn->ysize; i++){
	for(j = 0; j < scn->xsize; j++){
	  p = GetVoxelAddress(scn,j,i,k);
	  sum = scn->data[p];
	  
	  if(k-1>=0){
	    q = p-xysize;
	    sum += Iscn->data[q];
	  }
	  if(j-1>=0){
	    q = p-1;
	    sum += Iscn->data[q];
	  }
	  if(j-1>=0 && k-1>=0){
	    q = p-1-xysize;
	    sum -= Iscn->data[q];
	  }
	  if(i-1>=0){
	    q = p-scn->xsize;
	    sum += Iscn->data[q];
	  }
	  if(i-1>=0 && k-1>=0){
	    q = p-scn->xsize-xysize;
	    sum -= Iscn->data[q];
	  }
	  if(j-1>=0 && i-1>=0){
	    q = p-1-scn->xsize;
	    sum -= Iscn->data[q];
	  }
	  if(j-1>=0 && i-1>=0 && k-1>=0){
	    q = p-1-scn->xsize-xysize;
	    sum += Iscn->data[q];
	  }
	  Iscn->data[p] = sum;
	}
      }
    }
    return Iscn;
  }
    


  //The dimensions of the window (i.e., xsize,ysize,zsize) 
  //should be given in millimeters.
  void DrawWindow(Scene32 *scn,
		  int val,
		  float xsize, 
		  float ysize, 
		  float zsize,
		  Voxel u){
    Voxel v,lo,hi;
    int dx,dy,dz;
    int q;
    
    dx = ROUND(xsize/(scn->dx*2.0));
    dy = ROUND(ysize/(scn->dy*2.0));
    dz = ROUND(zsize/(scn->dz*2.0));
    
    //---------------------
    hi.c.x = MIN(u.c.x+dx, scn->xsize-1);
    hi.c.y = MIN(u.c.y+dy, scn->ysize-1);
    hi.c.z = MIN(u.c.z+dz, scn->zsize-1);
    
    lo.c.x = MAX(u.c.x-dx-1, 0);
    lo.c.y = MAX(u.c.y-dy-1, 0);
    lo.c.z = MAX(u.c.z-dz-1, 0);
    
    for(v.c.x = lo.c.x+1; v.c.x <= hi.c.x; v.c.x++){
      for(v.c.y = lo.c.y+1; v.c.y <= hi.c.y; v.c.y++){
	for(v.c.z = lo.c.z+1; v.c.z <= hi.c.z; v.c.z++){
	  q = GetVoxelAddress(scn,v.c.x,v.c.y,v.c.z);
	  scn->data[q] = val;
	}
      }
    }
  }


    //The dimensions of the window (i.e., xsize,ysize,zsize) 
    //should be given in millimeters.
    int WindowNvoxels(Scene32 *scn,
		      float xsize,
		      float ysize,
		      float zsize){
      int nw,dx,dy,dz;
      dx = ROUND(xsize/(scn->dx*2.0));
      dy = ROUND(ysize/(scn->dy*2.0));
      dz = ROUND(zsize/(scn->dz*2.0));
      nw = (2*dx + 1)*(2*dy + 1)*(2*dz + 1);
      return nw;
    }   
    

    
  } //end Scene32 namespace
} //end gft namespace

