/*****************************************************************************\
* ift.c
*
* AUTHOR    : Felipe Belem (Org.) and Alexandre Falcao et. al.
* DATE      : 2021-02-05
* LICENSE   : MIT License
* COPYRIGHT : Alexandre Xavier Falcao 2012-2021
* EMAIL     : felipe.belem@ic.unicamp.br
*             afalcao@ic.unicamp.br
\*****************************************************************************/
#include "ift.h"


// ---------- iftBasicDataTypes.c start 

void iftCopyVoxel(iftVoxel *src, iftVoxel *dst) 
{
    (*dst).x = (*src).x;
    (*dst).y = (*src).y;
    (*dst).z = (*src).z;
}

// ---------- iftBasicDataTypes.c end
// ---------- iftIntArray.c start

iftIntArray *iftCreateIntArray(long n) 
{
    iftIntArray *iarr = (iftIntArray*) iftAlloc(1, sizeof(iftIntArray));
    
    iarr->n = n;
    iarr->val = iftAllocIntArray(n);
    
    return iarr;
}

void iftDestroyIntArray(iftIntArray **iarr) 
{
    if (iarr != NULL && *iarr != NULL) {
        iftIntArray *iarr_aux = *iarr;
        
        if (iarr_aux->val != NULL)
            iftFree(iarr_aux->val);
        iftFree(iarr_aux);
        *iarr = NULL;
    }
}

void iftShuffleIntArray(int* array, int n) 
{
    // Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
    int j;
    for (int i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        j = rand() % (i+1);
        // Swap arr[i] with the element at random index
        iftSwap(array[i],array[j]);
    }
}

// ---------- iftIntArray.c end
// ---------- iftFloatArray.c start

iftFloatArray *iftCreateFloatArray(long n) 
{
    iftFloatArray *darr = (iftFloatArray *) iftAlloc(1, sizeof(iftFloatArray));
    
    darr->n   = n;
    darr->val = iftAllocFloatArray(n);
    
    return darr;
}

void iftDestroyFloatArray(iftFloatArray **darr) 
{
    if (darr != NULL && *darr != NULL) {
        iftFloatArray *darr_aux = *darr;
        
        if (darr_aux->val != NULL)
            iftFree(darr_aux->val);
        iftFree(darr_aux);
        *darr = NULL;
    }
}

// ---------- iftFloatArray.c end 
// ---------- iftCommon.c start 
float iftRandomUniform(float low, float high) {
    double d;
    d = ((double) rand()) / ((double) RAND_MAX);
    return low + d * (high - low);
}

int iftRandomInteger (int low, int high){
    int k;
    double d;

    d = (double) rand () / ((double) RAND_MAX + 0.5);
    k =  iftMin((int)(d * (high - low + 1.0) ) + low, high);

    return k;
}

int *iftRandomIntegers(int low, int high, int nelems) {
    char msg[512];

    if (low > high) {
        sprintf(msg, "Low is greater than High: (%d, %d)", low, high);
        iftError(msg, "iftRandomIntegers");
    }

    int total_of_elems = (high - low + 1);

    if (nelems > total_of_elems) {
        sprintf(msg, "Nelems = %d is greater than the total of integer number in the range: [%d, %d]",
                nelems, low, high);
        iftError(msg, "iftRandomIntegers");
    }

    int *selected = iftAllocIntArray(nelems);
    int *values   = iftAllocIntArray(total_of_elems);
    int *count    = iftAllocIntArray(total_of_elems);

    int t = 0;
    for (int i = low; i <= high; i++) {
        values[t] = i;
        count[t]  = 100;
        t++;
    }

    if (nelems == total_of_elems) {
        iftFree(count);
        iftFree(selected);

        return values;
    }

    // Randomly select samples
    t = 0;
    int roof = total_of_elems - 1;
    while (t < nelems) {
        int i = iftRandomInteger(0, roof);
        int v = values[i];

        if (count[i] == 0) {
            selected[t] = v;
            iftSwap(values[i], values[roof]);
            iftSwap(count[i], count[roof]);
            t++;
            roof--;
        } else {
            count[i]--;
        }
    }
    iftFree(values);
    iftFree(count);

    return selected;
}
void iftSkipComments(FILE *fp)
{
    //skip for comments
    while (fgetc(fp) == '#') {
        while (fgetc(fp) != '\n');
    }
    fseek(fp,-1,SEEK_CUR);
}

double iftLog(double val, double base) 
{
    return (log(val) / log(base));
}

long iftNormalizationValue(long maxval) 
{
    long norm_val = 1;
    
    if (maxval < 0)
        iftError("Input value %ld < 0", "iftNormalizationValue", maxval);
    else if (maxval <= 1)
        norm_val = 1;
    else if (maxval <= 255)
        norm_val = 255;
    else if (maxval <= 4095)
        norm_val = 4095;
    else if (maxval <= 65535)
        norm_val = 65535;
    else if (maxval <= 4294967295)
        norm_val = 4294967295;
    else iftError("Invalid maxval number %ld with number of bits > 32. It only supports values within [0, 2ˆn_bits -1], " \
                  "where n_bits in {1, 8, 12, 16, 32}", "iftNormalizationValue", maxval);
    
    return norm_val;
}

void iftRandomSeed(unsigned int seed)
{

    srand(seed);

}

timer *iftTic()
{
    timer *tic=NULL;
    tic = (timer *)iftAlloc(1, sizeof(timer));
    gettimeofday(tic,NULL);
    return(tic);
}

timer *iftToc()
{
    timer *toc=NULL;
    toc = (timer *)iftAlloc(1, sizeof(timer));
    gettimeofday(toc,NULL);
    return(toc);
}

float iftCompTime(timer *tic, timer *toc)
{
    float t=0.0;
    if ((tic!=NULL)&&(toc!=NULL)){
        t = (toc->tv_sec-tic->tv_sec)*1000.0 +
            (toc->tv_usec-tic->tv_usec)*0.001;
        iftFree(tic);iftFree(toc);
    }
    return(t);
}

// ---------- iftCommon.c end
// ---------- iftAdjacency.c start 

iftAdjRel  *iftCreateAdjRel(int n)
{
  iftAdjRel *A=(iftAdjRel *)iftAlloc(1,sizeof(iftAdjRel ));

  A->dx = (int *)iftAllocIntArray(n);
  A->dy = (int *)iftAllocIntArray(n);
  A->dz = (int *)iftAllocIntArray(n);
  A->n  = n;

  return(A);
}

void     iftDestroyAdjRel(iftAdjRel **A)
{
  iftAdjRel *aux = *A;

  if (aux != NULL){
    if (aux->dx != NULL) iftFree(aux->dx);
    if (aux->dy != NULL) iftFree(aux->dy);
    if (aux->dz != NULL) iftFree(aux->dz);
    iftFree(aux);
    *A = NULL;
  }
}

iftAdjRel *iftSpheric(float r) 
{
    int r0 = (int) r;
    float r2 = (int) (r*r + 0.5);

    int n = 0;
    for (int dz = -r0; dz <= r0; dz++)
        for (int dy = -r0; dy <= r0; dy++)
            for (int dx =- r0; dx <= r0; dx++)
                if ( ((dx*dx) + (dy*dy) + (dz*dz)) <= r2)
                    n++;


    iftAdjRel *A = iftCreateAdjRel(n);

    int i = 0;
    int i0 = 0;
    for (int dz = -r0; dz <= r0; dz++)
        for (int dy = -r0; dy <= r0; dy++)
            for (int dx = -r0; dx <= r0; dx++)
                if ( ((dx*dx) + (dy*dy) + (dz*dz)) <= r2) {
                    A->dx[i] = dx;
                    A->dy[i] = dy;
                    A->dz[i] = dz;
                
                if ((dx == 0) && (dy == 0) && (dz == 0))
                    i0 = i;
                i++;
            }

    // shift to right and place central voxel at first
    for (int i = i0; i > 0; i--) {
        int dx = A->dx[i];
        int dy = A->dy[i];
        int dz = A->dz[i];
        A->dx[i] = A->dx[i-1];
        A->dy[i] = A->dy[i-1];
        A->dz[i] = A->dz[i-1];
        A->dx[i-1] = dx;
        A->dy[i-1] = dy;
        A->dz[i-1] = dz;
    }


    // sort by radius, so the 6 closest neighbors will come first
    float *dr = iftAllocFloatArray(A->n);
    for (int i = 0; i < A->n; i++)
        dr[i] = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i];

    iftIntArray *idxs = iftIntRange(0, A->n-1, 1);
    iftFQuickSort(dr, idxs->val, 0, A->n-1, IFT_INCREASING);
    iftAdjRel *Asort = iftCreateAdjRel(A->n);

    for (int i = 0; i < A->n; i++) {
        int idx = idxs->val[i];
        Asort->dx[i] = A->dx[idx];
        Asort->dy[i] = A->dy[idx];
        Asort->dz[i] = A->dz[idx];
    }

    iftFree(dr);
    iftDestroyIntArray(&idxs);
    iftDestroyAdjRel(&A);

    return Asort;
}

iftAdjRel *iftCircular(float r)
{
    int r0 = (int) r;
    float r2 = (int) (r*r + 0.5);
    
    int n = 0;
    for (int dy = -r0; dy <= r0; dy++)
        for (int dx = -r0; dx <= r0; dx++)
            if (((dx*dx) + (dy*dy)) <= r2)
                n++;

    iftAdjRel *A = iftCreateAdjRel(n);
    int i = 0;
    int i0 = 0;
    for (int dy = -r0; dy <= r0; dy++)
        for (int dx = -r0; dx <= r0; dx++)
            if (((dx*dx) + (dy*dy)) <= r2) {
                A->dx[i] = dx;
                A->dy[i] = dy;
                A->dz[i] = 0;

                if ((dx==0) && (dy==0))
                    i0 = i;
                i++;
            }

    // shift to right and place central pixel at first
    for (int i = i0; i > 0; i--) {
        int dx     = A->dx[i];
        int dy     = A->dy[i];
        A->dx[i]   = A->dx[i-1];
        A->dy[i]   = A->dy[i-1];
        A->dx[i-1] = dx;
        A->dy[i-1] = dy;
    }


    // sort by radius, so the 4 closest neighbors will come first
    float *dr = iftAllocFloatArray(A->n);
    for (int i = 0; i < A->n; i++)
        dr[i] = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i];


    iftIntArray *idxs = iftIntRange(0, A->n-1, 1);
    iftFQuickSort(dr, idxs->val, 0, A->n-1, IFT_INCREASING);
    iftAdjRel *Asort = iftCreateAdjRel(A->n);

    for (int i = 0; i < A->n; i++) {
        int idx = idxs->val[i];
        Asort->dx[i] = A->dx[idx];
        Asort->dy[i] = A->dy[idx];
        Asort->dz[i] = A->dz[idx];
    }

    iftFree(dr);
    iftDestroyIntArray(&idxs);
    iftDestroyAdjRel(&A);

    return Asort;
}

iftAdjRel *iftCopyAdjacency(const iftAdjRel *A) 
{
  iftAdjRel *B = iftCreateAdjRel(A->n);
  int i;
  for (i=0; i < A->n; i++) {
    B->dx[i] = A->dx[i];
    B->dy[i] = A->dy[i];
    B->dz[i] = A->dz[i];
  }

  return(B);
}

inline iftVoxel iftGetAdjacentVoxel(const iftAdjRel *A, iftVoxel u, int adj)
{
  iftVoxel v;

  v.x = u.x + A->dx[adj];
  v.y = u.y + A->dy[adj];
  v.z = u.z + A->dz[adj];

  return(v);
}

// ---------- iftAdjacency.c end
// ---------- iftColor.c start 

iftColor iftRGBColor(int R, int G, int B)
{
    iftColor color;
    color.val[0] = R;
    color.val[1] = G;
    color.val[2] = B;

    return color;
}

iftColor iftRGBtoYCbCrBT2020(iftColor cin, const int rgbBitDepth, const int yCbCrBitDepth)
{
    int minLum, minChr, quantLum, quantChr;
    iftColor cout;

    switch (yCbCrBitDepth) {
        case 8:
            minLum = 16; // 16 * 2^(bitDepth-8)
            minChr = 128; // 128 * 2^(bitDepth-8)
            quantLum = 219.0; // 219 * 2^(bitDepth-8)
            quantChr = 224.0;  // 224 * 2^(bitDepth-8)
            break;
        case 10:
            minLum = 64; // 16 * 2^(bitDepth-8)
            minChr = 512; // 128 * 2^(bitDepth-8)
            quantLum = 876; // 219 * 2^(bitDepth-8)
            quantChr = 896;  // 224 * 2^(bitDepth-8)
            break;
        case 12:
            minLum = 256; // 16 * 2^(bitDepth-8)
            minChr = 2048; // 128 * 2^(bitDepth-8)
            quantLum = 3504; // 219 * 2^(bitDepth-8)
            quantChr = 3584;  // 224 * 2^(bitDepth-8)
            break;
        case 16:
            minLum = 4096; // 16 * 2^(bitDepth-8)
            minChr = 32768; // 128 * 2^(bitDepth-8)
            quantLum = 56064.0; // 219 * 2^(bitDepth-8)
            quantChr = 57344.0;  // 224 * 2^(bitDepth-8)
            break;
        default:
            iftError("Bit depth not specified in BT.2020", "iftRGBtoYCbCrBT2020");
            cout.val[0] = cout.val[1] = cout.val[2] = 0;
            return cout;
    }

    double maxRgbValue = (double) ((1 << rgbBitDepth) - 1);
    double r = cin.val[0] / maxRgbValue;
    double g = cin.val[1] / maxRgbValue;
    double b = cin.val[2] / maxRgbValue;

    double y = 0.2627 * r + 0.6780 * g + 0.0593 * b;
    double cb = (b - y) / 1.8814;
    double cr = (r - y) / 1.4746;

    // clip luminance to [0..1] and chrominance to [-0.5..0.5]
    if (y < 0.0) y = 0.0;
    else if (y > 1.0) y = 1.0;
    if (cb < -0.5) cb = -0.5;
    else if (cb > 0.5) cb = 0.5;
    if (cr < -0.5) cr = -0.5;
    else if (cr > 0.5) cr = 0.5;

    // perform quantization
    cout.val[0] = (int) (y * quantLum) + minLum;
    cout.val[1] = (int) (cb * quantChr) + minChr;
    cout.val[2] = (int) (cr * quantChr) + minChr;

    return cout;
}

iftColor iftRGBtoYCbCr(iftColor cin, int normalization_value)
{
    iftColor cout;
    float a = (16.0/255.0)*(float)normalization_value;
    float b = (128.0/255.0)*(float)normalization_value;

    cout.val[0]=(int)(0.256789062*(float)cin.val[0]+
                      0.504128906*(float)cin.val[1]+
                      0.09790625*(float)cin.val[2]+a);
    cout.val[1]=(int)(-0.148222656*(float)cin.val[0]+
                      -0.290992187*(float)cin.val[1]+
                      0.439214844*(float)cin.val[2]+b);
    cout.val[2]=(int)(0.439214844*(float)cin.val[0]+
                      -0.367789063*(float)cin.val[1]+
                      -0.071425781*(float)cin.val[2]+b);

    for(int i=0; i < 3; i++) {
        if (cout.val[i] < 0) cout.val[i] = 0;
        if (cout.val[i] > normalization_value) cout.val[i] = normalization_value;
    }

    return(cout);
}

iftColor iftYCbCrtoRGB(iftColor cin, int normalization_value)
{
    iftColor cout;
    float a = (16.0/255.0)*(float)normalization_value;
    float b = (128.0/255.0)*(float)normalization_value;

    cout.val[0]=(int)(1.164383562*((float)cin.val[0]-a)+
                      1.596026786*((float)cin.val[2]-b));

    cout.val[1]=(int)(1.164383562*((float)cin.val[0]-a)+
                      -0.39176229*((float)cin.val[1]-b)+
                      -0.812967647*((float)cin.val[2]-b));

    cout.val[2]=(int)(1.164383562*((float)cin.val[0]-a)+
                      2.017232143*((float)cin.val[1]-b));

    for(int i=0; i < 3; i++) {
        if (cout.val[i] < 0) cout.val[i] = 0;
        if (cout.val[i] > normalization_value) cout.val[i] = normalization_value;
    }

    return(cout);
}

iftColor iftYCbCrBT2020toRGB(iftColor cin, const int yCbCrBitDepth, const int rgbBitDepth)
{
    int minLum, minChr;
    double quantLum, quantChr;
    iftColor cout;

    switch (yCbCrBitDepth) {
        case 8:
            minLum = 16; // 16 * 2^(bitDepth-8)
            minChr = 128; // 128 * 2^(bitDepth-8)
            quantLum = 219.0; // 219 * 2^(bitDepth-8)
            quantChr = 224.0;  // 224 * 2^(bitDepth-8)
            break;
        case 10:
            minLum = 64; // 16 * 2^(bitDepth-8)
            minChr = 512; // 128 * 2^(bitDepth-8)
            quantLum = 876.0; // 219 * 2^(bitDepth-8)
            quantChr = 896.0;  // 224 * 2^(bitDepth-8)
            break;
        case 12:
            minLum = 256; // 16 * 2^(bitDepth-8)
            minChr = 2048; // 128 * 2^(bitDepth-8)
            quantLum = 3504.0; // 219 * 2^(bitDepth-8)
            quantChr = 3584.0;  // 224 * 2^(bitDepth-8)
            break;
        case 16:
            minLum = 4096; // 16 * 2^(bitDepth-8)
            minChr = 32768; // 128 * 2^(bitDepth-8)
            quantLum = 56064.0; // 219 * 2^(bitDepth-8)
            quantChr = 57344.0;  // 224 * 2^(bitDepth-8)
            break;
        default:
            iftError("Bit depth not specified in BT.2020", "iftYCbCrBT2020toRGB");
            cout.val[0] = cout.val[1] = cout.val[2] = 0;
            return cout;
    }

    double y = (cin.val[0] - minLum) / quantLum;
    double cb = (cin.val[1] - minChr) / quantChr;
    double cr = (cin.val[2] - minChr) / quantChr;

    double r = cr * 1.4746 + y;
    double b = cb * 1.8814 + y;
    double g = (y - 0.2627 * r - 0.0593 * b) / 0.6780;

    // clip rgb values to [0..1]
    if (r < 0.0) r = 0.0;
    else if (r > 1.0) r = 1.0;
    if (g < 0.0) g = 0.0;
    else if (g > 1.0) g = 1.0;
    if (b < 0.0) b = 0.0;
    else if (b > 1.0) b = 1.0;

    // perform quantization
    double maxRgbValue = (double) ((1 << rgbBitDepth) - 1);
    cout.val[0] = (int) (r * maxRgbValue);
    cout.val[1] = (int) (g * maxRgbValue);
    cout.val[2] = (int) (b * maxRgbValue);

    return cout;
}

iftFColor iftRGBtoLabNorm(iftColor rgb, int normalization_value)
{
    //RGB to XYZ

    float R = rgb.val[0]/(float)normalization_value;
    float G = rgb.val[1]/(float)normalization_value;
    float B = rgb.val[2]/(float)normalization_value;

    if(R <= 0.04045)	R = R/12.92;
    else	        R = pow((R+0.055)/1.055,2.4);

    if(G <= 0.04045)	G = G/12.92;
    else		G = pow((G+0.055)/1.055,2.4);

    if(B <= 0.04045)	B = B/12.92;
    else		B = pow((B+0.055)/1.055,2.4);

    float X = (0.4123955889674142161*R + 0.3575834307637148171*G + 0.1804926473817015735*B);
    float Y = (0.2125862307855955516*R + 0.7151703037034108499*G + 0.07220049864333622685*B);
    float Z = (0.01929721549174694484*R + 0.1191838645808485318*G + 0.9504971251315797660*B);

    //XYZ to lab
    X /= WHITEPOINT_X;
    Y /= WHITEPOINT_Y;
    Z /= WHITEPOINT_Z;
    X = LABF(X);
    Y = LABF(Y);
    Z = LABF(Z);
    float L = 116*Y - 16;
    float a = 500*(X - Y);
    float b = 200*(Y - Z);

    iftFColor lab;
    lab.val[0] = L;
    lab.val[1] = a;
    lab.val[2] = b;

    return lab;
}

iftColor iftRGBtoHSV(iftColor cin, int normalization_value) 
{
    float r = ((float)cin.val[0]/normalization_value),
            g = ((float)cin.val[1]/normalization_value),
            b = ((float)cin.val[2]/normalization_value), v, x, f;
    float a[3];
    int   i;
    iftColor cout;

    // RGB are each on [0, 1]. S and V are returned on [0, 1] and H is
    // returned on [0, 6].

    x = iftMin(iftMin(r, g), b);
    v = iftMax(iftMax(r, g), b);
    if (v == x) {
        a[0]=0.0;
        a[1]=0.0;
        a[2]=v;
    } else {
        f = (r == x) ? g - b : ((g == x) ? b - r : r - g);
        i = (r == x) ? 3 : ((g == x) ? 5 : 1);
        a[0]=((float)i)-f/(v-x);
        a[1]=(v-x)/v;
        a[2]=0.299*r+0.587*g+0.114*b;
    }

    // (un)normalize

    cout.val[0] = (int)(a[0]*60.0);
    cout.val[1] = (int)(a[1]*normalization_value);
    cout.val[2] = (int)(a[2]*normalization_value);

    return(cout);
}

iftColor iftHSVtoRGB(iftColor cin, int normalization_value) 
{
    // H is given on [0, 6]. S and V are given on [0, 1].
    // RGB are each returned on [0, 1].
    float h = ((float)cin.val[0]/60.0),
            s = ((float)cin.val[1]/normalization_value),
            v = ((float)cin.val[2]/normalization_value), m, n, f;
    float a[3]={0,0,0};
    int i;
    iftColor cout;

    if (s==0.0) {
        a[0]=a[1]=a[2]=v;
    } else {
        i = (int) floor(h);
        f = h - (float)i;
        if(!(i & 1)) f = 1 - f; // if i is even
        m = v * (1 - s);
        n = v * (1 - s * f);
        switch (i) {
            case 6:
            case 0: a[0]=v; a[1]=n; a[2]=m; break;
            case 1: a[0]=n; a[1]=v; a[2]=m; break;
            case 2: a[0]=m; a[1]=v; a[2]=n; break;
            case 3: a[0]=m; a[1]=n; a[2]=v; break;
            case 4: a[0]=n; a[1]=m; a[2]=v; break;
            case 5: a[0]=v; a[1]=m; a[2]=n; break;
        }
    }

    // (un)normalize
    for(i=0;i<3;i++)
        cout.val[i]=a[i]*normalization_value;

    return(cout);
}

// ---------- iftColor.c end
// ---------- iftDHeap.c start

iftDHeap *iftCreateDHeap(int n, double *value) 
{
    iftDHeap *H = NULL;
    int i;
    
    if (value == NULL) {
        iftError("Cannot create heap without priority value map", "iftCreateDHeap");
    }
    
    H = (iftDHeap *) iftAlloc(1, sizeof(iftDHeap));
    if (H != NULL) {
        H->n       = n;
        H->value   = value;
        H->color   = (char *) iftAlloc(sizeof(char), n);
        H->node    = (int *) iftAlloc(sizeof(int), n);
        H->pos     = (int *) iftAlloc(sizeof(int), n);
        H->last    = -1;
        H->removal_policy = MINVALUE;
        if (H->color == NULL || H->pos == NULL || H->node == NULL)
            iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateDHeap");
        for (i = 0; i < H->n; i++) {
            H->color[i] = IFT_WHITE;
            H->pos[i]   = -1;
            H->node[i] = -1;
        }
    }
    else
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateDHeap");
    
    return H;
}

void iftDestroyDHeap(iftDHeap **H) 
{
    iftDHeap *aux = *H;
    if (aux != NULL) {
        if (aux->node != NULL) iftFree(aux->node);
        if (aux->color != NULL) iftFree(aux->color);
        if (aux->pos != NULL)   iftFree(aux->pos);
        iftFree(aux);
        *H = NULL;
    }
}

char iftFullDHeap(iftDHeap *H) 
{
    if (H->last == (H->n - 1))
        return 1;
    else
        return 0;
}

char iftEmptyDHeap(iftDHeap *H) 
{
    if (H->last == -1){
        return 1;
    }else{
        return 0;
    }
}

char iftInsertDHeap(iftDHeap *H, int node) 
{
    
    if (!iftFullDHeap(H)) {
        H->last++;
        H->node[H->last] = node;
        H->color[node]   = IFT_GRAY;
        H->pos[node]     = H->last;
        iftGoUpDHeap(H, H->last);
        return 1;
    } else {
        iftWarning("DHeap is full","iftInsertDHeap");
        return 0;
    }
    
}

int iftRemoveDHeap(iftDHeap *H) 
{
    int node= IFT_NIL;
    
    if (!iftEmptyDHeap(H)) {
        node = H->node[0];
        H->pos[node]   = -1;
        H->color[node] = IFT_BLACK;
        H->node[0]     = H->node[H->last];
        H->pos[H->node[0]] = 0;
        H->node[H->last] = -1;
        H->last--;
        iftGoDownDHeap(H, 0);
    }else{
        iftWarning("DHeap is empty","iftRemoveDHeap");
    }
    
    return node;
    
}

void    iftRemoveDHeapElem(iftDHeap *H, int pixel)
{
    
    if(H->pos[pixel] == -1)
        iftError("Element is not in the Heap", "iftRemoveDHeapElem");
    
    double aux = H->value[pixel];
    
    if(H->removal_policy == MINVALUE)
        H->value[pixel] = IFT_INFINITY_DBL_NEG;
    else
        H->value[pixel] = IFT_INFINITY_DBL;
    
    iftGoUpDHeap(H, H->pos[pixel]);
    iftRemoveDHeap(H);
    
    H->value[pixel] = aux;
    H->color[pixel] = IFT_WHITE;
    
}

void  iftGoUpDHeap(iftDHeap *H, int i) 
{
    int j = iftDad(i);
    
    if(H->removal_policy == MINVALUE){
        
        while ((j >= 0) && (H->value[H->node[j]] > H->value[H->node[i]])) {
            iftSwap(H->node[j], H->node[i]);
            H->pos[H->node[i]] = i;
            H->pos[H->node[j]] = j;
            i = j;
            j = iftDad(i);
        }
    }
    else{ /* removal_policy == MAXVALUE */
        
        while ((j >= 0) && (H->value[H->node[j]] < H->value[H->node[i]])) {
            iftSwap(H->node[j], H->node[i]);
            H->pos[H->node[i]] = i;
            H->pos[H->node[j]] = j;
            i = j;
            j = iftDad(i);
        }
    }
}

void iftGoDownDHeap(iftDHeap *H, int i) 
{
    int j, left = iftLeftSon(i), right = iftRightSon(i);
    
    j = i;
    if(H->removal_policy == MINVALUE){
        
        if ((left <= H->last) &&
            (H->value[H->node[left]] < H->value[H->node[i]]))
            j = left;
        if ((right <= H->last) &&
            (H->value[H->node[right]] < H->value[H->node[j]]))
            j = right;
    }
    else{ /* removal_policy == MAXVALUE */
        
        if ((left <= H->last) &&
            (H->value[H->node[left]] > H->value[H->node[i]]))
            j = left;
        if ((right <= H->last) &&
            (H->value[H->node[right]] > H->value[H->node[j]]))
            j = right;
    }
    
    if(j != i) {
        iftSwap(H->node[j], H->node[i]);
        H->pos[H->node[i]] = i;
        H->pos[H->node[j]] = j;
        iftGoDownDHeap(H, j);
    }
}

void iftResetDHeap(iftDHeap *H)
{
    int i;
    
    for (i=0; i < H->n; i++) {
        H->color[i] = IFT_WHITE;
        H->pos[i]   = -1;
        H->node[i] = -1;
    }
    H->last = -1;
}

// ---------- iftDHeap.c end
// ---------- iftFile.c start

bool iftFileExists(const char *pathname) 
{
    return (iftPathnameExists(pathname) && !iftDirExists(pathname));
}

const char *iftFileExt(const char *pathname) 
{
    if (pathname == NULL)
        iftError("Pathname is NULL", "iftFileExt");
    
    const char *dot = strrchr(pathname, '.'); // returns a pointer to the last occurrence of '.'
    
    if ( (!dot) || (dot == pathname)) {
        return ("");
    } else {
        if (iftRegexMatch(pathname, "^.*\\.tar\\.(gz|bz|bz2)$") || iftRegexMatch(pathname, "^.*\\.(scn|nii)\\.gz$")) {
            dot -= 4; // points to the penultimate dot '.'
        }
        
        return dot; // returns the extension with '.'
    }
}

char *iftJoinPathnames(long n, ...) 
{
    if (n <= 0)
        iftError("Number of pathnames to be concatenated is <= 0", "iftJoinPathnames");
    
    long out_str_size = 1; // '\0'
    
    // Counts the size of the concatenated string
    va_list path_list;
    va_start(path_list, n);
    for (int i = 0; i < n; i++)
        out_str_size += strlen(va_arg(path_list, char*)) + 1; // one char for '/' (separation char)
    va_end(path_list);
    
    char *joined_path = iftAllocCharArray(out_str_size);
    char *aux = iftAllocCharArray(out_str_size);
    
    va_start(path_list, n);
    strcpy(joined_path, va_arg(path_list, char*));
    
    for (int i = 1; i < n; i++) {
        char *path = va_arg(path_list, char*);
        if (iftStartsWith(path, IFT_SEP_C))
            path++; // skip the first char, which is the directory separator
        
        if (iftEndsWith(joined_path, IFT_SEP_C))
            sprintf(aux, "%s%s", joined_path, path);
        else
            sprintf(aux, "%s%s%s", joined_path, IFT_SEP_C, path);

        iftFree(joined_path);
        joined_path = iftCopyString(aux);
    }
    iftFree(aux);
    
    return joined_path;
}

char *iftFilename(const char *pathname, const char *suffix) 
{
    if (pathname == NULL)
        iftError("Pathname is NULL", "iftFilename");
    
    char *base = iftSplitStringAt(pathname, IFT_SEP_C, -1);
    
    if ((suffix != NULL) && (!iftCompareStrings(suffix, ""))) {
        char *out_base = iftRemoveSuffix(base, suffix);
        iftFree(base);
        base = out_base;
    }
    
    return base;
}

void iftDestroyFile(iftFile **f) {
    if (f != NULL) {
        iftFile *f_aux = *f;
        
        if (f_aux != NULL) {
            if (f_aux->path != NULL) {
                iftFree(f_aux->path);
                f_aux->path = NULL;
            }
            if(f_aux->suffix != NULL) {
                iftFree(f_aux->suffix);
                f_aux->suffix = NULL;
            }
            iftFree(f_aux);
            *f = NULL;
        }
    }
}
bool iftIsImageFile(const char *img_pathname) {
    if (img_pathname == NULL)
        iftError("Image Pathname is NULL", "iftIsImageFile");
    
    return (iftFileExists(img_pathname) && iftIsImagePathnameValid(img_pathname));
}


bool iftIsImagePathnameValid(const char *img_pathname) {
    if (img_pathname == NULL)
        iftError("Image Pathname is NULL", "iftIsImagePathnameValid");
    
    char *lower_pathname = iftLowerString(img_pathname);
    bool is_valid        = iftRegexMatch(lower_pathname, "^.+\\.(jpg|jpeg|pgm|ppm|scn|png|scn\\.gz|zscn|hdr|nii|nii\\.gz)$");
    iftFree(lower_pathname);
    
    return is_valid;
}

iftFile *iftCreateFile(const char *format, ...) {
    va_list args;
    char pathname[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(pathname, format, args);
    va_end(args);
    
    // it is a Directory instead of a File
    if (iftDirExists(pathname))
        iftError("Pathname \"%s\" is a directory", "iftCreateFile", pathname);
    
    iftFile *f = (iftFile*) iftAlloc(1, sizeof(iftFile));
    f->path = iftCopyString(pathname);
    f->suffix = NULL;
    
    
    return f;
}

char *iftExpandUser(const char *path) {
    return iftReplaceString(path,"~", getenv("HOME"));
}

iftFile *iftCopyFile(const iftFile* file) {
    iftFile* f = NULL;
    
    if (file == NULL)
        iftError("The file to be copied is NULL", "iftCopyFile");
    else {
        f         = (iftFile*) iftAlloc(1, sizeof(iftFile));
        f->path   = iftCopyString(file->path);
        f->sample = file->sample;
        f->label  = file->label;
        f->status = file->status;
        f->suffix = NULL;
        if(file->suffix != NULL)
            f->suffix   = iftCopyString(file->suffix);
    }
    
    return f;
}

// ---------- iftFile.c end
// ---------- iftBMap.c start

iftBMap *iftCreateBMap(int n) 
{
    iftBMap *b;
    b= (iftBMap *) iftAlloc(1,sizeof(iftBMap));
    b->n        = n;
    b->nbytes   = n/8;
    if (n%8) b->nbytes++;
    b->val = (char *) iftAlloc(b->nbytes,sizeof(char));
    if (b->val==NULL){
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateBMap");
    }
    return b;
}

void iftDestroyBMap(iftBMap **bmap) 
{
    iftBMap *aux=*bmap;
    
    if (aux != NULL) {
        iftFree(aux->val);
        iftFree(aux);
        *bmap=NULL;
    }   
}
// ---------- iftBMap.c end
// ---------- iftImage.c start
#if IFT_LIBPNG
#include <png.h>
#endif
#if IFT_LIBJPEG
#include <jpeglib.h>
#endif

iftImage *iftReadImageByExt(const char *format, ...) 
{
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    if (!iftFileExists(filename))
        iftError("Image %s does not exist", "iftReadImageByExt", filename);

    iftImage *img = NULL;
    char *ext = iftLowerString(iftFileExt(filename));

    if(iftCompareStrings(ext, ".png")) {
        img = iftReadImagePNG(filename);
    }
    else if (iftCompareStrings(ext, ".pgm")){
        FILE *fp = fopen(filename,"r");
        char type[10];
        if(fscanf(fp,"%s",type)!=1) iftError("Reading Error", "iftReadImageByExt");
        if (iftCompareStrings(type,"P5")){
            fclose(fp);
            img   = iftReadImageP5(filename);
        } else {
            fclose(fp);
            img   = iftReadImageP2(filename);
        }
    } else if (iftCompareStrings(ext, ".ppm")){
        img   = iftReadImageP6(filename);
    } else if (iftCompareStrings(ext, ".scn")){
        img   = iftReadImage(filename);
    } else if (iftCompareStrings(ext, ".jpg") || iftCompareStrings(ext, ".jpeg")){
        img = iftReadImageJPEG(filename);
    } else {
        iftError("Invalid image format: \"%s\" - Try .scn, .ppm, .pgm, .jpg, .png",
                 "iftReadImageByExt", ext);
    }

    iftFree(ext);
    return(img);
}

iftImage  *iftCreateImage(int xsize,int ysize,int zsize) 
{
    int *val = iftAllocIntArray(xsize*ysize*zsize);

    return iftCreateImageFromBuffer(xsize, ysize, zsize, val);
}

iftImage *iftCopyImage(const iftImage *img) 
{
    if (img == NULL)
        return NULL;

    iftImage *imgc=iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftCopyImageInplace(img, imgc);

    return(imgc);
}

void iftDestroyImage(iftImage **img) 
{
    if(img != NULL) {
        iftImage *aux = *img;

        if (aux != NULL) {
            if (aux->val != NULL) iftFree(aux->val);
            if (aux->Cb != NULL) iftFree(aux->Cb);
            if (aux->Cr != NULL) iftFree(aux->Cr);
            if (aux->alpha != NULL) iftFree(aux->alpha);
            if (aux->tby != NULL) iftFree(aux->tby);
            if (aux->tbz != NULL) iftFree(aux->tbz);
            iftFree(aux);
            *img = NULL;
        }
    }
}

void iftWriteImageByExt(const iftImage *img, const char *format, ...) 
{
    if (img == NULL)
        iftWarning("Image is NULL... Nothing to write", "iftWriteImageByExt");
    else {
        char command[400];

        va_list args;
        char filename[IFT_STR_DEFAULT_SIZE];

        va_start(args, format);
        vsprintf(filename, format, args);
        va_end(args);

        char *parent_dir = iftParentDir(filename);
        if (!iftDirExists(parent_dir))
            iftMakeDir(parent_dir);
        iftFree(parent_dir);

        char *ext = iftLowerString(iftFileExt(filename));

        if(iftCompareStrings(ext, ".png")) {
            iftWriteImagePNG(img,filename);
        } else if (iftCompareStrings(ext, ".scn")) {
            iftWriteImage(img, filename);
        }else if (iftCompareStrings(ext, ".pgm")) {
            if (iftMaximumValue(img)>255)
                iftWriteImageP2(img,filename);
            else
                iftWriteImageP5(img,filename);
        } else if (iftCompareStrings(ext, ".ppm")){
            iftWriteImageP6(img,filename);
        } else if (iftIsColorImage(img)){
            iftWriteImageP6(img,"temp.ppm");
            sprintf(command,"convert temp.ppm %s",filename);
            if (system(command)==-1)
                iftError("Program convert failed or is not installed", "iftWriteImageByExt");
            if (system("rm -f temp.ppm")==-1)
                iftError("Cannot remore temp.ppm", "iftWriteImageByExt");
        } else if(iftCompareStrings(ext, ".jpg") || iftCompareStrings(ext, ".jpeg")) {
            iftWriteImageJPEG(img,filename);
        } else {
            printf("Invalid image format: %s. Please select among the accepted ones: .scn, .ppm, .pgm, .png\n",ext);
            exit(-1);
        }


        iftFree(ext);
    }
}

int iftMaximumValue(const iftImage *img) 
{
    if (img == NULL)
        iftError("Image is NULL", "iftMaximumValue");

    iftBoundingBox bb;
    bb.begin.x = bb.begin.y = bb.begin.z = 0;
    bb.end.x   = img->xsize-1;
    bb.end.y   = img->ysize-1;
    bb.end.z   = img->zsize-1;

    return iftMaximumValueInRegion(img, bb);
}

int iftMinimumValue(const iftImage *img) 
{
    int img_min_val = IFT_INFINITY_INT;

    for (int p = 0; p < img->n; p++)
        if (img_min_val > img->val[p])
            img_min_val = img->val[p];

    return img_min_val;
}

inline iftVoxel iftGetVoxelCoord(const iftImage *img, int p)
{
    /* old
     * u.x = (((p) % (((img)->xsize)*((img)->ysize))) % (img)->xsize)
     * u.y = (((p) % (((img)->xsize)*((img)->ysize))) / (img)->xsize)
     * u.z = ((p) / (((img)->xsize)*((img)->ysize)))
     */
    iftVoxel u;
    div_t res1 = div(p, img->xsize * img->ysize);
    div_t res2 = div(res1.rem, img->xsize);

    u.x = res2.rem;
    u.y = res2.quot;
    u.z = res1.quot;

    return u;
}

iftImage *iftSelectImageDomain(int xsize, int ysize, int zsize)
{
    iftImage *mask=iftCreateImage(xsize,ysize,zsize);

    iftSetImage(mask,1);
    return(mask);
}

iftBoundingBox iftMinBoundingBox(const iftImage *img, iftVoxel *gc_out) 
{
    if (img == NULL)
        iftError("Image is NULL", "iftMinBoundingBox");

    long n = 0; // number of spels non-background (non-zero)
    iftVoxel gc = {0.0, 0.0, 0.0};
    iftBoundingBox mbb;
    mbb.begin.x = mbb.begin.y = mbb.begin.z = IFT_INFINITY_INT;
    mbb.end.x = mbb.end.y = mbb.end.z = IFT_INFINITY_INT_NEG;

    for (long p = 0; p < img->n; p++) {
        if (img->val[p] != 0) {
            iftVoxel v = iftGetVoxelCoord(img, p);

            mbb.begin.x = iftMin(mbb.begin.x, v.x);
            mbb.begin.y = iftMin(mbb.begin.y, v.y);
            mbb.begin.z = iftMin(mbb.begin.z, v.z);

            mbb.end.x = iftMax(mbb.end.x, v.x);
            mbb.end.y = iftMax(mbb.end.y, v.y);
            mbb.end.z = iftMax(mbb.end.z, v.z);

            gc.x += v.x;
            gc.y += v.y;
            gc.z += v.z;
            n++;
        }
    }

    if (mbb.begin.x == IFT_INFINITY_INT) {
        mbb.begin.x = mbb.begin.y = mbb.begin.z = -1;
        mbb.end.x   = mbb.end.y   = mbb.end.z   = -1;
        gc.x        = gc.y        = gc.z        = -1.0;
    } else {
        gc.x /= n;
        gc.y /= n;
        gc.z /= n;
    }

    if (gc_out != NULL)
        *gc_out = gc;

    return mbb;
}

iftImage *iftReadImage(const char *format, ...) 
{
    iftImage *img    = NULL;
    FILE     *fp     = NULL;
    uchar    *data8  = NULL;
    ushort   *data16 = NULL;
    int      *data32 = NULL;
    char     type[10];
    int      p, v, xsize, ysize, zsize;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "rb");

    if (fp == NULL) {
        iftError("Cannot open file: \"%s\"", "iftReadImage", filename);
    }
    if (fscanf(fp, "%s\n", type) != 1)
        iftError("Reading error: Image type", "iftReadImage");

    if (iftCompareStrings(type, "SCN")) {

        //iftSkipComments(fp);

        if (fscanf(fp, "%d %d %d\n", &xsize, &ysize, &zsize) != 3)
            iftError("Reading error: Image resolution/size", "iftReadImage");
        img = iftCreateImage(xsize, ysize, zsize);
        if (fscanf(fp, "%f %f %f\n", &img->dx, &img->dy, &img->dz) != 3) {
            iftError("Reading error: Pixel/Voxel size", "iftReadImage");
        }
        if (fscanf(fp, "%d", &v) != 1)
            iftError("Reading error", "iftReadImage");

        while (fgetc(fp) != '\n');

        if (v == 8) {
            data8 = iftAllocUCharArray(img->n);
            if (fread(data8, sizeof(uchar), img->n, fp) != (uint)img->n)
                iftError("Reading error", "iftReadImage");
            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data8[p];
            iftFree(data8);
        } else if (v == 16) {
            data16 = iftAllocUShortArray(img->n);
            if (fread(data16, sizeof(ushort), img->n, fp) != (uint)img->n)
                iftError("Reading error 16 bits", "iftReadImage");
            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data16[p];
            iftFree(data16);
        } else if (v == 32) {
            data32 = iftAllocIntArray(img->n);
            if (fread(data32, sizeof(int), img->n, fp) != (uint)img->n)
                iftError("Reading error", "iftReadImage");
            for (p = 0; p < img->n; p++)
                img->val[p] = data32[p];
            iftFree(data32);
        } else {
            iftError("Input scene must be 8, 16, or 32 bit", "iftReadImage");
        }
    } else {
        iftError("Invalid file type", "iftReadImage");
    }

    fclose(fp);
    return (img);
}

#if IFT_LIBPNG
png_bytep* iftReadPngImageAux(const char *file_name, png_structp *png_ptr, png_infop *info_ptr)
{
    png_byte header[8];    // 8 is the maximum size that can be checked

    /* open file and test for it being a png */
    FILE *fp = fopen(file_name, "rb");
    if (!fp)
        iftError("File %s could not be opened for reading", "iftReadPngImageAux", file_name);
    if (fread(header, 1, 8, fp)!=8) iftError("Reading error", "iftReadPngImageAux");
    if (png_sig_cmp(header, 0, 8))
        iftError("File %s is not recognized as a PNG file", "iftReadPngImageAux", file_name);

    int height;
    png_bytep * row_pointers;

    /* initialize stuff */
    png_structp ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!ptr)
        iftError("Internal error: png_create_read_struct failed", "iftReadImagePNG");

    *png_ptr = ptr;
    *info_ptr = png_create_info_struct(*png_ptr);
    int depth = png_get_bit_depth((*png_ptr), (*info_ptr));
    if(depth < 8){
        png_set_expand_gray_1_2_4_to_8(ptr);
    }


    if (!(*info_ptr))
        iftError("Internal error: png_create_info_struct failed", "iftReadImagePNG");

    if (setjmp(png_jmpbuf(*png_ptr)))
        iftError("Internal error: Error during init_io", "iftReadImagePNG");


    png_init_io(*png_ptr, fp);
    png_set_sig_bytes(*png_ptr, 8);

    png_read_info(*png_ptr, *info_ptr);
	// reduces the pixels back down to the original bit depth
	//png_color_8p sig_bit = NULL;
	//if (png_get_sBIT(*png_ptr, *info_ptr, &sig_bit)) {
	//	png_set_shift(*png_ptr, sig_bit);
	//}
    height = png_get_image_height(*png_ptr, *info_ptr);
    png_read_update_info(*png_ptr, *info_ptr);


    /* read file */
    if (setjmp(png_jmpbuf(*png_ptr)))
        iftError("Internal error: Error during read_image", "iftReadImagePNG");

    row_pointers = (png_bytep*) iftAlloc(height, sizeof(png_bytep));
    for (int y=0; y<height; y++)
        row_pointers[y] = (png_byte*) iftAlloc(png_get_rowbytes(*png_ptr, *info_ptr), 1);


    png_read_image(*png_ptr, row_pointers);

    fclose(fp);

    return row_pointers;
}
#endif

iftImage* iftReadImagePNG(const char* format, ...) 
{
    #if IFT_LIBPNG
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    png_infop info_ptr;
    png_structp png_ptr;
    png_bytep *row_pointers;

    row_pointers = iftReadPngImageAux(filename, &png_ptr, &info_ptr);

    int width, height, color_type, depth;

    width = png_get_image_width(png_ptr, info_ptr);
    height = png_get_image_height(png_ptr, info_ptr);
    color_type = png_get_color_type(png_ptr, info_ptr);
    depth = png_get_bit_depth(png_ptr, info_ptr);
    iftImage* img = iftCreateImage(width, height, 1);
    unsigned int numberChannels = png_get_channels(png_ptr, info_ptr);

    int byteshift = depth/8;

    int x, y;

    int p = 0;

    if(color_type==PNG_COLOR_TYPE_GRAY)//gray image
    {
        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];
            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                img->val[p] = ptr[0];
                if(depth==16) {
                    img->val[p] = (img->val[p]<<8)+ptr[1];
                }
                p++;
            }
        }
    }else if(color_type==PNG_COLOR_TYPE_GRAY_ALPHA ){
        if(img->alpha == NULL){
            iftSetAlpha(img,0);
        }
        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];
            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                if(depth == 8){
                    img->val[p] = ptr[0];
                    img->alpha[p] = ptr[1];
                }
                else if(depth==16) {
                    img->val[p] = ptr[0];
                    img->val[p] = (img->val[p]<<8)+ptr[1];
                    img->alpha[p] = ptr[2];
                    img->alpha[p] = (img->alpha[p]<<8)+ptr[3];
                }
                p++;
            }
        }
    }
    else if(color_type == PNG_COLOR_TYPE_RGB){//color image

        iftSetCbCr(img, 128);
        iftColor rgb, ycbcr;

        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];

            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                rgb.val[0] = ptr[0*byteshift];
                rgb.val[1] = ptr[1*byteshift];
                rgb.val[2] = ptr[2*byteshift];

                if(depth==16) { //read second byte in case of 16bit images
                    rgb.val[0] = (rgb.val[0]<<8) + ptr[1];
                    rgb.val[1] = (rgb.val[1]<<8) + ptr[3];
                    rgb.val[2] = (rgb.val[2]<<8) + ptr[5];
                }

                ycbcr = iftRGBtoYCbCrBT2020(rgb, depth, depth);

                img->val[p] = ycbcr.val[0];
                img->Cb[p]  = ycbcr.val[1];
                img->Cr[p]  = ycbcr.val[2];

                p++;
            }
        }
    }else if(color_type == PNG_COLOR_TYPE_RGB_ALPHA){
        iftSetCbCr(img, 128);
        iftColor rgb, ycbcr;
        if(img->alpha == NULL){
            iftSetAlpha(img,0);
        }

        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];
            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                rgb.val[0] = ptr[0*byteshift];
                rgb.val[1] = ptr[1*byteshift];
                rgb.val[2] = ptr[2*byteshift];
                ushort alpha = ptr[3*byteshift];

                if(depth==16) { //read second byte in case of 16bit images
                    rgb.val[0] = (rgb.val[0]<<8) + ptr[1];
                    rgb.val[1] = (rgb.val[1]<<8) + ptr[3];
                    rgb.val[2] = (rgb.val[2]<<8) + ptr[5];
                    alpha = (alpha<<8) +  ptr[7];
                }

                ycbcr = iftRGBtoYCbCr(rgb, depth==8?255:65535);

                img->val[p] = ycbcr.val[0];
                img->Cb[p] = ycbcr.val[1];
                img->Cr[p] = ycbcr.val[2];
                img->alpha[p] = alpha;

                p++;
            }
        }
    }

    for (y = 0; y < height; ++y) {
        iftFree(row_pointers[y]);
    }

    iftFree(row_pointers);

    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

    img->dz = 0.0;

    return img;
    #else
    iftError("LibPNG support was not enabled!","iftReadImagePNG");
    return NULL;
    #endif
}

iftImage* iftReadImageJPEG(const char* format, ...) 
{
    #if IFT_LIBJPEG
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    iftImage* image = NULL;
    //code based on externals/libjpeg/source/example.c
    /* This struct contains the JPEG decompression parameters and pointers to
* working space (which is allocated as needed by the JPEG library).
*/
    struct jpeg_decompress_struct cinfo;
    /* We use our private extension JPEG error handler.
* Note that this struct must live as long as the main JPEG parameter
* struct, to avoid dangling-pointer problems.
*/
    struct jpeg_error_mgr jerr;

    /* More stuff */
    FILE * infile;		/* source file */
    JSAMPARRAY buffer;		/* Output row buffer */
    int row_stride;		/* physical row width in output buffer */
    /* In this example we want to open the input file before doing anything else,
 * so that the setjmp() error recovery below can assume the file is open.
 * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
 * requires it in order to read binary files.
 */

    if ((infile = fopen(filename, "rb")) == NULL) {
        printf("[readImageJPEG] can't open %s\n",filename);
        return NULL;
    }

    /* Step 1: allocate and initialize JPEG decompression object */

    /* We set up the normal JPEG error routines, then override error_exit. */
    cinfo.err = jpeg_std_error(&jerr);
    //jerr.pub.error_exit = my_error_exit;

    /* Establish the setjmp return context for my_error_exit to use. */
    jmp_buf setjmp_buffer;
    if (setjmp(setjmp_buffer)) {
        /* If we get here, the JPEG code has signaled an error.
         * We need to clean up the JPEG object, close the input file, and return.
         */
        jpeg_destroy_decompress(&cinfo);
        printf("[readImageJPEG] code has signaled an error\n");
        fclose(infile);
        return NULL;
    }

    /* Now we can initialize the JPEG decompression object. */
    jpeg_create_decompress(&cinfo);
    /* Step 2: specify data source (eg, a file) */
    jpeg_stdio_src(&cinfo, infile);
    /* Step 3: read file parameters with jpeg_read_header() */
    (void) jpeg_read_header(&cinfo, TRUE);
    /* We can ignore the return value from jpeg_read_header since
     *   (a) suspension is not possible with the stdio data source, and
     *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
     * See libjpeg.txt for more info.
     */

    /* Step 4: set parameters for decompression */

    /* In this example, we don't need to change any of the defaults set by
     * jpeg_read_header(), so we do nothing here.
     */

    /* Step 5: Start decompressor */
    (void) jpeg_start_decompress(&cinfo);
    /* We can ignore the return value since suspension is not possible
 * with the stdio data source.
 */

    /* We may need to do some setup of our own at this point before reading
 * the data.  After jpeg_start_decompress() we have the correct scaled
 * output image dimensions available, as well as the output colormap
 * if we asked for color quantization.
 * In this example, we need to make an output work buffer of the right size.
 */
    /* JSAMPLEs per row in output buffer */
    row_stride = cinfo.output_width * cinfo.output_components;
    /* Make a one-row-high sample array that will go away when done with image */
    buffer = (*cinfo.mem->alloc_sarray)
            ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


    /* Step 6: while (scan lines remain to be read) */
    /*           jpeg_read_scanlines(...); */

    /* Here we use the library's state variable cinfo.output_scanline as the
     * loop counter, so that we don't have to keep track ourselves.
     */
    image = iftCreateImage(cinfo.output_width,cinfo.output_height,1);

    //0 - JCS_GRAYSCALE
    //1 - JCS_RGB
    //2 - JCS_YCbCr
    //3 - JCS_CMYK
    //4 - JCS_YCCK
    //5 - JCS_BG_RGB
    //6 - JCS_BG_YCC
    unsigned  int imageRow = 0;
    unsigned int imageCol = 0;
    iftColor rgb;
    iftColor YCbCr;
    float scalingFactor = pow(2,cinfo.data_precision)-1;
    switch (cinfo.out_color_space){
        case JCS_GRAYSCALE:

            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    iftImgVal(image,imageCol,imageRow,0) = buffer[0][shift];
                    imageCol++;
                }
                imageRow++;
            }
            break;

        case JCS_RGB:

            iftSetCbCr(image,128);
            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    rgb.val[0] = buffer[0][shift+0];
                    rgb.val[1] = buffer[0][shift+1];
                    rgb.val[2] = buffer[0][shift+2];
                    YCbCr = iftRGBtoYCbCr(rgb,scalingFactor);
                    iftImgVal(image,imageCol,imageRow,0) = YCbCr.val[0];
                    iftImgCb(image,imageCol,imageRow,0) = YCbCr.val[1];
                    iftImgCr(image,imageCol,imageRow,0) = YCbCr.val[2];
                    imageCol++;
                }
                imageRow++;
            }

            break;

        case JCS_YCbCr:

            iftSetCbCr(image,128);
            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    iftImgVal(image,imageCol,imageRow,0) = buffer[0][shift+0];
                    iftImgCb(image,imageCol,imageRow,0) = buffer[0][shift+1];
                    iftImgCr(image,imageCol,imageRow,0) = buffer[0][shift+2];
                    imageCol++;
                }
                imageRow++;
            }

            break;
        case JCS_CMYK:

            iftSetCbCr(image,128);
            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    //convert CMYK to RGB (reference: http://www.rapidtables.com/convert/color/cmyk-to-rgb.htm)
                    rgb.val[0] = 255*(100-buffer[0][shift+0])*(100-buffer[0][shift+3]);
                    rgb.val[1] = 255*(100-buffer[0][shift+1])*(100-buffer[0][shift+3]);;
                    rgb.val[2] = 255*(100-buffer[0][shift+2])*(100-buffer[0][shift+3]);;
                    YCbCr = iftRGBtoYCbCr(rgb,scalingFactor);
                    iftImgVal(image,imageCol,imageRow,0) = YCbCr.val[0];
                    iftImgCb(image,imageCol,imageRow,0) = YCbCr.val[1];
                    iftImgCr(image,imageCol,imageRow,0) = YCbCr.val[2];
                    imageCol++;
                }
                imageRow++;
            }

            break;
        case JCS_YCCK:

            iftSetCbCr(image,128);
            imageRow = 0;
            iftWarning("Image is Y/Cb/Cr/K color space. The channel K is ignored", "iftReadImageJPEG");
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    iftImgVal(image,imageCol,imageRow,0) = buffer[0][shift+0];
                    iftImgCb(image,imageCol,imageRow,0) = buffer[0][shift+1];
                    iftImgCr(image,imageCol,imageRow,0) = buffer[0][shift+2];
                    imageCol++;
                }
                imageRow++;
            }

            break;
        //case JCS_BG_RGB:
        //    iftError("Big gamut red/green/blue color space not supported", "iftReadImageJPEG");
        //    break;
        //case JCS_BG_YCC:
        //    iftError("Big gamut red/green/blue color space not supported", "iftReadImageJPEG");
        //    break;
        default:
            iftError("Unkwon color space", "iftReadImageJPEG");
            break;
    }

    /* Step 7: Finish decompression */
    (void) jpeg_finish_decompress(&cinfo);

    /* This is an important step since it will release a good deal of memory. */
    jpeg_destroy_decompress(&cinfo);

    /* After finish_decompress, we can close the input file.
     * Here we postpone it until after no more JPEG errors are possible,
     * so as to simplify the setjmp error logic above.  (Actually, I don't
     * think that jpeg_destroy can do an error exit, but why assume anything...)
     */

    fclose(infile);
    /* At this point you may want to check to see whether any corrupt-data
     * warnings occurred (test whether jerr.pub.num_warnings is nonzero).
     */

    //jerr.num_warnings; //useful to know about corrupted data
    //printf("%ld\n",jerr.num_warnings);

    return image;
    #else
    iftError("LibJPEG support was not enabled!","iftReadImageJPEG");
    return NULL;
    #endif
}

iftImage *iftReadImageP5(const char *format, ...) 
{
    iftImage *img    = NULL;
    FILE     *fp     = NULL;
    uchar    *data8  = NULL;
    ushort   *data16 = NULL;
    char     type[10];
    int      p, v, xsize, ysize, zsize, hi, lo;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        iftError(MSG_FILE_OPEN_ERROR, "iftReadImageP5", filename);
    }

    if (fscanf(fp, "%s\n", type) != 1) {
        iftError("Reading error", "iftReadImageP5");
    }

    if (iftCompareStrings(type, "P5")) {

        iftSkipComments(fp);

        if (fscanf(fp, "%d %d\n", &xsize, &ysize) != 2)
            iftError("Reading error", "iftReadImageP5");
        zsize = 1;

        img = iftCreateImage(xsize, ysize, zsize);
        img->dz = 0.0;

        if (fscanf(fp, "%d", &v) != 1)
            iftError("Reading error", "iftReadImageP5");

        while (fgetc(fp) != '\n');

        if ((v <= 255) && (v > 0)) {
            data8 = iftAllocUCharArray(img->n);

            if (fread(data8, sizeof(uchar), img->n, fp) != (uint)img->n)
                iftError("Reading error", "iftReadImageP5");

            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data8[p];

            iftFree(data8);

        } else if ((v <= 65535) && (v > 255)) {
            data16 = iftAllocUShortArray(img->n);

            for (p = 0; p < img->n; p++) {
                if ((hi = fgetc(fp)) == EOF)
                    iftError("Reading error", "iftReadImageP5");
                if ((lo = fgetc(fp)) == EOF)
                    iftError("Reading error", "iftReadImageP5");

                data16[p] = (hi << 8) + lo;
            }

            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data16[p];

            iftFree(data16);

        } else {
            iftError("Invalid maximum value", "iftReadImageP5");
        }
    } else {
        iftError("Invalid image type", "iftReadImageP5");
    }

    fclose(fp);
    return (img);
}

iftImage *iftReadImageP6(const char *format, ...)
{
    iftImage  *img=NULL;
    FILE    *fp=NULL;
    char    type[10];
    int     p,v,xsize,ysize,zsize;
    ushort rgb16[3];
    iftColor RGB,YCbCr;

    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename,"r");
    if (fp == NULL){
        iftError(MSG_FILE_OPEN_ERROR, "iftReadImageP6", filename);
    }

    if(fscanf(fp,"%s\n",type)!=1)
        iftError("Reading error", "iftReadImageP6");
    if(iftCompareStrings(type,"P6")){

        iftSkipComments(fp);

        if(fscanf(fp,"%d %d\n",&xsize,&ysize)!=2)
            iftError("Reading error", "iftReadImageP6");

        zsize = 1;
        img = iftCreateImage(xsize,ysize,zsize);
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
        img->dz = 0.0;
        if (fscanf(fp,"%d",&v)!=1)
            iftError("Reading error", "iftReadImageP6");

        while(fgetc(fp) != '\n');

        if (v >= 0 && v < 256) {
            for (p=0; p < img->n; p++) {
                RGB.val[0] = fgetc(fp);
                RGB.val[1] = fgetc(fp);
                RGB.val[2] = fgetc(fp);
                YCbCr      = iftRGBtoYCbCr(RGB,255);
                img->val[p]=YCbCr.val[0];
                img->Cb[p] =(ushort)YCbCr.val[1];
                img->Cr[p] =(ushort)YCbCr.val[2];
            }
        } else if (v >= 256 && v <= 65536) {

            int rgbBitDepth = ceil(iftLog(v, 2));
            int ycbcrBitDepth = rgbBitDepth;

            if(ycbcrBitDepth<10)
                ycbcrBitDepth = 10;
            else if(ycbcrBitDepth < 12)
                ycbcrBitDepth = 12;
            else if(ycbcrBitDepth < 16)
                ycbcrBitDepth = 16;

            for (p=0; p < img->n; p++) {
                // read 6 bytes for each image pixel
                if (fread(rgb16, 2, 3, fp) == 3) {
                    // the PPM format specifies 2-byte integers as big endian,
                    // so we need to swap the bytes if the architecture is little endian
                    RGB.val[0]  = ((rgb16[0] & 0xff) << 8) | ((ushort) rgb16[0] >> 8);
                    RGB.val[1]  = ((rgb16[1] & 0xff) << 8) | ((ushort) rgb16[1] >> 8);
                    RGB.val[2]  = ((rgb16[2] & 0xff) << 8) | ((ushort) rgb16[2] >> 8);
                    YCbCr       = iftRGBtoYCbCrBT2020(RGB, rgbBitDepth, ycbcrBitDepth);
//                    YCbCr = iftRGBtoYCbCr(RGB, v);
                    img->val[p] = YCbCr.val[0];
                    img->Cb[p]  = (ushort)YCbCr.val[1];
                    img->Cr[p]  = (ushort)YCbCr.val[2];
                }
            }
        } else {
            iftError("Invalid maximum value", "iftReadImageP6");
        }
    }else{
        iftError("Invalid image type", "iftReadImageP6");
    }

    fclose(fp);
    return(img);
}

iftImage *iftReadImageP2(const char *format, ...) 
{
    iftImage *img = NULL;
    FILE     *fp  = NULL;
    char     type[10];
    int      p, v, xsize, ysize, zsize;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "r");

    if (fp == NULL) {
        iftError(MSG_FILE_OPEN_ERROR, "iftReadImageP2", filename);
    }

    if (fscanf(fp, "%s\n", type) != 1)
        iftError("Reading error", "iftReadImageP2");

    if (iftCompareStrings(type, "P2")) {

        iftSkipComments(fp);

        if (fscanf(fp, "%d %d\n", &xsize, &ysize) != 2)
            iftError("Reading error", "iftReadImageP2");
        zsize = 1;
        img   = iftCreateImage(xsize, ysize, zsize);
        img->dz = 0.0;
        if (fscanf(fp, "%d", &v) != 1)
            iftError("Reading error", "iftReadImageP2");

        while (fgetc(fp) != '\n');

        for (p = 0; p < img->n; p++)
            if (fscanf(fp, "%d", &img->val[p]) != 1)
                iftError("Reading error", "iftReadImageP2");

    } else {
        iftError("Invalid image type", "iftReadImageP2");
    }

    fclose(fp);
    return (img);
}

iftImage *iftCreateImageFromBuffer(int xsize, int ysize, int zsize, int *val) 
{
    iftImage *img = NULL;
    int      y, z, xysize;

    img = (iftImage *) iftAlloc(1, sizeof(iftImage));
    if (img == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateImage");
    }

    img->val   = val;
    img->Cb    = img->Cr = NULL;
    img->alpha = NULL;
    img->xsize = xsize;
    img->ysize = ysize;
    img->zsize = zsize;
    img->dx    = 1.0;
    img->dy    = 1.0;
    img->dz    = 1.0;
    img->tby   = iftAllocIntArray(ysize);
    img->tbz   = iftAllocIntArray(zsize);
    img->n     = xsize * ysize * zsize;

    if (img->val == NULL || img->tbz == NULL || img->tby == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateImage");
    }

    img->tby[0] = 0;
    for (y = 1; y < ysize; y++)
        img->tby[y] = img->tby[y - 1] + xsize;

    img->tbz[0] = 0;
    xysize = xsize * ysize;
    for (z = 1; z < zsize; z++)
        img->tbz[z] = img->tbz[z - 1] + xysize;

    return (img);
}

void iftCopyImageInplace(const iftImage *src, iftImage *dest) 
{
    int p;

    iftVerifyImageDomains(src, dest, "iftCopyImageInplace");

    iftCopyVoxelSize(src, dest);

    for (p=0; p < src->n; p++)
        dest->val[p]= src->val[p];

    if (src->Cb != NULL) {
        if(dest->Cb == NULL)
            dest->Cb = iftAllocUShortArray(src->n);
        if(dest->Cr == NULL)
            dest->Cr = iftAllocUShortArray(src->n);
        for (p=0; p < src->n; p++) {
            dest->Cb[p]= src->Cb[p];
            dest->Cr[p]= src->Cr[p];
        }
    }
}

void iftWriteImage(const iftImage *img, const char *format, ...) 
{
    FILE   *fp     = NULL;
    int    p;
    uchar  *data8  = NULL;
    ushort *data16 = NULL;
    int    *data32 = NULL;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    int img_min_val = iftMinimumValue(img);
    int img_max_val = iftMaximumValue(img);

    if (img_min_val < 0) {
        char msg[200];
        sprintf(msg, "Shifting image values from [%d,%d] to [%d,%d] on the original image\n",
                img_min_val, img_max_val, 0, img_max_val - img_min_val);
        iftWarning(msg, "iftWriteImage");
        for (p = 0; p < img->n; p++)
            img->val[p] = img->val[p] - img_min_val;
        img_max_val = img_max_val - img_min_val;
    }

    fp = fopen(filename, "wb");
    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftWriteImage", filename);

    fprintf(fp, "SCN\n");
    fprintf(fp, "%d %d %d\n", img->xsize, img->ysize, img->zsize);
    fprintf(fp, "%f %f %f\n", img->dx, img->dy, img->dz);


    if (img_max_val < 256) {
        fprintf(fp, "%d\n", 8);
        data8 = iftAllocUCharArray(img->n);
        for (p = 0; p < img->n; p++)
            data8[p] = (uchar) img->val[p];
        fwrite(data8, sizeof(uchar), img->n, fp);
        iftFree(data8);
    } else if (img_max_val < 65536) {
        fprintf(fp, "%d\n", 16);
        data16 = iftAllocUShortArray(img->n);
        for (p = 0; p < img->n; p++)
            data16[p] = (ushort) img->val[p];
        fwrite(data16, sizeof(ushort), img->n, fp);

        iftFree(data16);
    } else if (img_max_val < IFT_INFINITY_INT) {
        fprintf(fp, "%d\n", 32);
        data32 = iftAllocIntArray(img->n);
        for (p = 0; p < img->n; p++)
            data32[p] = img->val[p];
        fwrite(data32, sizeof(int), img->n, fp);
        iftFree(data32);
    }

    fclose(fp);
}

void iftWriteImageP5(const iftImage *img, const char *format, ...) 
{
    FILE   *fp     = NULL;
    int    p, hi, lo;
    uchar  *data8  = NULL;
    ushort *data16 = NULL;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "wb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteImageP5", filename);

    fprintf(fp, "P5\n");
    fprintf(fp, "%d %d\n", img->xsize, img->ysize);

    int img_max_val = iftMaximumValue(img);
    int img_min_val = iftMinimumValue(img);

    if ((img_max_val < 256) && (img_min_val >= 0)) {
        fprintf(fp, "%d\n", 255);
        data8 = iftAllocUCharArray(img->n);
        for (p = 0; p < img->n; p++)
            data8[p] = (uchar) img->val[p];
        fwrite(data8, sizeof(uchar), img->n, fp);
        iftFree(data8);
    } else if (img_max_val < 65536) {
        fprintf(fp, "%d\n", 65535);
        data16 = iftAllocUShortArray(img->n);
        for (p = 0; p < img->n; p++)
            data16[p] = (ushort) img->val[p];

        {
#define HI(num) (((num) & 0x0000FF00) >> 8)
#define LO(num) ((num) & 0x000000FF)
            for (p = 0; p < img->n; p++) {
                hi = HI(data16[p]);
                lo = LO(data16[p]);
                fputc(hi, fp);
                fputc(lo, fp);
            }
        }

        iftFree(data16);
    } else {
        char msg[200];
        sprintf(msg, "Cannot write image as P5 (%d/%d)", img_max_val, img_min_val);
        iftError(msg, "iftWriteImageP5");
    }
    fclose(fp);
}

void iftWriteImageP6(const iftImage *img, const char *format, ...) 
{
    FILE     *fp = NULL;
    int      p;
    ushort   rgb16[3];
    iftColor YCbCr, RGB;

    if (!iftIsColorImage(img))
        iftError("Image is not colored", "iftWriteImageP6");

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "w");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteImageP6", filename);

    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", img->xsize, img->ysize);

    int img_max_val = iftMaximumValue(img);
    int img_min_val = iftMinimumValue(img);

    if (img_min_val < 0) {
        iftError("Cannot write image as P6", "iftWriteImageP6");
    }
    if (img_max_val < 256) {
        fprintf(fp, "%d\n", 255);
        for (p = 0; p < img->n; p++) {
            YCbCr.val[0] = img->val[p];
            YCbCr.val[1] = img->Cb[p];
            YCbCr.val[2] = img->Cr[p];

            RGB = iftYCbCrtoRGB(YCbCr, 255);

            fputc(((uchar) RGB.val[0]), fp);
            fputc(((uchar) RGB.val[1]), fp);
            fputc(((uchar) RGB.val[2]), fp);
        }
    } else if (img_max_val < 65536) {
//        int rgbBitDepth = 9;
//        // find the bit depth for the maximum value img_max_val
//        while ((1 << rgbBitDepth) <= img_max_val) {
//            rgbBitDepth++;
//        }

        int rgbBitDepth = ceil(iftLog(img_max_val, 2));

        fprintf(fp, "%d\n", (1 << rgbBitDepth) - 1);
        for (p = 0; p < img->n; p++) {
            YCbCr.val[0] = img->val[p];
            YCbCr.val[1] = img->Cb[p];
            YCbCr.val[2] = img->Cr[p];
            RGB = iftYCbCrBT2020toRGB(YCbCr, rgbBitDepth, rgbBitDepth);
//            RGB = iftYCbCrtoRGB(YCbCr, img_max_val);
            // the PPM format specifies 2-byte integers as big endian,
            // so we need to swap the bytes if the architecture is little endian
            rgb16[0] = ((RGB.val[0] & 0xff) << 8) | ((ushort) RGB.val[0] >> 8);
            rgb16[1] = ((RGB.val[1] & 0xff) << 8) | ((ushort) RGB.val[1] >> 8);
            rgb16[2] = ((RGB.val[2] & 0xff) << 8) | ((ushort) RGB.val[2] >> 8);
            // write 6 bytes for each image pixel
            if (fwrite(rgb16, 2, 3, fp) != 3) {
                iftError("Cannot write 16-bit image as P6", "iftWriteImageP6");
            }
        }
    } else {
        iftError("Cannot write image as P6", "iftWriteImageP6");
    }
    fclose(fp);
}

void iftWriteImageP2(const iftImage *img, const char *format, ...) 
{
    FILE *fp = NULL;
    int  p;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];
    int     depth = iftImageDepth(img);
      
    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "w");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteImageP2", filename);

    fprintf(fp, "P2\n");
    fprintf(fp, "%d %d\n", img->xsize, img->ysize);

    int img_max_val = iftMaxImageRange(depth);
    fprintf(fp, "%d\n", img_max_val);
    for (p = 0; p < img->n; p++) {
        fprintf(fp, "%d ", img->val[p]);
        if (iftGetXCoord(img, p) == (img->xsize - 1)) fprintf(fp, "\n");
    }

    fclose(fp);
}

#if IFT_LIBPNG
void iftWritePngImageAux(const char *file_name, png_bytep *row_pointers, int width, int height, int bit_depth, int color_type) 
{

    /* create file */
    FILE *fp = fopen(file_name, "wb");
    if (!fp)
        iftError("Internal Error: File %s could not be opened for writing", "iftWritePngImageAux", file_name);


    /* initialize stuff */
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr)
        iftError("Internal Error: png_create_write_struct failed", "iftWriteImagePNG");

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
        iftError("Internal Error: png_create_info_struct failed", "iftWriteImagePNG");

    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during init_io", "iftWriteImagePNG");

    png_init_io(png_ptr, fp);


    /* write header */
    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during writing header", "iftWriteImagePNG");

    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);


    /* write bytes */
    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during writing bytes", "iftWriteImagePNG");

    png_write_image(png_ptr, row_pointers);


    /* end write */
    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during end of write", "iftWriteImagePNG");

    png_write_end(png_ptr, NULL);

    png_destroy_write_struct(&png_ptr, &info_ptr);

    /* cleanup heap allocation */
    for (int y=0; y<height; y++)
        iftFree(row_pointers[y]);
    iftFree(row_pointers);

    fclose(fp);
}
#endif

void iftWriteImagePNG(const iftImage* img, const char* format, ...) 
{
    #if IFT_LIBPNG
    png_bytep *row_pointers;
    int width, height, depth, byteshift;

    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    width = img->xsize;
    height = img->ysize;
    png_byte color_type;
    depth = iftImageDepth(img);

    if(depth<=8) {
        depth = 8;
    } else {
        depth = 16;
    }

    byteshift = depth/8;
    //int offset = depth==16?1:0;//to read the second byte first in cases of 16bit images

    size_t numberOfChannels=1;
    if(iftIsColorImage(img)){
        if(img->alpha == NULL){
            numberOfChannels = 3;//RGB
            color_type = PNG_COLOR_TYPE_RGB;
        }else{
            numberOfChannels = 4;//RGB_ALPHA
            color_type = PNG_COLOR_TYPE_RGB_ALPHA;
        }
    }else{
        if(img->alpha == NULL){
            numberOfChannels = 1;//GRAY
            color_type = PNG_COLOR_TYPE_GRAY;
        }else{
            numberOfChannels = 2;//GRAY_ALPHA
            color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
        }
    }

    //size_t pixel_size = (iftIsColorImage(img)?3:1 ) * byteshift;
    row_pointers = (png_bytep*) iftAlloc(height, sizeof(png_bytep));
    for (int y=0; y<height; y++)
        row_pointers[y] = (png_byte*) iftAlloc(width, numberOfChannels*byteshift);

    if(color_type == PNG_COLOR_TYPE_GRAY){
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                ptr[0] = img->val[p] & 0xFF;//get first byte

                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[1] = ptr[0];
                    ptr[0] = (img->val[p]>>8) & 0xFF;//get second byte
                }

                p++;
            }
        }
    }else if(color_type == PNG_COLOR_TYPE_GRAY_ALPHA){
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                if(depth==8){
                    ptr[0] = img->val[p] & 0xFF;//get first byte
                    ptr[1] = img->alpha[p] & 0xFF;//get second byte
                }


                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[0] = img->val[p]>>8;//get first byte
                    ptr[1] = img->val[p] & 0xFF;//get second byte


                    ptr[2] = img->alpha[p]>>8;//get first byte;
                    ptr[3] = img->alpha[p] & 0xFF;;//get second byte
                }
                p++;
            }
        }
    }else if(color_type == PNG_COLOR_TYPE_RGB){
        iftColor rgb, ycbcr;
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                ycbcr.val[0] = img->val[p];
                ycbcr.val[1] = img->Cb[p];
                ycbcr.val[2] = img->Cr[p];

                rgb = iftYCbCrBT2020toRGB(ycbcr, depth, depth);

                ptr[0*byteshift] = rgb.val[0] & 0xFF;//get first byte
                ptr[1*byteshift] = rgb.val[1] & 0xFF;
                ptr[2*byteshift] = rgb.val[2] & 0xFF;

                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[(0*byteshift)+1] = ptr[0*byteshift];
                    ptr[(1*byteshift)+1] = ptr[1*byteshift];
                    ptr[(2*byteshift)+1] = ptr[2*byteshift];

                    ptr[0*byteshift] = ((rgb.val[0]>>8) & 0xFF);//get second byte
                    ptr[1*byteshift] = ((rgb.val[1]>>8) & 0xFF);
                    ptr[2*byteshift] = ((rgb.val[2]>>8) & 0xFF);
                }

                p++;
            }
        }

    }else if(color_type == PNG_COLOR_TYPE_RGB_ALPHA){
        iftColor rgb, ycbcr;
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                ycbcr.val[0] = img->val[p];
                ycbcr.val[1] = img->Cb[p];
                ycbcr.val[2] = img->Cr[p];
                ushort alpha = img->alpha[p];

                rgb = iftYCbCrBT2020toRGB(ycbcr, depth, depth);

                ptr[0*byteshift] = rgb.val[0] & 0xFF;//get first byte
                ptr[1*byteshift] = rgb.val[1] & 0xFF;
                ptr[2*byteshift] = rgb.val[2] & 0xFF;
                ptr[3*byteshift] = alpha & 0xFF;

                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[(0*byteshift)+1] = ptr[0*byteshift];
                    ptr[(1*byteshift)+1] = ptr[1*byteshift];
                    ptr[(2*byteshift)+1] = ptr[2*byteshift];
                    ptr[(3*byteshift)+1] = ptr[(3*byteshift)];

                    ptr[0*byteshift] = ((rgb.val[0]>>8) & 0xFF);//get second byte
                    ptr[1*byteshift] = ((rgb.val[1]>>8) & 0xFF);
                    ptr[2*byteshift] = ((rgb.val[2]>>8) & 0xFF);
                    ptr[(3*byteshift)] = ((alpha>>8) & 0xFF);
                }
                p++;
            }
        }

    }else{
        iftError("Unknwon color scape", "iftWriteImagePNG");
    };


    iftWritePngImageAux(filename, row_pointers, width, height, depth, color_type);
    #else
    iftError("LibPNG support was not enabled!","iftWriteImagePNG");
    #endif
}

void iftWriteImageJPEG(const iftImage* img, const char* format, ...)
{
    #if IFT_LIBJPEG
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    //code based on external/libjpeg/source/example.c
    /* This struct contains the JPEG compression parameters and pointers to
 * working space (which is allocated as needed by the JPEG library).
 * It is possible to have several such structures, representing multiple
 * compression/decompression processes, in existence at once.  We refer
 * to any one struct (and its associated working data) as a "JPEG object".
 */
    struct jpeg_compress_struct cinfo;

    /* This struct represents a JPEG error handler.  It is declared separately
 * because applications often want to supply a specialized error handler
 * (see the second half of this file for an example).  But here we just
 * take the easy way out and use the standard error handler, which will
 * print a message on stderr and call exit() if compression fails.
 * Note that this struct must live as long as the main JPEG parameter
 * struct, to avoid dangling-pointer problems.
 */
    struct jpeg_error_mgr jerr;

    /* More stuff */
    FILE * outfile;		/* target file */
    JSAMPARRAY buffer;
    /* Step 1: allocate and initialize JPEG compression object */

    /* We have to set up the error handler first, in case the initialization
     * step fails.  (Unlikely, but it could happen if you are out of memory.)
     * This routine fills in the contents of struct jerr, and returns jerr's
     * address which we place into the link field in cinfo.
     */
    cinfo.err = jpeg_std_error(&jerr);

    /* Now we can initialize the JPEG compression object. */
    jpeg_create_compress(&cinfo);
    /* Step 2: specify data destination (eg, a file) */
    /* Note: steps 2 and 3 can be done in either order. */

    /* Here we use the library-supplied code to send compressed data to a
     * stdio stream.  You can also write your own code to do something else.
     * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
     * requires it in order to write binary files.
     */
    if ((outfile = fopen(filename, "wb")) == NULL) {
        fprintf(stderr, "can't open %s\n", filename);
        exit(1);
    }
    jpeg_stdio_dest(&cinfo, outfile);

    /* First we supply a description of the input image.
* Four fields of the cinfo struct must be filled in:
*/


    cinfo.image_width = img->xsize; 	/* image width and height, in pixels */
    cinfo.image_height = img->ysize;
    cinfo.data_precision = iftImageDepth(img);

    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_YCbCr;
    cinfo.jpeg_color_space = JCS_YCbCr;



    /* Now use the library's routine to set default compression parameters.
* (You must set at least cinfo.in_color_space before calling this,
* since the defaults depend on the source color space.)
*/
    jpeg_set_defaults(&cinfo);

    /* Now you can set any non-default parameters you wish to.
* Here we just illustrate the use of quality (quantization table) scaling:
*/
    int quality = 100;
    jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

    /* Step 4: Start compressor */
    /* TRUE ensures that we will write a complete interchange-JPEG file.
     * Pass TRUE unless you are very sure of what you're doing.
     */
    jpeg_start_compress(&cinfo, TRUE);

    /* Step 5: while (scan lines remain to be written) */
    /*           jpeg_write_scanlines(...); */

    /* Here we use the library's state variable cinfo.next_scanline as the
     * loop counter, so that we don't have to keep track ourselves.
     * To keep things simple, we pass one scanline per call; you can pass
     * more if you wish, though.
     */
    int row_stride = cinfo.image_width * cinfo.num_components;
    buffer = (*cinfo.mem->alloc_sarray)
            ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
    unsigned int imageRow = 0;
    while (cinfo.next_scanline < cinfo.image_height) {
        /* jpeg_write_scanlines expects an array of pointers to scanlines.
         * Here the array is only one element long, but you could pass
         * more than one scanline at a time if that's more convenient.
         */
        unsigned int imageCol = 0;
        for (unsigned int i = 0; i < (unsigned int)cinfo.image_width; i++) {
            int shift = i*cinfo.num_components;
            buffer[0][(shift+0)] = (unsigned char)iftImgVal(img,imageCol,imageRow,0);
            buffer[0][(shift+1)] = (unsigned char)iftImgCb(img,imageCol,imageRow,0);
            buffer[0][(shift+2)] = (unsigned char)iftImgCr(img,imageCol,imageRow,0);

            imageCol++;
        }
        imageRow++;
        (void) jpeg_write_scanlines(&cinfo, buffer, 1);
    }

    /* Step 6: Finish compression */

    jpeg_finish_compress(&cinfo);
    /* After finish_compress, we can close the output file. */
    fclose(outfile);

    /* Step 7: release JPEG compression object */

    /* This is an important step since it will release a good deal of memory. */
    jpeg_destroy_compress(&cinfo);
    #else
    iftError("LibPNG support was not enabled!","iftWriteImageJPEG");
    #endif
}

int iftMaximumValueInRegion(const iftImage *img, iftBoundingBox bb) 
{
    // checkers
    if (img == NULL)
        iftError("Image is NULL", "iftMaximumValueInRegion");
    if (!iftValidVoxel(img, bb.begin))
        iftError("Beginning voxel (%d, %d, %d) from Region (Bound. Box) is not in the Image Domain\n" \
                 "Img (xsize, ysize, zsize): (%d, %d, %d)", "iftMaximumValueInRegion",
                 bb.begin.x, bb.begin.y, bb.begin.z, img->xsize, img->ysize, img->zsize);
    if (!iftValidVoxel(img, bb.end))
        iftError("Ending voxel (%d, %d, %d) from Region (Bound. Box) is not in the Image Domain\n" \
                 "Img (xsize, ysize, zsize): (%d, %d, %d)", "iftMaximumValueInRegion",
                 bb.end.x, bb.end.y, bb.end.z, img->xsize, img->ysize, img->zsize);

    int img_max_val = IFT_INFINITY_INT_NEG;

    iftVoxel v;
    for (v.z = bb.begin.z; v.z <= bb.end.z; v.z++) {
        for (v.y = bb.begin.y; v.y <= bb.end.y; v.y++) {
            for (v.x = bb.begin.x; v.x <= bb.end.x; v.x++) {
                int p = iftGetVoxelIndex(img, v);
                if (img_max_val < img->val[p]) {
                    img_max_val = img->val[p];
                }
            }
        }
    }

    return img_max_val;
}

void iftSetImage(iftImage *img, int value) 
{
    for (int p = 0; p < img->n; p++)
        img->val[p] = value;
}

void  iftSetAlpha(iftImage *img, ushort value)
{
    int p;
    if(img->alpha == NULL){
        img->alpha = iftAllocUShortArray(img->n);
    }

    for (p=0; p < img->n; p++) {
        img->alpha[p] = value;
    }
}

void    iftSetCbCr(iftImage *img, ushort value)
{
    int p;

    if (!iftIsColorImage(img)){
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
    }
    for (p=0; p < img->n; p++) {
        img->Cb[p] = value;
        img->Cr[p] = value;
    }
}

void iftVerifyImageDomains(const iftImage *img1, const iftImage *img2, const char *function) 
{
    if ((img1==NULL)||(img2==NULL)||(img1->xsize!=img2->xsize) || (img1->ysize!=img2->ysize) || (img1->zsize!=img2->zsize)) {
        iftError("Images with different domains:\n" \
                "img1 (xsize, ysize, zsize): (%d, %d, %d)\n" \
                "img2 (xsize, ysize, zsize): (%d, %d, %d)\n",
                 function,
                 img1->xsize, img1->ysize, img1->zsize, img2->xsize, img2->ysize, img2->zsize);
    }
}

uchar iftImageDepth(const iftImage *img) 
{
    int img_min, img_max;
    iftMinMaxValues(img, &img_min, &img_max);
    
    long max_range;

    if (img_min >= 0)
        max_range = iftNormalizationValue(img_max) + 1;
    else
        max_range = iftNormalizationValue(img_max - img_min) + 1;
    
    return (uchar) iftLog(max_range, 2);
}

void iftMinMaxValues(const iftImage *img, int *min, int *max) 
{
    *min = *max = img->val[0];

    for (int p = 1; p < img->n; p++) {
        if (img->val[p] < *min)
            *min = img->val[p];
        else if (img->val[p] > *max)
            *max = img->val[p];
    }
}

iftImage *iftCreateImageFromImage(const iftImage *src) 
{
    iftImage *out = NULL;

    if (src != NULL) {
        if (iftIsColorImage(src)) {
            out = iftCreateColorImage(src->xsize, src->ysize, src->zsize, iftImageDepth(src));
        }
        else {
            out = iftCreateImage(src->xsize, src->ysize, src->zsize);
        }
        iftCopyVoxelSize(src, out);
    }

    return out;
}

iftImage  *iftCreateColorImage(int xsize,int ysize,int zsize, int depth)
{
    iftImage *img=NULL;
    img = iftCreateImage(xsize, ysize, zsize);

    iftSetCbCr(img, (iftMaxImageRange(depth)+1)/2);

    return(img);
}

void iftPutXYSlice(iftImage *img, const iftImage *slice, int zcoord)
{
    iftVoxel  u;
    int       p,q;

    if ( (zcoord < 0) || (zcoord >= img->zsize))
        iftError("Invalid z coordinate", "iftPutXYSlice");

    if ( (img->ysize!=slice->ysize)||(img->xsize!=slice->xsize) )
        iftError("Image and slice are incompatibles", "iftPutXYSlice");

    u.z   = zcoord;
    p     = 0;
    #if IFT_OMP
    #pragma omp parallel for private(p,q)
    #endif
    for (int y = 0; y < img->ysize; y++)
        for (int x = 0; x < img->xsize; x++)
        {
            u.x = x; u.y = y;
            q = iftGetVoxelIndex(img,u);
            img->val[q] = slice->val[p];
            if(iftIsColorImage(img))
            {
                img->Cb[q] = slice->Cb[p];
                img->Cr[q] = slice->Cr[p];
            }
            p++;
        }
}

iftImage *iftGetXYSlice(const iftImage *img, int zcoord)
{
    iftImage *slice;
    iftVoxel  u;
    int       p,q;

    if ( (zcoord < 0) || (zcoord >= img->zsize))
        iftError("Invalid z coordinate", "iftGetXYSlice");

    if(iftIsColorImage(img))
        slice = iftCreateColorImage(img->xsize,img->ysize,1, iftImageDepth(img));
    else
        slice = iftCreateImage(img->xsize,img->ysize,1);

    u.z   = zcoord;
    q     = 0;
    #if IFT_OMP
    #pragma omp parallel for private(p,q)
    #endif
    for (int y = 0; y < img->ysize; y++)
        for (int x = 0; x < img->xsize; x++)
        {
            u.x = x; u.y = y;
            p = iftGetVoxelIndex(img,u);
            slice->val[q] = img->val[p];
            if(iftIsColorImage(img))
            {
                slice->Cb[q] = img->Cb[p];
                slice->Cr[q] = img->Cr[p];
            }
            q++;
        }
    iftCopyVoxelSize(img,slice);

    return(slice);
}

iftImage *iftBMapToBinImage(const iftBMap *bmap, int xsize, int ysize, int zsize) {
    if (bmap == NULL)
        iftError("Bin Map is NULL", "iftBMapToBinImage");

    int n = xsize * ysize * zsize;

    if (bmap->n != n)
        iftError("Image Domain is != of the Bin Map Size\n" \
                 "Img Domain: (%d, %d, %d) = %d spels\nBin Map: %d elems", "iftBMapToBinImage",
                 xsize, ysize, zsize, n, bmap->n);

    iftImage *bin = iftCreateImage(xsize, ysize, zsize);

    for (int p = 0; p < bmap->n; p++)
        bin->val[p] = iftBMapValue(bmap, p);

    return bin;
}


iftBMap *iftBinImageToBMap(const iftImage *bin_img) {
    if (bin_img == NULL)
        iftError("Bin Image is NULL", "iftBinImageToBMap");

    iftBMap *bmap = iftCreateBMap(bin_img->n);

    for (int i = 0; i < bin_img->n; i++)
        if (bin_img->val[i])
            iftBMapSet1(bmap, i);

    return bmap;
}

// ---------- iftImage.c end 
// ---------- iftMatrix.c start 

iftMatrix *iftCreateMatrix(int ncols, int nrows) 
{
    iftMatrix *M = (iftMatrix *) iftAlloc(1, sizeof(iftMatrix));
    
    M->ncols = ncols;
    M->nrows = nrows;
    M->tbrow = (long*)iftAlloc(nrows,sizeof(long));
    //M->tbrow = iftAllocIntArray(nrows);
    for (long r = 0; r < (long)nrows; r++) {
        M->tbrow[r] = r * ncols;
    }
    M->n = (long) ncols * (long) nrows;
    M->allocated = true;
    
    M->val = iftAllocFloatArray(M->n);
    
    return (M);
}

iftMatrix *iftCopyMatrix(const iftMatrix *A) 
{
    iftMatrix *B;
    int i;

    B = iftCreateMatrix(A->ncols, A->nrows);

    for (i = 0; i < A->n; i++)
        B->val[i] = A->val[i];

    return (B);
}

void iftDestroyMatrix(iftMatrix **M) 
{
    if (M != NULL) {
        iftMatrix *aux = *M;
        
        if (aux != NULL) {
            if (aux->allocated && aux->val != NULL)
                iftFree(aux->val);
            iftFree(aux->tbrow);
            iftFree(aux);
        }
        *M = NULL;
    }
}

// ---------- iftMatrix.c end
// ---------- iftMImage.c start 

iftMImage * iftCreateMImage(int xsize,int ysize,int zsize, int nbands)
{
  iftMImage *img=NULL;
  int        i,y,z,xysize;

  img = (iftMImage *) iftAlloc(1,sizeof(iftMImage));
  if (img == NULL){
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateMImage");
  }

  img->n       = xsize*ysize*zsize;
  img->m       = nbands;
  img->data    = iftCreateMatrix(img->m, img->n);
  img->val     = iftAlloc(img->n, sizeof *img->val);
  for (i = 0; i < img->n; i++)
      img->val[i] = iftMatrixRowPointer(img->data, i);
  img->xsize   = xsize;
  img->ysize   = ysize;
  img->zsize   = zsize;
  img->dx      = 1.0;
  img->dy      = 1.0;
  img->dz      = 1.0;
  img->tby     = iftAllocIntArray(ysize);
  img->tbz     = iftAllocIntArray(zsize);

  img->tby[0]=0;
  for (y=1; y < ysize; y++)
    img->tby[y]=img->tby[y-1] + xsize;

  img->tbz[0]=0; xysize = xsize*ysize;
  for (z=1; z < zsize; z++)
    img->tbz[z]=img->tbz[z-1] + xysize;

  return(img);
}

void iftDestroyMImage(iftMImage **img)
{
    iftMImage *aux;

    aux = *img;
    if (aux != NULL) {
        if (aux->val != NULL)
            iftFree(aux->val);
        if (aux->data != NULL)
            iftDestroyMatrix(&aux->data);
        if (aux->tby != NULL)
            iftFree(aux->tby);
        if (aux->tbz != NULL)
            iftFree(aux->tbz);
        iftFree(aux);
        *img = NULL;
    }
}

iftMImage * iftImageToMImage(const iftImage *img1, char color_space)
{
  iftMImage *img2=NULL;
  int normalization_value = iftNormalizationValue(iftMaximumValue(img1));

  switch (color_space) {

  case YCbCr_CSPACE:
    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);
#if IFT_OMP
#pragma omp parallel for shared(img1, img2, normalization_value)
#endif
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)img1->val[p]);
      img2->val[p][1]=((float)img1->Cb[p]);
      img2->val[p][2]=((float)img1->Cr[p]);
    }
    break;

  case YCbCrNorm_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);

#if IFT_OMP
#pragma omp parallel for shared(img1, img2, normalization_value)
#endif
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)img1->val[p])/(float)normalization_value;
      img2->val[p][1]=((float)img1->Cb[p])/(float)normalization_value;
      img2->val[p][2]=((float)img1->Cr[p])/(float)normalization_value;
    }
    break;

  case LAB_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);

#if IFT_OMP
#pragma omp parallel for shared(img1, img2, normalization_value)
#endif
    for (int p=0; p < img2->n; p++) {
      iftColor  YCbCr,RGB;
      iftFColor Lab;
      YCbCr.val[0] = img1->val[p];
      YCbCr.val[1] = img1->Cb[p];
      YCbCr.val[2] = img1->Cr[p];
      RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
      Lab = iftRGBtoLab(RGB,normalization_value);
      img2->val[p][0]=Lab.val[0];
      img2->val[p][1]=Lab.val[1];
      img2->val[p][2]=Lab.val[2];
    }
    break;

    case LABNorm_CSPACE:

      img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);

#if IFT_OMP
#pragma omp parallel for shared(img1, img2, normalization_value)
#endif
      for (int p=0; p < img2->n; p++) {
        iftColor  YCbCr,RGB;
        iftFColor Lab;
        YCbCr.val[0] = img1->val[p];
        YCbCr.val[1] = img1->Cb[p];
        YCbCr.val[2] = img1->Cr[p];
        RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
        Lab = iftRGBtoLabNorm(RGB,normalization_value);
        img2->val[p][0]=Lab.val[0];
        img2->val[p][1]=Lab.val[1];
        img2->val[p][2]=Lab.val[2];
      }
      break;

  case LABNorm2_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);

#if IFT_OMP
#pragma omp parallel for shared(img1, img2, normalization_value)
#endif
    for (int p=0; p < img2->n; p++) {
      iftColor  YCbCr,RGB;
      iftFColor LabNorm;
      YCbCr.val[0] = img1->val[p];
      YCbCr.val[1] = img1->Cb[p];
      YCbCr.val[2] = img1->Cr[p];
      RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
      LabNorm = iftRGBtoLabNorm2(RGB,normalization_value);
      img2->val[p][0]=LabNorm.val[0];
      img2->val[p][1]=LabNorm.val[1];
      img2->val[p][2]=LabNorm.val[2];
    }
    break;

  case RGB_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);

#if IFT_OMP
#pragma omp parallel for shared(img1, img2, normalization_value)
#endif
    for (int p=0; p < img2->n; p++) {
      iftColor YCbCr,RGB;
      YCbCr.val[0] = img1->val[p];
      YCbCr.val[1] = img1->Cb[p];
      YCbCr.val[2] = img1->Cr[p];
      RGB          = iftYCbCrtoRGB(YCbCr,normalization_value);
      img2->val[p][0]=((float)RGB.val[0]);
      img2->val[p][1]=((float)RGB.val[1]);
      img2->val[p][2]=((float)RGB.val[2]);
    }
    break;

  case RGBNorm_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);

#if IFT_OMP
#pragma omp parallel for shared(img1, img2, normalization_value)
#endif
    for (int p=0; p < img2->n; p++) {
      iftColor YCbCr,RGB;
      YCbCr.val[0] = img1->val[p];
      YCbCr.val[1] = img1->Cb[p];
      YCbCr.val[2] = img1->Cr[p];
      RGB          = iftYCbCrtoRGB(YCbCr,normalization_value);
      img2->val[p][0]=((float)RGB.val[0])/(float)normalization_value;
      img2->val[p][1]=((float)RGB.val[1])/(float)normalization_value;
      img2->val[p][2]=((float)RGB.val[2])/(float)normalization_value;
    }
    break;

  case GRAY_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,1);

#if IFT_OMP
#pragma omp parallel for shared(img1, img2, normalization_value)
#endif
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)img1->val[p]);
    }
    break;

  case GRAYNorm_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,1);

#if IFT_OMP
#pragma omp parallel for shared(img1, img2, normalization_value)
#endif
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)img1->val[p])/(float)normalization_value;
    }

    break;

  case WEIGHTED_YCbCr_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);


#if IFT_OMP
#pragma omp parallel for shared(img1, img2, normalization_value)
#endif
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=(0.2/2.2)*((float)img1->val[p]/(float)normalization_value);
      img2->val[p][1]=(1.0/2.2)*((float)img1->Cb[p]/(float)normalization_value);
      img2->val[p][2]=(1.0/2.2)*((float)img1->Cr[p]/(float)normalization_value);
    }
    break;
  
  case HSV_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);

#if IFT_OMP
#pragma omp parallel for shared(img1, img2, normalization_value)
#endif
    for (int p=0; p < img2->n; p++) {
      iftColor  YCbCr,RGB;
      iftColor HSV;
      YCbCr.val[0] = img1->val[p];
      YCbCr.val[1] = img1->Cb[p];
      YCbCr.val[2] = img1->Cr[p];
      RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
      HSV = iftRGBtoHSV(RGB,normalization_value);
      img2->val[p][0]=HSV.val[0];
      img2->val[p][1]=HSV.val[1];
      img2->val[p][2]=HSV.val[2];
    }
    break;

  default:
      iftError("Invalid color space (see options in iftColor.h)", "iftImageToMImage");
  }

  img2->dx = img1->dx;
  img2->dy = img1->dy;
  img2->dz = img1->dz;

  return(img2);

}

iftImage * iftMImageToImage(const iftMImage *img1, int Imax, int band)
{
  iftImage *img2=iftCreateImage(img1->xsize,img1->ysize,img1->zsize);
  int p,b=band;
  double min = IFT_INFINITY_FLT, max = IFT_INFINITY_FLT_NEG;

  if ((band < 0)||(band >= img1->m))
      iftError("Invalid band", "iftMImageToImage");

  for (p=0; p < img1->n; p++) {
    if (img1->val[p][b] < min)
      min = img1->val[p][b];
    if (img1->val[p][b] > max)
      max = img1->val[p][b];
  }

  //printf("min %lf max %lf\n",min,max);

  if (max > min){
    for (p=0; p < img2->n; p++) {
      img2->val[p]=(int)(Imax*(img1->val[p][b]-min)/(max-min));
    }
  }else{
    char msg[100];
    sprintf(msg,"Image is empty: max = %f and min = %f\n",max,min);
    // iftWarning(msg,"iftMImageToImage");
  }

  img2->dx = img1->dx;
  img2->dy = img1->dy;
  img2->dz = img1->dz;

  return(img2);
}

inline iftVoxel iftMGetVoxelCoord(const iftMImage *img, int p)
{
    /* old
     * u.x = (((p) % (((img)->xsize)*((img)->ysize))) % (img)->xsize)
     * u.y = (((p) % (((img)->xsize)*((img)->ysize))) / (img)->xsize)
     * u.z = ((p) / (((img)->xsize)*((img)->ysize)))
     */
    iftVoxel u;
    div_t res1 = div(p, img->xsize * img->ysize);
    div_t res2 = div(res1.rem, img->xsize);

    u.x = res2.rem;
    u.y = res2.quot;
    u.z = res1.quot;

  return(u);
}

inline char iftMValidVoxel(const iftMImage *img, iftVoxel v) 
{
  if ((v.x >= 0)&&(v.x < img->xsize)&&
      (v.y >= 0)&&(v.y < img->ysize)&&
      (v.z >= 0)&&(v.z < img->zsize))
    return(1);
  else
    return(0);
}

float iftMMaximumValue(const iftMImage *img, int band) {
  int b, i;
  float max_val = IFT_INFINITY_FLT_NEG;

  if(band < 0) {
    for(b = 0; b < img->m; b++) {
      for(i = 0; i < img->n; i++) {
  max_val = iftMax(img->val[i][b], max_val);
      }
    }
  } else {
    for(i = 0; i < img->n; i++) {
      max_val = iftMax(img->val[i][band], max_val);
    }
  }

  return max_val;
}

// ---------- iftMImage.c end
// ---------- iftKernel.c start
iftKernel *iftCreateKernel(iftAdjRel *A)
{
  iftKernel *K=(iftKernel *) iftAlloc(1,sizeof(iftKernel));

  K->A      = iftCopyAdjacency(A); 
  K->weight = iftAllocFloatArray(K->A->n);

  return(K);
}

void iftDestroyKernel(iftKernel **K)
{
  iftKernel *aux=*K;

  if (aux != NULL) {
    iftDestroyAdjRel(&aux->A);
    iftFree(aux->weight);
    iftFree(aux);
    *K = NULL;
  }
}
// ---------- iftKernel.c end
// ---------- iftMemory.c start 

#ifdef __linux__
#include <sys/sysinfo.h>
#include <malloc.h>
#endif

#ifdef __APPLE__
#include <mach/task.h>
#include <mach/mach_init.h>
#endif

#ifdef _WINDOWS
#include <windows.h>
#else
#include <sys/resource.h>
#endif

int *iftAllocIntArray(long n) 
{
    int *v = NULL;

    v = (int *) iftAlloc(n, sizeof(int));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocIntArray");
    return(v);
}

void iftCopyIntArray(int *array_dst, const int *array_src, int nelems) 
{
    #if IFT_OMP
    #pragma omp parallel for
    #endif
    
    for (int i = 0; i < nelems; i++) {
        array_dst[i] = array_src[i];
    }
}

#ifndef  __cplusplus
long long *iftAllocLongLongIntArray(long n) 
{
    long long *v = NULL;

    v = (long long *) iftAlloc(n, sizeof(long long));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocLongLongIntArray");

    return v;
}

void iftCopyLongLongIntArray(long long *array_dst, const long long *array_src, int nelems) 
{
    #if IFT_OMP
    #pragma omp parallel for
    #endif
    
    for (int i = 0; i < nelems; ++i) {
        array_dst[i] = array_src[i];
    }
}
#endif

float *iftAllocFloatArray(long n) 
{
    float *v = NULL;
    v = (float *) iftAlloc(n, sizeof(float));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocFloatArray");
    return(v);
}

void iftCopyFloatArray(float *array_dst, float *array_src, int nelems) 
{

    memmove(array_dst, array_src, nelems*sizeof(float));
}

double *iftAllocDoubleArray(long n) 
{
    double *v = NULL;

    v = (double *) iftAlloc(n, sizeof(double));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocDoubleArray");
    return (v);
}

void iftCopyDoubleArray(double *array_dst, double *array_src, int nelems) 
{
    #if IFT_OMP
    #pragma omp parallel for
    #endif
    
    for (int i = 0; i < nelems; i++)
        array_dst[i] = array_src[i];
}

ushort *iftAllocUShortArray(long n) 
{
    ushort *v = NULL;

    v = (ushort *) iftAlloc(n, sizeof(ushort));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocUShortArray");
    return(v);
}

uchar *iftAllocUCharArray(long n) 
{
    uchar *v = NULL;

    v = (uchar *) iftAlloc(n, sizeof(uchar));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocUCharArray");
    return (v);
}

char *iftAllocCharArray(long n) 
{
    char *v = NULL;

    v = (char *) iftAlloc(n, sizeof(char));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocCharArray");
    return (v);
}

char *iftAllocString(long n) 
{
    return (iftAllocCharArray(n+1));
}

// ---------- iftMemory.c end
// ---------- iftSet.c start 

void iftInsertSet(iftSet **S, int elem)
{
    iftSet *p=NULL;
    
    p = (iftSet *) iftAlloc(1,sizeof(iftSet));
    if (p == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftInsertSet");
    if (*S == NULL){
        p->elem  = elem;
        p->next  = NULL;
    }else{
        p->elem  = elem;
        p->next  = *S;
    }
    *S = p;
}

int iftRemoveSet(iftSet **S)
{
    iftSet *p;
    int elem = IFT_NIL;
    
    if (*S != NULL){
        p    =  *S;
        elem = p->elem;
        *S   = p->next;
        iftFree(p);
    }
    
    return(elem);
}

void iftRemoveSetElem(iftSet **S, int elem)
{
    if (S == NULL || *S == NULL)
        return;
    
    iftSet *tmp = *S;
    
    if (tmp->elem == elem) {
        *S = tmp->next;
        iftFree(tmp);
    } else {
        while (tmp->next != NULL && tmp->next->elem != elem)
            tmp = tmp->next;
        if (tmp->next == NULL)
            return;
        iftSet *remove = tmp->next;
        tmp->next = remove->next;
        iftFree(remove);
    }
}

void iftDestroySet(iftSet **S)
{
    iftSet *p;
    while(*S != NULL){
        p = *S;
        *S = p->next;
        iftFree(p);
    }
    *S = NULL;
}

iftSet* iftSetUnion(iftSet *S1,iftSet *S2)
{
    iftSet *S = 0;
    
    iftSet *s = S1;
    while(s){
        iftInsertSet(&S,s->elem);
        s = s->next;
    }
    
    s = S2;
    while(s){
        iftUnionSetElem(&S,s->elem);
        s = s->next;
    }
    
    return S;
}

iftSet* iftSetConcat(iftSet *S1,iftSet *S2)
{
    iftSet *S = 0;
    
    iftSet *s = S1;
    while(s){
        iftInsertSet(&S,s->elem);
        s = s->next;
    }
    
    s = S2;
    while(s){
        iftInsertSet(&S,s->elem);
        s = s->next;
    }
    
    return S;
}

char iftUnionSetElem(iftSet **S, int elem)
{
    iftSet *aux=*S;
    
    while (aux != NULL) {
        if (aux->elem == elem)
            return(0);
        aux = aux->next;
    }
    iftInsertSet(S,elem);
    return(1);
}

void iftInvertSet(iftSet **S)
{
    iftSet *set=NULL;
    
    while (*S != NULL)
        iftInsertSet(&set,iftRemoveSet(S));
    *S = set;
}

int iftSetSize(const iftSet* S)
{
    const iftSet *s = S;
    
    int i = 0;
    while (s != NULL){
        i++;
        s = s->next;
    }
    
    return i;
    
}

iftSet* iftSetCopy(iftSet* S)
{
    return iftSetUnion(S,0);
}

int iftSetHasElement(iftSet *S, int elem)
{
    iftSet *s = S;
    while(s){
        if(s->elem == elem)
            return 1;
        
        s = s->next;
    }
    
    return 0;
}

iftIntArray *iftSetToArray(iftSet *S) {
    int n_elems = iftSetSize(S);
    iftIntArray *array = iftCreateIntArray(n_elems);
    
    iftSet *sp = S;
    int i      = 0;
    while (sp != NULL) {
        array->val[i++] = sp->elem;
        sp = sp->next;
    }
    
    return array;
}

// ---------- iftSet.c end
// ---------- iftSList.c start

iftSList *iftCreateSList() 
{
    iftSList *SL = (iftSList *) iftAlloc(1, sizeof(iftSList));

    // just to really force
    SL->n    = 0;
    SL->head = NULL;
    SL->tail = NULL;

    return SL;
}

void iftDestroySList(iftSList **SL) 
{
    if (SL != NULL) {
        iftSList *SL_aux = *SL;

        if (SL_aux != NULL) {
            iftSNode *snode = SL_aux->head;
            iftSNode *tnode = NULL;

            while (snode != NULL) {
                tnode = snode;
                snode = snode->next;

                if (tnode->elem != NULL)
                    iftFree(tnode->elem);
                iftFree(tnode);
            }
            iftFree(SL_aux);
            *SL = NULL;
        }
    }
}

void iftInsertSListIntoHead(iftSList *SL, const char *elem) 
{
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftInsertSListIntoHead");

    iftSNode *snode = (iftSNode*) iftAlloc(1, sizeof(iftSNode));
    snode->elem     = iftAllocCharArray(512);
    snode->prev = NULL;
    snode->next     = NULL;
    strcpy(snode->elem, elem);


    // The String Linked List is empty
    if (SL->head == NULL) {
        SL->head = snode;
        SL->tail = snode;
        SL->n    = 1;
    }
    else {
        snode->next        = SL->head;
        SL->head->prev = snode;
        SL->head           = snode;
        SL->n++;
    }
}

void iftInsertSListIntoTail(iftSList *SL, const char *elem) 
{
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftInsertSListIntoTail");
    if (elem == NULL)
        iftError("The Element to be Inserted is NULL", "iftInsertSListIntoTail");

    iftSNode *snode = (iftSNode*) iftAlloc(1, sizeof(iftSNode));
    snode->prev     = NULL;
    snode->next     = NULL;
    snode->elem     = iftAllocCharArray(strlen(elem) + 1);
    strcpy(snode->elem, elem);


    // The String Linked List is empty
    if (SL->head == NULL) {
        SL->head = snode;
        SL->tail = snode;
        SL->n    = 1;
    }
    else {
        SL->tail->next = snode;
        snode->prev    = SL->tail;
        SL->tail       = snode;
        SL->n++;
    }
}

char *iftRemoveSListHead(iftSList *SL) 
{
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftRemoveSListHead");
    
    char *elem      = NULL;
    iftSNode *snode = NULL;

    // if there are elements
    if (SL->head != NULL) {
        snode    = SL->head;
        SL->head = SL->head->next;

        // checks if the list is empty now
        if (SL->head == NULL)
            SL->tail = NULL;
        else
            SL->head->prev = NULL;

        SL->n--;
        elem = snode->elem;
        snode->elem = NULL;
        iftFree(snode); // it does not deallocates snode->free
    }

    return elem;
}

char *iftRemoveSListTail(iftSList *SL) 
{
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftRemoveSListTail");

    char *elem      = NULL;
    iftSNode *snode = NULL;

    // if there are elements
    if (SL->head != NULL) {
        snode    = SL->tail;
        SL->tail = SL->tail->prev;

        // checks if the list is empty now
        if (SL->tail == NULL)
            SL->head = NULL;
        else
            SL->tail->next = NULL;

        SL->n--;
        elem = snode->elem;
        snode->elem = NULL;
        iftFree(snode); // it does not deallocates snode->free
    }

    return elem;
}

// ---------- iftSList.c end
// ---------- iftDialog.c start 

void iftError(const char *msg, const char *func, ...) 
{
    va_list args;
    char    final_msg[4096];
    
    va_start(args, func);
    vsprintf(final_msg, msg, args);
    va_end(args);
    
    fprintf(stderr, "\nError in %s: \n%s\n", func, final_msg);
    fflush(stdout);
    exit(-1);
}

void iftWarning(const char *msg, const char *func, ...) 
{
    va_list args;
    char    final_msg[4096];
    
    va_start(args, func);
    vsprintf(final_msg, msg, args);
    va_end(args);
    
    fprintf(stdout, "\nWarning in %s: \n%s\n", func, final_msg);
}

// ---------- iftDialog.c end
// ---------- iftStream.c start 
char *iftGetLine(FILE *stream) {
    if (stream == NULL || feof(stream))
        return NULL;
    
    size_t buffer_size = 0;
    char *line = NULL;
    
    if (getline(&line, &buffer_size, stream) == -1) {
        if (feof(stream)) {
            if (line != NULL)
                free(line);
            return NULL;
        } else
            iftError("Error with getline command", "iftGetLine");
    }
    
    size_t str_len = strlen(line);
    if (line[str_len-1] == '\n')
        line[str_len-1] = '\0'; // without this, it reads the \n and does not put the \0 at the end
    
    return line;
}
// ---------- iftStream.c end
// ---------- iftString.c start 

void iftRightTrim(char* s, char c) 
{

    int idx = strlen(s) - 1;

    while(s[idx] == c) {
        idx--;
    }

    s[idx+1] = '\0';

}

iftSList *iftSplitString(const char *phrase, const char *delimiter) 
{
    if (phrase == NULL)
        iftError("String to be splitted is NULL", "iftSplitString");
    if (delimiter == NULL)
        iftError("Delimiter is NULL", "iftSplitString");

    char *buf       = iftAllocString(strlen(phrase)+1);
    const char *pr  = phrase;
    const char *loc = strstr(pr, delimiter); // pointer to first delimiter occurrence in the string

    size_t length = strlen(delimiter);
    size_t bytes;

    iftSList *SL = iftCreateSList();

    // build a list of sub-strings
    while (loc != NULL) {
        bytes = loc - pr;
        strncpy(buf, pr, bytes);
        buf[bytes] = '\0'; // \0 character must be added manually because strncpy does not do that, as opposed to other functions such as strcpy and sprintf

        iftInsertSListIntoTail(SL, buf);

        pr = loc + length;
        loc = strstr(pr, delimiter);
    }

    // Copies the last substring to the left of the last delimiter found OR
    // Copies the whole string if it doesn't have the delimiter 
    strcpy(buf, pr);
    iftInsertSListIntoTail(SL, buf);

    iftFree(buf);

    return SL;
}

char *iftLowerString(const char *str) 
{
    if (str == NULL)
        iftError("Input string is NULL", "iftLowerString");

    char *out_str = iftAllocCharArray(strlen(str)+1);

    for (size_t c = 0; c < strlen(str); c++)
        out_str[c] = tolower(str[c]);

    return out_str;
}

bool iftCompareStrings(const char *str1, const char *str2) 
{
    if (str1 == NULL)
        iftError("First String is NULL", "iftCompareStrings");
    if (str2 == NULL)
        iftError("Second String is NULL", "iftCompareStrings");

    return (strcmp(str1, str2) == 0);
}

char *iftSplitStringAt(const char *phrase, const char *delimiter, long position) 
{
    if (phrase == NULL)
        iftError("String to be splitted is NULL", "iftSplitStringAt");
    if (delimiter == NULL)
        iftError("Delimiter is NULL", "iftSplitStringAt");

    iftSList *SL = iftSplitString(phrase, delimiter);
    iftSNode *snode = NULL;

    // Copies the split sub-string of the position
    if (position >= 0) {
        if ((position+1) <= SL->n) {

            snode = SL->head;
            for (size_t i = 0; i < (ulong)position; i++)
                snode = snode->next;
        }
        else {
            iftError("Invalid Position %ld\n-> Position Index must be < %ld\n",
                     "iftSplitStringAt", position, SL->n);
        }
    } else {
        if (labs(position) <= SL->n) {
            long real_pos = SL->n + position;
            snode = SL->tail;
            for (size_t i = SL->n-1; i > (ulong)real_pos; i--) {
                snode = snode->prev;
            }
        }
        else {
            iftError("Invalid Negative Position %ld\n-> Negative Position Index must be >= %ld\n",
                     "iftSplitStringAt", position, -1 * SL->n);
        }
    }

    char *str = iftCopyString(snode->elem);

    iftDestroySList(&SL);

    return str;
}

char *iftCopyString(const char *format, ...) 
{
    va_list args;
    char str[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(str, format, args);
    va_end(args);
    
    char *copy = iftAllocCharArray(strlen(str) + 1);
    strcpy(copy, str);

    return copy;
}

char *iftRemoveSuffix(const char *str, const char *suffix) 
{
    if (str == NULL)
        iftError("String is NULL", "iftRemoveSuffix");
    if (suffix == NULL)
        iftError("Suffix is NULL", "iftRemoveSuffix");

    if (!iftCompareStrings(suffix, "") && iftEndsWith(str, suffix)) {
        size_t shift = strlen(str) - strlen(suffix);
        char *out_str = iftCopyString(str);
        out_str[shift] = '\0';

        return (out_str);   
    }
    else {
        return iftCopyString(str);
    }
}

bool iftEndsWith(const char *str, const char *suffix) 
{
    if (str == NULL)
        iftError("String is NULL", "iftEndsWith");
    if (suffix == NULL)
        iftError("Suffix is NULL", "iftEndsWith");

    size_t len_suffix = strlen(suffix);
    size_t len_str    = strlen(str);

    if (len_suffix <= len_str) {
        size_t shift = len_str - len_suffix;
        return (strncmp(str+shift, suffix, len_suffix) == 0);
    }
    else
        return false;
}

bool iftStartsWith(const char *str, const char *prefix) 
{
    if (str == NULL)
        iftError("String is NULL", "iftStartsWith");
    if (prefix == NULL)
        iftError("Prefix is NULL", "iftStartsWith");

    size_t len_prefix = strlen(prefix);
    size_t len_str    = strlen(str);

    if (len_prefix <= len_str)
        return (strncmp(str, prefix, len_prefix) == 0);        
    else
        return false;
}

char *iftConcatStrings(int n, ...) 
{
    if (n <= 0)
        iftError("Number of Strings to be concatenated is <= 0", "iftConcatStrings");

    size_t out_str_size = 1; // '\0'

    // Counts the size of the concatenated string
    va_list strings;
    va_start(strings, n);
    for (int i = 0; i < n; i++)
        out_str_size += strlen(va_arg(strings, char*));
    va_end(strings);

    char *concat_str = iftAllocCharArray(out_str_size);

    va_start(strings, n);
    for (int i = 0; i < n; i++)
        strcat(concat_str, va_arg(strings, char*));
    va_end(strings);

    return concat_str;
}

char *iftRemovePrefix(const char *str, const char *prefix) 
{
    if (str == NULL)
        iftError("String is NULL", "iftRemovePrefix");
    if (prefix == NULL)
        iftError("Prefix is NULL", "iftRemovePrefix");

    if (!iftCompareStrings(prefix, "") && iftStartsWith(str, prefix)) {
        size_t shift = strlen(prefix);
        return (iftCopyString(str + shift));   
    }
    else 
        return (iftCopyString(str));
}

char *iftReplaceString(const char *str, const char *old_sub, const char *new_sub) {
    if (str == NULL)
        iftError("String is NULL", "iftReplaceString");
    if (old_sub == NULL)
        iftError("Old Substring is NULL", "iftReplaceString");
    if (new_sub == NULL)
        iftError("New Substring is NULL", "iftReplaceString");

    iftSList *SL = iftSplitString(str, old_sub);

    long n_sub   = SL->n-1; // number of old subtrings found in str
    size_t str_len = strlen(str) + (n_sub * strlen(new_sub)) + 1; // size of the replaced string (with '\0')

    // adds the number of chars of the new substring + '\0'
    str_len += (n_sub * strlen(new_sub)) + 1;
    char *rep_str = iftAllocCharArray(str_len);

    // builds the replaced String - the list has always at least one element
    char *elem = iftRemoveSListHead(SL);
    while ((elem != NULL) && (SL->n >= 1)) {
        strcat(rep_str, elem);
        strcat(rep_str, new_sub);
        elem = iftRemoveSListHead(SL);
    }
    strcat(rep_str, elem); // copies the last element from the List

    iftDestroySList(&SL);

    return rep_str;
}

// ---------- iftString.c end
// ---------- iftDir.c start 

int _iftCmpFiles(const void *a, const void *b) 
{
    iftFile **f1 = (iftFile**) a;
    iftFile **f2 = (iftFile**) b;
    
    return strcmp((*f1)->path, (*f2)->path);
}

int _iftCmpDirs(const void *a, const void *b) 
{
    iftDir **dir1 = (iftDir**) a;
    iftDir **dir2 = (iftDir**) b;
    
    return strcmp((*dir1)->path, (*dir2)->path);
}

void _iftCountFilesInDirectory(const char *dir_pathname, long *nfiles, long *nsubdirs) 
{
    //http://pubs.opengroup.org/onlinepubs/007908799/xsh/dirent.h.html
    //http://www.delorie.com/gnu/docs/glibc/libc_270.html
    DIR *system_dir;
    struct dirent *entry;
    char msg[512];
    char *pathname = NULL;
    *nfiles = 0;
    *nsubdirs = 0;
    
    system_dir = opendir(dir_pathname);
    if (system_dir) {
        while ((entry = readdir(system_dir)) != NULL)
            // it excludes the system_dir . and ..
            if ((strcmp(entry->d_name, ".") != 0) && (strcmp(entry->d_name, "..") != 0)) {
                pathname = iftJoinPathnames(2, dir_pathname, entry->d_name);
                
                if (iftDirExists(pathname)) {
                    (*nsubdirs)++;
                }
                else {
                    (*nfiles)++;
                }
                
                iftFree(pathname);
                pathname = NULL;
            }
        closedir(system_dir);
    }
    else {
        sprintf(msg, "Error opening directory path: \"%s\"", dir_pathname);
        iftError(msg, "_iftCountFilesInDirectory");
    }
}

void _iftListDirectoryRec(iftDir *dir, long hier_levels, long curr_level) 
{
    DIR *system_dir;
    struct dirent *entry;
    char *pathname = NULL;
    
    dir->files   = NULL;
    dir->subdirs = NULL;
    
    if (curr_level <= hier_levels) {
        system_dir = opendir(dir->path);
        if (system_dir) {
            _iftCountFilesInDirectory(dir->path, &dir->nfiles, &dir->nsubdirs);
            if (dir->nfiles != 0)
                dir->files   = (iftFile**) iftAlloc(dir->nfiles, sizeof(iftFile*));
            if (dir->nsubdirs != 0)
                dir->subdirs = (iftDir**) iftAlloc(dir->nsubdirs, sizeof(iftDir*));
            
            long i = 0, j = 0;
            while ((entry = readdir(system_dir)) != NULL) {
                // it excludes the dir . and ..
                if ((strcmp(entry->d_name, ".") != 0) && (strcmp(entry->d_name, "..") != 0)) {
                    pathname = iftJoinPathnames(2, dir->path, entry->d_name);
                    
                    if (iftDirExists(pathname)) { // it is a directory
                        iftDir *subdir = (iftDir*) iftAlloc(1, sizeof(iftDir));
                        subdir->path = pathname;
                        
                        subdir->nfiles   = 0;
                        subdir->nsubdirs = 0;
                        subdir->files    = NULL;
                        subdir->subdirs  = NULL;
                        
                        _iftListDirectoryRec(subdir, hier_levels, curr_level+1);
                        dir->subdirs[j++] = subdir;
                    }
                    else { // it is a File
                        iftFile *f = (iftFile*) iftAlloc(1, sizeof(iftFile));
                        f->path = pathname;
                        dir->files[i++] = f;
                        f->suffix = NULL;
                    }
                }
            }
            closedir(system_dir);
            
            /* sorts the pathnames using qsort functions */
            qsort(dir->files, dir->nfiles, sizeof(iftFile*), _iftCmpFiles);
            qsort(dir->subdirs, dir->nsubdirs, sizeof(iftDir*), _iftCmpDirs);
        }
        else {
            char msg[512];
            sprintf(msg, "Error opening directory path: \"%s\"", dir->path);
            iftError(msg, "_iftListDirectoryRec");
        }
    }
}

void _iftListDirectory(iftDir *root_dir, long hier_levels) 
{
    if (root_dir == NULL) {
        iftError("Directory is NULL", "_iftListDirectory");
    }
    else {
        if (hier_levels == 0)
            hier_levels = IFT_INFINITY_INT; // trick to set the hier_levels as the possible maximum
        
        root_dir->nfiles   = 0;
        root_dir->nsubdirs = 0;
        root_dir->files    = NULL;
        root_dir->subdirs  = NULL;
        
        long curr_level = 1;
        _iftListDirectoryRec(root_dir, hier_levels, curr_level);
    }
}

bool iftDirExists(const char *format, ...) 
{
    va_list args;
    char pathname[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(pathname, format, args);
    va_end(args);
    
    struct stat st;
    
    if (stat(pathname, &st) == 0) {
        if (S_ISDIR(st.st_mode))
            return true; //it's a directory
    }
    return false ;
}

char *iftParentDir(const char *pathname) 
{
    char *filename   = NULL;
    char *parent_dir = NULL;
    
    filename   = iftSplitStringAt(pathname, IFT_SEP_C, -1);
    parent_dir = iftRemoveSuffix(pathname, filename);
    iftFree(filename);
    
    // if the parent_dir is empty
    if (strcmp(parent_dir, "") == 0) {
        strcpy(parent_dir, ".");
    }
    else { // eliminates the slash (dir. separator char) at the end
        parent_dir[strlen(parent_dir)-1] = '\0';
    }
    
    return (parent_dir);
}

void iftMakeDir(const char *dir_path) 
{
    if (!iftDirExists(dir_path)) {
        char *parent_dir = iftAllocCharArray(IFT_STR_DEFAULT_SIZE);
        strcpy(parent_dir, "");
        
        iftSList *SL    = iftSplitString(dir_path, IFT_SEP_C);
        char *inter_dir = iftRemoveSListHead(SL);
        
        while (inter_dir != NULL) {
            strcat(parent_dir, inter_dir);
            strcat(parent_dir, IFT_SEP_C);
            
            if (!iftDirExists(parent_dir)) {
            #if defined(__linux) || defined(__APPLE__)
                if (mkdir(parent_dir, 0777) == -1) // Create the directory
                    iftError("Problem to create the directory: %s", "iftMakeDir", dir_path);
            #else
                if (!CreateDirectory(parent_dir, NULL))
                        iftError("Problem to create the directory", "iftMakeDir");
            #endif
            }
            
            iftFree(inter_dir);
            inter_dir = iftRemoveSListHead(SL);
        }
        iftDestroySList(&SL);
        iftFree(parent_dir);
    }
}

iftDir *iftLoadFilesFromDirByRegex(const char *dir_pathname, const char *regex)
{
    if (dir_pathname == NULL)
        iftError("Dir's Pathname is NULL", "iftLoadFilesFromDirByRegex");
    if (regex == NULL)
        iftError("Regex is NULL", "iftLoadFilesFromDirByRegex");
    
    iftDir *dir          = iftLoadDir(dir_pathname, 1);
    long n_all_files   = dir->nfiles;
    long n_final_files = 0;
    
    for (long f = 0; f < n_all_files; f++) {
        char *filename = iftFilename(dir->files[f]->path, NULL);
        
        if (iftRegexMatch(filename, regex))
            n_final_files++;
        
        iftFree(filename);
    }
    
    iftFile **file_array = dir->files;
    dir->files           = NULL;
    
    dir->files  = (iftFile**) iftAlloc(n_final_files, sizeof(iftFile*));
    dir->nfiles = n_final_files;
    
    long i = 0;
    for (long f = 0; f < n_all_files; f++) {
        char *filename = iftFilename(file_array[f]->path, NULL);
        
        if (iftRegexMatch(filename, regex)) {
            dir->files[i++] = file_array[f];
            file_array[f]   = NULL;
        }
        else iftDestroyFile(&file_array[f]);
        
        iftFree(filename);
    }
    iftFree(file_array);
    
    return dir;
}

iftDir *iftLoadDir(const char *dir_pathname, long hier_levels) 
{
    char msg[512];
    iftDir *dir = NULL;
    
    
    if (iftPathnameExists(dir_pathname)) {
        // it is really a directory and it exists
        if (iftDirExists(dir_pathname)) {
            dir = (iftDir*) iftAlloc(1, sizeof(iftDir));
            dir->path = iftAllocCharArray(strlen(dir_pathname) + 2); // one more char to put the separation '/'
            strcpy(dir->path, dir_pathname);
            
            // puts the '/' at the end of the pathname
            if (dir->path[strlen(dir->path) - 1] != IFT_SEP_C[0])
                strcat(dir->path, IFT_SEP_C);
            
            _iftListDirectory(dir, hier_levels);
        }
            // it is a File instead of a Directory
        else {
            sprintf(msg, "Pathname \"%s\" is a File", dir_pathname);
            iftError(msg, "iftLoadDir");
        }
    }
    else {
        sprintf(msg, "Pathname \"%s\" does not exist!", dir_pathname);
        iftError(msg, "iftLoadDir");
    }
    
    return dir;
}

void iftDestroyDir(iftDir **dir) 
{
    if (dir != NULL) {
        iftDir *dir_aux = *dir;
        
        if (dir_aux != NULL) {
            if (dir_aux->path != NULL)
                iftFree(dir_aux->path);
            
            // deallocates the files
            if (dir_aux->files != NULL) {
                for (long i = 0; i < dir_aux->nfiles; i++)
                    iftDestroyFile(&dir_aux->files[i]);
                
                iftFree(dir_aux->files);
                dir_aux->files = NULL;
            }
            
            if (dir_aux->subdirs != NULL) {
                // deallocates the subdirs
                for (long j = 0; j < dir_aux->nsubdirs; j++)
                    iftDestroyDir(&dir_aux->subdirs[j]);
                
                iftFree(dir_aux->subdirs);
                dir_aux->subdirs = NULL;
            }
            iftFree(dir_aux);
            *dir = NULL;
        }
    }
}

// ---------- iftDir.c end
// ---------- iftRegex.c start

bool iftRegexMatch(const char *str, const char *regex_pattern, ...) {
    if (str == NULL)
        iftError("String is NULL", "iftRegexMatch");
    if (regex_pattern == NULL)
        iftError("Regular Expression is NULL", "iftRegexMatch");
    
    char error[IFT_STR_DEFAULT_SIZE];
    regex_t regex;
    int reti;
    
    va_list args;
    char final_regex_pattern[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, regex_pattern);
    vsprintf(final_regex_pattern, regex_pattern, args);
    va_end(args);
    
    // Compile Regex
    if ((reti = regcomp(&regex, final_regex_pattern, REG_EXTENDED|REG_NOSUB)) != 0) {
        regerror(reti, &regex, error, sizeof(error));
        iftError("Regex Compilation Failed: \"%s\"\n" \
                 "IFT_ERROR: %s", "iftRegexMatch", final_regex_pattern, error);
    }
    
    // Execute Regex
    reti = regexec(&regex, str, (size_t) 0, NULL, 0);
    regfree(&regex);
    
    return (reti == 0);
}

// ---------- iftRegex.c end 
// ---------- iftNumerical.c start

iftIntArray *iftIntRange(int begin, int end, int inc) 
{
    int n = ((end - begin) / inc) + 1;
    
    iftIntArray *space = iftCreateIntArray(n);
    
    for (int i = 0; i < n; ++i) {
        space->val[i] = begin + (inc*i);
    }
    
    return space;
}

// ---------- iftNumerical.c end
// ---------- iftSort.c start

void iftFQuickSort( float *value, int *index, int i0, int i1, uchar order ) 
{
  int m, d;


	if( i0 < i1 ) {
		/* random index to avoid bad pivots.*/
		// d = iftRandomInteger( i0, i1 );
    // to guarantee the same behavior on ties 
    d = (i0 + i1) / 2;
		iftSwap( value[ d ], value[ i0 ] );
		iftSwap( index[ d ], index[ i0 ] );
		m = i0;

		if(order == IFT_INCREASING ) {
			for( d = i0 + 1; d <= i1; d++ ) {
				if( value[ d ] < value[ i0 ] ) {
					m++;
					iftSwap( value[ d ], value[ m ] );
					iftSwap( index[ d ], index[ m ] );
				}
			}
		}
		else {
			for( d = i0 + 1; d <= i1; d++ ) {
				if( value[ d ] > value[ i0 ] ) {
					m++;
					iftSwap( value[ d ], value[ m ] );
					iftSwap( index[ d ], index[ m ] );
				}
			}
		}
		iftSwap( value[ m ], value[ i0 ] );
		iftSwap( index[ m ], index[ i0 ] );
		iftFQuickSort( value, index, i0, m - 1, order );
		iftFQuickSort( value, index, m + 1, i1, order );
	}
}

// ---------- iftSort.c end
// ---------- iftVideo.c start
int iftCountNumberOfFrameFolders(const char *path)
{
    int file_count = 0;

#ifdef __linux
    DIR * dirp;
    struct dirent * entry;

    dirp = opendir(path); /* There should be error handling after this */
    if(dirp != NULL)
    {
        while ((entry = readdir(dirp)) != NULL) {
            if (entry->d_type == DT_DIR) { /* If the entry is a regular file */
                file_count++;
            }
        }

        closedir(dirp);
    }

    return file_count-2; // eliminating direcotires '.' and '..'
#else
    char command[IFT_STR_DEFAULT_SIZE], frame[IFT_STR_DEFAULT_SIZE];
    FILE *f = NULL;

    sprintf(command, "ls -v %s > _tmp_.txt", path);
    system(command);
    
    f = fopen("_tmp_.txt", "r");
    while(!feof(f)){
      fscanf(f, "%s", frame);
      file_count++;
    }
    fclose(f);
    file_count--; // Removing '\n'
    
    system("rm -f _tmp_.txt");

    return file_count;
#endif
}

char* iftFramePath(const char *folder, const  char *name, int frame_id)
{
    char path[IFT_STR_DEFAULT_SIZE], *output_path = NULL;

    sprintf(path, "%s/%0*d/%s", folder, IFT_VIDEO_FOLDER_FRAME_NZEROES,
            frame_id, name);

    output_path = (char*)iftAlloc(strlen(path), sizeof(char));

    strcpy(output_path, path);

    return output_path;
}

iftImage* iftReadFrame(const char *path, char *name, int frame_id)
{
    char *filename = iftFramePath(path, name, frame_id);
    iftImage *frame = NULL;

    frame = iftReadImageByExt(filename);

    if(filename == NULL)
        iftError("Invalid frame file format", "iftReadFrame");

    iftFree(filename);

    return frame;
}

void iftWriteFrame(iftImage *frame, const char *path, char *name, int frame_id)
{
    int success = 0;
    char *pos, ext[10];
    char *filename = iftFramePath(path, name, frame_id);
    char *folder_name = iftFramePath(path,"", frame_id);
    char cmd[16+strlen(folder_name)];

    sprintf(cmd,"mkdir -p %s", folder_name);

    if(system(cmd) == -1) iftError("Command error", "iftWriteFrame");

    pos = strrchr(name,'.');
    pos++;
    sscanf(pos,"%s",ext);

    if(!iftIs3DImage(frame))
    {
        if(strcmp(ext,"ppm") == 0)
        {
            if(iftIsColorImage(frame))
            {
                iftWriteImageP6(frame, filename);
                success = 1;
            }
        }
        else if(strcmp(ext,"pgm") == 0)
        {
            if(iftMaximumValue(frame) > 255)
                iftWriteImageP2(frame, filename);
            else
                iftWriteImageP5(frame, filename);
            success = 1;
        }
    }

    if(success == 0)
        iftError("Invalid frame file format", "iftWriteFrame");

    iftFree(filename);
    iftFree(folder_name);
}

iftImage* iftReadVideoFolderAsVolume(const char* folder_name, int start_frame, int end_frame, char *frame_name)
{
    int i;
    char *pos, ext[10];
    iftImage *volume = NULL, *frame = NULL;

    if(start_frame < 0)
        start_frame = 1;

    if(end_frame < 0)
        end_frame = iftCountNumberOfFrameFolders(folder_name);

    if(end_frame < start_frame)
    {
        iftError("Start frame must be smaller than end frame", "iftReadVideoFolderAsVolume");
    }


    frame = iftReadFrame(folder_name, frame_name, start_frame);

    pos = strrchr(frame_name,'.') + 1;
    sscanf(pos,"%s",ext);

    if(strcmp(ext,"ppm") == 0)
    {
        volume = iftCreateColorImage(frame->xsize, frame->ysize, end_frame - start_frame+1, 8);
    }
    else if(strcmp(ext,"pgm") == 0)
    {
        volume = iftCreateImage(frame->xsize, frame->ysize, end_frame - start_frame+1);
    }
    else
    {
        iftError("Invalid frame file format", "iftReadVideoFolderAsVolume");
    }

    iftPutXYSlice(volume, frame, 0);
    iftDestroyImage(&frame);

    for(i = start_frame+1; i <= end_frame; i++)
    {
        frame = iftReadFrame(folder_name, frame_name, i);

        iftPutXYSlice(volume, frame, i-start_frame);

        iftDestroyImage(&frame);
    }

    return volume;
}

iftImage* iftReadImageFolderAsVolume(const char* folder_name)
{
    int z;
    iftImage *volume = NULL;
    iftFileSet *files = NULL;

    files = iftLoadFileSetFromDirOrCSV(folder_name, 0, true);

    // Verifying if all
    for(z = 0; z < files->n; z++)
    {
        const char *path = files->files[z]->path;
        if(!iftIsImageFile(path) || iftCompareStrings(iftFileExt(path),".scn")) {
            iftError("File %s cannot be read as an image for constructing a 3D volume!", "iftReadImageFolderAsVolume",
                     path);
        }
    }

    for(z = 0; z < files->n; z++)
    {
        const char *path = files->files[z]->path;
        iftImage *slice = iftReadImageByExt(path);

        if(volume == NULL) {
            if(iftIsColorImage(slice))
                volume = iftCreateColorImage(slice->xsize, slice->ysize, files->n, iftImageDepth(slice));
            else
                volume = iftCreateImage(slice->xsize, slice->ysize, files->n);
        }

        if(slice->xsize != volume->xsize || slice->ysize != volume->ysize)
            iftError("The XY dimensions of the current image being loaded (%s) differs from the first image selected" \
            "for initializing the volume! First image: %d %d %d, current image: %d %d %d\n",
                     "iftReadImageFolderAsVolume",
                     path, volume->xsize, volume->ysize, volume->zsize,
                     slice->xsize, slice->ysize, slice->zsize);

        if(iftIsColorImage(slice) != iftIsColorImage(volume)) {
            iftError(
                    "Both the current slice (%s) image being loaded and the first one must either be colored or grayscaled!"
                    "iftReadImageFolderAsVolume", path);
        }

        iftPutXYSlice(volume, slice, z);

        iftDestroyImage(&slice);
    }

    return volume;
}


void iftWriteVolumeAsVideoFolder(iftImage *video, const char *folder_name, char *frame_name)
{
    int i;
    char *pos, ext[10];
    iftImage *frame = NULL;

    pos = strrchr(frame_name,'.');
    pos++;
    sscanf(pos,"%s",ext);

    if((strcmp(ext,"ppm") == 0 && !iftIsColorImage(video))
       || (strcmp(ext,"pgm") == 0 && iftIsColorImage(video)))
    {
        iftError("Invalid frame file format", "iftWriteVolumeAsVideoFolder");
    }

    for(i = 0; i < video->zsize; i++)
    {
        frame = iftGetXYSlice(video, i);
        iftWriteFrame(frame, folder_name, frame_name, i+1);

        iftDestroyImage(&frame);
    }
}

void iftConvertVideoFileToVideoFolder(const char *video_path, const char *output_folder, int rotate)
{
    iftConvertVideoFramesToImages(video_path, "_tmp", "frame", "ppm", rotate);

    iftStoreFramesInVideoFolder("_tmp","ppm", output_folder, "frame");

    if(system("rm -rf _tmp") == -1) iftError("Command error", "iftConvertVideoFileToVideoFolder");
}

void iftConvertVideoFramesToImages(const char *video_path, const char *output_folder, const char *frame_name, const char *extension, int rotate)
{
    char cmd[65536];
    char output_fname[65536];

    sprintf(output_fname,"%s%%0%dd.%s",frame_name, IFT_VIDEO_FOLDER_FRAME_NZEROES, extension);

    sprintf(cmd,"mkdir -p %s", output_folder);
    if(system(cmd) == -1) iftError("Command error", "iftConvertVideoFramesToImages");

    if(rotate < 0)
        sprintf(cmd, "ffmpeg -i %s -f image2 %s", video_path, output_fname);
    else
        sprintf(cmd, "ffmpeg -i %s -f image2 -vf \"transpose=%d\" %s", video_path,
                rotate, output_fname);

    if(system(cmd) == -1) iftError("Command error", "iftConvertVideoFramesToImages");

    sprintf(cmd, "mv %s*.%s %s",frame_name, extension, output_folder);
    if(system(cmd) == -1) iftError("Command error", "iftConvertVideoFramesToImages");
}

void iftStoreFramesInVideoFolder(const char *input_folder, const char *input_extension, const char *output_folder, const char *output_filename)
{
    char cmd[65536];
    char expr[] = "j=1;" \
                "for i in $(ls %s/*.%s); do " \
                    "number=`printf \"%%0%dd\" $j`;"
            "mkdir -p %s/$number; " \
                    "mv $i %s/$number/%s.%s;" \
                    "j=`expr $j + 1`; " \
                "done";
    sprintf(cmd, expr, input_folder, input_extension, IFT_VIDEO_FOLDER_FRAME_NZEROES,
            output_folder, output_folder, output_filename,  input_extension);

    if(system(cmd) == -1) iftError("Command error", "iftStoreFramesInVideoFolder");
}

// ---------- iftVideo.c end
// ---------- iftFileSet.c start
int _iftCmpFilesSortFileSet(const void *a, const void *b) {
    iftFile **f1 = (iftFile**) a;
    iftFile **f2 = (iftFile**) b;
    
    return strcmp((*f1)->path, (*f2)->path);
}

void _iftGetFilesFromDirRec(iftDir *dir, iftSList *SL) {
    for (long i = 0; i < dir->nfiles; i++)
        iftInsertSListIntoTail(SL, dir->files[i]->path);
    
    for (long i = 0; i < dir->nsubdirs; i++)
        _iftGetFilesFromDirRec(dir->subdirs[i], SL);
}

iftFileSet *iftLoadFileSetFromDirOrCSV(const char *file_entry, long hier_levels, bool sort_pathnames) {
    iftFileSet *fset = NULL;
    if (iftDirExists(file_entry))
        fset = iftLoadFileSetFromDir(file_entry, hier_levels); // it also returns a sorted list
    else {
        char *lower_file_entry = iftLowerString(file_entry);
        
        if (iftFileExists(file_entry) && iftEndsWith(lower_file_entry, ".csv")) {
            fset = iftLoadFileSetFromCSV(file_entry, sort_pathnames);
            iftFree(lower_file_entry);
            
            // if (sort_pathnames)
            //     iftSortFileSet(fset);
        }
        else
            iftError("Invalid File Entry: %s\nIt is neither a directory nor a CSV file",
                     "iftLoadFileSetFromDirOrCSV", file_entry);
    }
    
    return fset;
}

iftFileSet *iftLoadFileSetFromCSV(const char *csv_pathname, bool sort_pathnames) {
    iftCSV *csv   = iftReadCSV(csv_pathname, ',');
    long nfiles = csv->nrows*csv->ncols;
    
    iftFileSet *farr = iftCreateFileSet(nfiles);

    #if IFT_OMP
    #pragma omp parallel for
    #endif
    for (long i = 0; i < csv->nrows; i++)
        for (long j = 0; j < csv->ncols; j++) {
            long p = j + i*csv->ncols;
            char *aux = iftExpandUser(csv->data[i][j]);
            farr->files[p] = iftCreateFile(aux);
            iftFree(aux);
        }
    iftDestroyCSV(&csv);
    
    if (sort_pathnames)
        iftSortFileSet(farr);
    
    return farr;
}

iftFileSet *iftLoadFileSetFromDir(const char *dir_pathname, long hier_level) {
    if (dir_pathname == NULL)
        iftError("Directory \"%s\" is NULL", "iftLoadFileSetFromDir");
    if (!iftDirExists(dir_pathname))
        iftError("Directory \"%s\" does not exist", "iftLoadFileSetFromDir", dir_pathname);
    
    iftDir *dir        = NULL; // loads all files from all dir levels
    iftFileSet *farr   = NULL;
    iftSList *SL       = NULL;
    iftSNode *snode    = NULL, *tnode = NULL;
    
    dir = iftLoadDir(dir_pathname, hier_level);
    SL  = iftCreateSList();
    
    // Gets all files in the entire from the directory and its subdirs
    _iftGetFilesFromDirRec(dir, SL);
    
    /**** Converts the String List into a File Array ****/
    farr        = iftCreateFileSet(SL->n);
    
    // copies each node and destroys it
    snode = SL->head;
    long i = 0;
    while (snode!= NULL) {
        farr->files[i] = iftCreateFile(snode->elem);
        i++;
        
        tnode = snode;
        snode = snode->next;
        iftFree(tnode->elem);
        iftFree(tnode);
    }
    SL->head = SL->tail = NULL;
    iftDestroySList(&SL);
    /*******************/
    
    iftDestroyDir(&dir);
    
    return farr;
}
iftFileSet *iftCreateFileSet(long nfiles) {
    iftFileSet *farr = NULL;
    farr        = (iftFileSet *) iftAlloc(1, sizeof(iftFileSet));
    farr->files = (iftFile**) iftAlloc(nfiles, sizeof(iftFile*));
    farr->n = nfiles;
    
    for(long i = 0; i < farr->n; i++)
        farr->files[i] = NULL;
    
    return farr;
}

iftFileSet *iftLoadFileSetFromDirByRegex(const char *dir_pathname, const char *regex, bool sort_pathnames) {
    if (dir_pathname == NULL)
        iftError("Dir's pathname is NULL", "iftLoadFilesFromDirByRegex");
    if (regex == NULL)
        iftError("Regex is NULL", "iftLoadFilesFromDirByRegex");
    if (!iftDirExists(dir_pathname))
        iftError("Directory \"%s\" does not exist", "iftReadFilesFromdDirectory", dir_pathname);
    
    iftDir *dir      = iftLoadFilesFromDirByRegex(dir_pathname, regex);
    iftFileSet *farr = iftCreateFileSet(dir->nfiles);
    
    for (long i = 0; i < dir->nfiles; i++)
        farr->files[i] = iftCopyFile(dir->files[i]);
    
    if (sort_pathnames)
        iftSortFileSet(farr);
    
    iftDestroyDir(&dir);
    
    return farr;
}

void iftDestroyFileSet(iftFileSet **farr) {
    if (farr != NULL) {
        iftFileSet *faux = *farr;
        
        if (faux != NULL) {
            if (faux->files != NULL) {
                for (long i = 0; i < faux->n; i++)
                    iftDestroyFile(&(faux->files[i]));
                iftFree(faux->files);
            }
            iftFree(faux);
            *farr = NULL;
        }
    }
}

void iftSortFileSet(iftFileSet *files) {
    qsort(files->files, files->n, sizeof(iftFile*), _iftCmpFilesSortFileSet);
}
// ---------- iftFileSet.c end
// ---------- iftCSV.c start
bool _iftHasCSVHeader(const char *csv_pathname, char separator) {
    bool has_header = false;
    
    FILE *fp = fopen(csv_pathname, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "_iftHasCSVHeader", csv_pathname);
    
    // reads the first line for checking if it is a CSV header
    char *line = iftGetLine(fp);
    if (line != NULL) {
        bool is_first_line_ok = true;
        
        // if all columns of first line have letters in the beginning of the string, the row can be the header
        char strSeparator[2] = {separator, '\0'};
        iftSList *SL = iftSplitString(line, strSeparator);
        while (!iftIsSListEmpty(SL)) {
            char *column = iftRemoveSListHead(SL);
            
            if (!iftRegexMatch(column, "^[a-zA-Z]+.*$", separator)) {
                is_first_line_ok = false;
                iftFree(column);
                break;
            }
            iftFree(column);
        }
        iftDestroySList(&SL);
        
        if (is_first_line_ok) {
            iftFree(line);
            
            line = iftGetLine(fp);
            if (line != NULL) {
                iftSList *SL = iftSplitString(line, strSeparator);
                iftFree(line);
                
                while (SL->n != 0) {
                    char *column = iftRemoveSListHead(SL);
                    
                    // if at least one column of the second row is a number (integer or real)
                    // the first row is a header
                    if (iftRegexMatch(column, "^[0-9]+(.[0-9]+)?$", separator)) {
                        iftFree(column);
                        has_header = true;
                        break;
                    }
                    iftFree(column);
                }
                iftDestroySList(&SL);
            }
        }
    }
    fclose(fp);
    
    return has_header;
}

void _iftCountNumOfRowsAndColsFromCSVFile(const char *csv_pathname, long *nrows, long *ncols, char separator) {
    char strSeparator[2] = {separator, '\0'};
    
    FILE *fp = fopen(csv_pathname, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "_iftCountNumOfRowsAndColsFromCSVFile", csv_pathname);
    
    *nrows = 0;
    *ncols = 0;
    
    // reads the first line from the file to get the number of cols
    iftSList *SL = NULL;
    char *line = iftGetLine(fp);
    
    // gets the number of columns from the first line, because such number must be the same for
    // the entire csv
    if (line != NULL) {
        SL = iftSplitString(line, strSeparator);
        (*nrows)++;
        *ncols = SL->n;
        iftDestroySList(&SL);
    }
    
    iftFree(line);
    line = iftGetLine(fp);
    
    // gets each line of the file
    while (line != NULL) {
        SL = iftSplitString(line, strSeparator);
        
        if (*ncols != SL->n)
            iftError("Number of Columns is different in the lines: %d - %d",
                     "_iftCountNumOfRowsAndColsFromCSVFile", *ncols, SL->n);
        
        iftDestroySList(&SL);
        (*nrows)++;
        
        iftFree(line);
        line = iftGetLine(fp);
    }
    
    fclose(fp);
}

iftCSV *_iftCreateCSVWithoutStringAllocation(long nrows, long ncols) {
    iftCSV *csv = (iftCSV*) iftAlloc(1, sizeof(iftCSV));
    
    csv->nrows = nrows;
    csv->ncols = ncols;
    
    // allocates the CSV string matrix
    csv->data = (char***) iftAlloc(nrows, sizeof(char**));
    for (long i = 0; i < nrows; i++) {
        csv->data[i] = (char**) iftAlloc(ncols, sizeof(char*));
    }
    
    return csv;
}

iftCSV *iftReadCSV(const char *csv_pathname, const char separator) {
    if (!iftFileExists(csv_pathname))
        iftError("The CSV file pathname \"%s\" does not exists!", "iftReadCSV", csv_pathname);
    
    char strSeparator[2] = {separator, '\0'};
    
    bool has_header = _iftHasCSVHeader(csv_pathname, separator);
    
    long nrows, ncols;
    _iftCountNumOfRowsAndColsFromCSVFile(csv_pathname, &nrows, &ncols, separator);
    if (has_header)
        nrows--;
    
    iftCSV *csv = _iftCreateCSVWithoutStringAllocation(nrows, ncols);
    
    FILE *fp = fopen(csv_pathname, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "_iftCountNumOfRowsAndColsFromCSVFile", csv_pathname);
    
    // copies the values from the CSV file
    iftSList *SL = NULL;
    char *line = iftGetLine(fp);
    
    if (has_header) {
        csv->header = iftAlloc(csv->ncols, sizeof(char*));
        
        SL = iftSplitString(line, strSeparator);
        
        for (long j = 0; j < csv->ncols; j++) {
            csv->header[j] = iftRemoveSListHead(SL); // just points to string
            // removes the '\n' and '\r' from the paths
            iftRightTrim(csv->header[j], '\n');
            iftRightTrim(csv->header[j], '\r');
        }
        iftDestroySList(&SL);
        
        iftFree(line);
        line = iftGetLine(fp);
    }
    
    long i = 0;
    while (line != NULL) {
        SL = iftSplitString(line, strSeparator);
        
        for (long j = 0; j < csv->ncols; j++) {
            csv->data[i][j] = iftRemoveSListHead(SL); // just points to string
            // removes the '\n' and '\r' from the paths
            iftRightTrim(csv->data[i][j], '\n');
            iftRightTrim(csv->data[i][j], '\r');
        }
        i++;
        iftDestroySList(&SL);
        
        iftFree(line);
        line = iftGetLine(fp);
    }
    
    fclose(fp);
    
    return csv;
}

void iftDestroyCSV(iftCSV **csv) {
    iftCSV *csv_aux = *csv;
    
    if (csv_aux != NULL) {
        if (csv_aux->data != NULL) {
            // deallocates the CSV string matrix
            for (long i = 0; i < csv_aux->nrows; i++) {
                if (csv_aux->data[i] != NULL)
                    for (long j = 0; j < csv_aux->ncols; j++) {
                        iftFree(csv_aux->data[i][j]);
                    }
                iftFree(csv_aux->data[i]);
            }
        }
        iftFree(csv_aux->data);
        
        if (csv_aux->header != NULL) {
            for (int c = 0; c < csv_aux->ncols; c++)
                iftFree(csv_aux->header[c]);
            iftFree(csv_aux->header);
        }
        iftFree(csv_aux);
        *csv = NULL;
    }
}
// ---------- iftCSV.c end

//===========================================================================//
// ADDED BY FELIPE
//===========================================================================//
void iftConvertNewBitDepth(iftImage **img, int new_depth)
{
    #if IFT_DEBUG
    assert(img != NULL && *img != NULL);
    assert(new_depth > 0);
    #endif
    int old_depth, old_norm, new_norm;
    float ratio;

    old_depth = iftImageDepth(*img);
    old_norm = iftMaxImageRange(old_depth);
    new_norm = iftMaxImageRange(new_depth);
    ratio = new_norm / (float) old_norm;

    #if IFT_OMP
    #pragma omp parallel for
    #endif
    for(int p = 0; p < (*img)->n; ++p)
    {
        iftColor rgb, new_ycbcr;

        if(iftIsColorImage(*img) == true)
        {
            iftColor old_ycbcr;

            old_ycbcr.val[0] = (*img)->val[p];
            old_ycbcr.val[1] = (*img)->Cb[p];
            old_ycbcr.val[2] = (*img)->Cr[p];

            rgb = iftYCbCrtoRGB(old_ycbcr, old_norm);
        }
        else rgb.val[0] = rgb.val[1] = rgb.val[2] = (*img)->val[p];

        rgb.val[0] *= ratio; rgb.val[1] *= ratio; rgb.val[2] *= ratio;
        new_ycbcr = iftRGBtoYCbCr(rgb, new_norm);

        (*img)->val[p] = new_ycbcr.val[0];
        if(iftIsColorImage(*img) == true)
        {
            (*img)->Cb[p] = new_ycbcr.val[1];
            (*img)->Cr[p] = new_ycbcr.val[2];
        }
    }
}

void iftWriteVolumeAsSingleVideoFolder
(const iftImage *video, const char *path)
{
    int i;
    iftImage *frame = NULL;
    const char *EXT = iftFileExt(path);
    char *base = iftRemoveSuffix(path, EXT);

    for(i = 0; i < video->zsize; i++)
    {
        char mod_path[IFT_STR_DEFAULT_SIZE];

        sprintf(mod_path, "%s_%0*d%s", base, IFT_VIDEO_FOLDER_FRAME_NZEROES, 
                                       i, EXT);

        frame = iftGetXYSlice(video, i);
        iftWriteImageByExt(frame, mod_path);

        iftDestroyImage(&frame);
    }
    free(base);
}

iftBMap *iftGetBorderMap
(const iftImage *label_img)
{
  #if IFT_DEBUG
  assert(label_img != NULL);
  #endif
  iftAdjRel *A;
  iftBMap *border_map;

  if(iftIs3DImage(label_img) == true) A = iftSpheric(1.0);
  else A = iftCircular(1.0);

  border_map = iftCreateBMap(label_img->n);

  #if IFT_OMP //-------------------------------------------------------------//
  #pragma omp parallel for
  #endif //------------------------------------------------------------------//
  for(int p = 0; p < label_img->n; ++p)
  {
    bool is_border;
    int i;
    iftVoxel p_vxl;

    is_border = false;
    p_vxl = iftGetVoxelCoord(label_img, p);

    i = 0;
    while(is_border == false && i < A->n)
    {
      iftVoxel adj_vxl;

      adj_vxl = iftGetAdjacentVoxel(A, p_vxl, i);

      if(iftValidVoxel(label_img, adj_vxl) == true)
      {
        int adj_idx;

        adj_idx = iftGetVoxelIndex(label_img, adj_vxl);

        if(label_img->val[p] != label_img->val[adj_idx]) 
          is_border = true;
      }

      ++i;
    }

    if(is_border == true) iftBMapSet1(border_map, p);
  }
  iftDestroyAdjRel(&A);

  return border_map;
}