/*
 * @file ift.h
 * @brief Bundle of relevant data types and functions from the IFT library.
 */

#ifndef _IFT_H_
#define _IFT_H_

#ifdef __cplusplus
extern "C" {
#endif

// ---------- iftGQueue.h start 

#define MINVALUE   0
#define MAXVALUE   1
#define FIFOBREAK  0
#define LIFOBREAK  1
#define IFT_QSIZE      65536
#define iftSetTieBreak(a,b) a->C.tiebreak=b
#define iftSetRemovalPolicy(a,b) a->C.removal_policy=b

typedef struct ift_gqnode {
    int  next;
    int  prev;
    char color;
} iftGQNode;

typedef struct ift_gdoublylinkedlists {
    iftGQNode *elem;
    int nelems;
    int *value;
} iftGDoublyLinkedLists;

typedef struct ift_gcircularqueue {
    int  *first;
    int  *last;
    int  nbuckets;
    int  minvalue;
    int  maxvalue;
    char tiebreak;
    char removal_policy;
} iftGCircularQueue;

typedef struct ift_gqueue {
    iftGCircularQueue C;
    iftGDoublyLinkedLists L;
} iftGQueue;

iftGQueue *iftCreateGQueue(int nbuckets, int nelems, int *value);
void   iftDestroyGQueue(iftGQueue **Q);
int    iftEmptyGQueue(iftGQueue *Q);
void   iftInsertGQueue(iftGQueue **Q, int elem);
int    iftRemoveGQueue(iftGQueue *Q);
void   iftRemoveGQueueElem(iftGQueue *Q, int elem);
void iftResetGQueue(iftGQueue *Q);
iftGQueue *iftGrowGQueue(iftGQueue **Q, int nbuckets);

// ---------- iftGQueue.h end 
// ---------- iftBasicDataTypes.h start 

#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <regex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define IFT_INFINITY_INT       INT_MAX
#define IFT_INFINITY_INT_NEG   INT_MIN
#define IFT_INFINITY_LONG      LONG_MAX
#define IFT_INFINITY_LONG_NEG  LONG_MIN
#define IFT_INFINITY_ULONG     ULONG_MAX
#define IFT_INFINITY_FLT       FLT_MAX
#define IFT_INFINITY_FLT_NEG  -FLT_MAX
#define IFT_INFINITY_DBL       DBL_MAX
#define IFT_INFINITY_DBL_NEG  -DBL_MAX
#define IFT_INFINITY_LDBL      LDBL_MAX
#define IFT_INFINITY_LDBL_NEG -LDBL_MAX

#define IFT_STR_DEFAULT_SIZE 4096
#define IFT_ULONG_NIL IFT_INFINITY_ULONG

typedef struct ift_dict iftDict;
typedef struct ift_char_array iftCharArray;
typedef struct ift_int_array iftIntArray;
typedef struct ift_dbl_array iftDblArray;
typedef struct ift_str_array iftStrArray;
typedef struct ift_int_matrix iftIntMatrix;
typedef struct ift_matrix iftMatrix;
typedef struct ift_str_matrix iftStrMatrix;

typedef struct timeval timer;
typedef unsigned char  uchar;
typedef unsigned short ushort;
typedef unsigned int   uint;
typedef unsigned long  ulong;
#ifndef  __cplusplus
typedef long long llong;
typedef unsigned long long ullong;
#endif

typedef enum ift_cdata_type {
    IFT_UNTYPED, IFT_BOOL_TYPE, IFT_CHAR_TYPE, IFT_UCHAR_TYPE, IFT_STR_TYPE, IFT_SHORT_TYPE, IFT_USHORT_TYPE,
    IFT_INT_TYPE, IFT_UINT_TYPE, IFT_LONG_TYPE, IFT_ULONG_TYPE, IFT_FLT_TYPE, IFT_DBL_TYPE, IFT_INT_ARRAY_TYPE,
    IFT_DBL_ARRAY_TYPE, IFT_STR_ARRAY_TYPE, IFT_INT_MATRIX_TYPE, IFT_DBL_MATRIX_TYPE, IFT_STR_MATRIX_TYPE,
    IFT_DICT_TYPE, IFT_PTR_TYPE
} iftCDataType;

typedef struct ift_band {
    float *val;
} iftBand;

typedef struct ift_vector {
    float x, y, z;
} iftVector, iftPoint;

typedef struct ift_voxel {
    int x, y, z;
} iftVoxel;

typedef struct ift_size {
    int x, y, z;
} iftSize;

typedef struct ift_interval {
    float begin;
    float end;
} iftInterval;

typedef struct ift_image_domain {
    int xsize;
    int ysize;
    int zsize;
} iftImageDomain;

typedef struct ift_voxel_size {
    float dx;
    float dy;
    float dz;
} iftVoxelSize;

typedef struct ift_complex {
    double r;
    double i;
} iftComplex;

typedef struct ift_gval {
    iftCDataType type;
    union {
        bool   bool_val;
        char   char_val;
        uchar  uchar_val;
        char*  str_val;
        long   long_val;
        ulong  ulong_val;
        double dbl_val;
        iftIntArray *int_array_val;
        iftDblArray* dbl_array_val;
        iftStrArray* str_array_val;
        iftIntMatrix* int_matrix_val;
        iftMatrix* dbl_matrix_val;
        iftStrMatrix* str_matrix_val;
        iftDict* dict_val;
        void*  ptr_val;
    };
} iftGVal;

#define iftBoolAsString(b) b ? "True" : "False"

#define iftCreateGVal(VAL) _Generic((VAL), \
    bool:        iftInitBoolGVal, \
    char:        iftInitCharGVal, \
    uchar:       iftInitUCharGVal, \
    char*:       iftInitStrGVal, \
    const char*: iftInitStrGVal, \
    short:       iftInitLongGVal, \
    ushort:      iftInitULongGVal, \
    int:         iftInitLongGVal, \
    uint:        iftInitULongGVal, \
    long:        iftInitLongGVal, \
    ulong:       iftInitULongGVal, \
    float:       iftInitDblGVal, \
    double:      iftInitDblGVal, \
    iftIntArray*: iftInitIntArrayGVal, \
    iftDblArray*: iftInitDblArrayGVal, \
    iftStrArray*: iftInitStrArrayGVal, \
    iftIntMatrix*: iftInitIntMatrixGVal, \
    iftMatrix*: iftInitDblMatrixGVal, \
    iftStrMatrix*: iftInitStrMatrixGVal, \
    iftDict*: iftInitDictGVal, \
    void*:       iftInitPtrGVal, \
    iftGVal:     iftCopyGVal, \
    default:     iftInitPtrGVal \
    )(VAL)

iftGVal iftInitBoolGVal(bool val);
iftGVal iftInitCharGVal(char val);
iftGVal iftInitUCharGVal(uchar val);
iftGVal iftInitStrGVal(const char *val);
iftGVal iftInitLongGVal(long val);
iftGVal iftInitULongGVal(ulong val);
iftGVal iftInitDblGVal(double val);
iftGVal iftInitIntArrayGVal(iftIntArray *array);
iftGVal iftInitDblArrayGVal(iftDblArray *array);
iftGVal iftInitStrArrayGVal(iftStrArray *array);
iftGVal iftInitIntMatrixGVal(iftIntMatrix *mat);
iftGVal iftInitDblMatrixGVal(iftMatrix *mat);
iftGVal iftInitStrMatrixGVal(iftStrMatrix *mat);
iftGVal iftInitDictGVal(iftDict *dict);
iftGVal iftInitPtrGVal(void *val);
void iftFreeGVal(iftGVal gv);
const char *iftCDataTypeToString(iftCDataType datatype);
char *iftGValToString(iftGVal gval);
bool iftCompareGVal(iftGVal val1, iftGVal val2);
iftGVal iftCopyGVal(iftGVal gval);
bool iftGetBoolVal(iftGVal gval);
char iftGetCharVal(iftGVal gval);
uchar iftGetUCharVal(iftGVal gval);
char *iftGetStrVal(iftGVal gval);
const char *iftGetConstStrVal(iftGVal gval);
long iftGetLongVal(iftGVal gval);
ulong iftGetULongVal(iftGVal gval);
double iftGetDblVal(iftGVal gval);
iftIntArray *iftGetIntArrayVal(iftGVal gval);
iftDblArray *iftGetDblArrayVal(iftGVal gval);
iftStrArray *iftGetStrArrayVal(iftGVal gval);
iftIntMatrix *iftGetIntMatrixVal(iftGVal gval);
iftMatrix *iftGetDblMatrixVal(iftGVal gval);
iftStrMatrix *iftGetStrMatrixVal(iftGVal gval);
iftDict *iftGetDictVal(iftGVal gval);
void *iftGetPtrVal(iftGVal gval);
void iftCopyVoxel(iftVoxel *src, iftVoxel *dst);

// ---------- iftBasicDataTypes.h end
// ---------- iftIntArray.h start

struct ift_int_array {
    long n;
    int *val;
};

iftIntArray *iftCreateIntArray(long n);
void iftDestroyIntArray(iftIntArray **iarr);
void iftShuffleIntArray(int* array, int n);

// ---------- iftIntArray.h end
// ---------- iftFloatArray.h start 

typedef struct ift_flt_array {
    long n;
    float *val;
} iftFloatArray;

iftFloatArray *iftCreateFloatArray(long n);
void iftDestroyFloatArray(iftFloatArray **darr);

// ---------- iftFloatArray.h end
// ---------- iftDblArray.h start

struct ift_dbl_array {
    /** Number of elements. */
    long n;
    /** Array of double values. */
    double *val;
};

iftDblArray *iftCreateDblArray(long n);
iftDblArray *iftCopyDblArray(const double* array, long n);
void iftDestroyDblArray(iftDblArray **darr);

// ---------- iftDblArray.h end
// ---------- iftStrArray.h start

struct ift_str_array {
    long n;
    char **val;
};

iftStrArray *iftCreateStrArray(long n);
iftStrArray *iftCopyStrArray(char **str_arr, long n);
void iftDestroyStrArray(iftStrArray **sarr);

// ---------- iftStrArray.h end
// ---------- iftIntMatrix.h start

struct ift_int_matrix {
    int *val;
    int ncols, nrows, *tbrow;
    int n;
};

iftIntMatrix *iftCreateIntMatrix(int ncols, int nrows);
iftIntMatrix *iftCopyIntMatrix(int *val, int nrows, int ncols);
void iftDestroyIntMatrix(iftIntMatrix **iM);

// ---------- iftIntMatrix.h end 
// ---------- iftStrMatrix.h start 
// ---------- iftStrMatrix.h end
// ---------- iftCommon.h start 

#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include <float.h>
#include <libgen.h>
#include <limits.h>
#include <math.h>
#include <mm_malloc.h>
#include <regex.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#define IFT_INFINITY_FLT_NEG  -FLT_MAX

#define IFT_PI      3.14159265358979323846
#define IFT_WHITE       0
#define IFT_GRAY        1
#define IFT_BLACK       2
#define IFT_NIL        -1
#define IFT_INCREASING  1
#define IFT_DECREASING  0
#define IFT_EPSILON     1E-07

#ifndef iftMax
#define iftMax(x,y) (((x) > (y))?(x):(y))
#endif

#ifndef iftMin
#define iftMin(x,y) (((x) < (y))?(x):(y))
#endif

#define iftRound(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))

#define iftSmoothEuclideanDistance(squared_dist) ((fabs(squared_dist-1)<1.0e-6)? 0.9016 : (fabs(squared_dist-2)<1.0e-6)? 1.289 : 1.615 )
#define iftSquaredVoxelDistance(u,v) ((u.x-v.x)*(u.x-v.x)+(u.y-v.y)*(u.y-v.y)+(u.z-v.z)*(u.z-v.z))

int iftFArgmax(const float *x, int n);
float iftFeatDistance(float *A, float *B, int n);
int iftRandomInteger (int low, int high);
void iftSkipComments(FILE *fp);
double iftLog(double val, double base);
long iftNormalizationValue(long maxval);
int iftAlmostZero(double x);
void iftRandomSeed(unsigned int);
int iftSafeMod(int a, int n);
bool iftIsPrime(long n);
void iftNormalizeFeatures(float *feats, int nelems);
timer *iftTic(void);
timer *iftToc(void);
float iftCompTime(timer *tic, timer *toc);

// ---------- iftCommon.h end
// ---------- iftAdjacency.h start 

typedef struct ift_adjrel {
    int *dx, *dy, *dz;
    int n;
} iftAdjRel;

iftAdjRel *iftCreateAdjRel(int n);
void iftDestroyAdjRel(iftAdjRel **A);
iftAdjRel *iftSpheric(float r);
iftAdjRel *iftCircular(float r);
iftAdjRel *iftCopyAdjacency(const iftAdjRel *A);
iftVoxel iftGetAdjacentVoxel(const iftAdjRel *A, iftVoxel u, int adj);
void iftMaxAdjShifts(const iftAdjRel *A, int *dx, int *dy, int *dz);
iftAdjRel *iftAdjacencyBoundaries(const iftAdjRel *A, const iftAdjRel *B);

// ---------- iftAdjacency.h end
// ---------- iftCSV.h start 

typedef struct ift_csv {
    long nrows;
    long ncols;
    char **header;
    char ***data;
} iftCSV;

iftCSV *iftCreateCSV(long nrows, long ncols);
void iftDestroyCSV(iftCSV **csv);

// ---------- iftCSV.h end
// ---------- iftColor.h start 

#define WHITEPOINT_X	0.950456
#define WHITEPOINT_Y	1.0
#define WHITEPOINT_Z	1.088754

#define LABF(t)	\
	((t >= 8.85645167903563082e-3) ? \
	pow(t,0.333333333333333) : (841.0/108.0)*(t) + (4.0/29.0))

typedef enum ift_color_space {
    YCbCr_CSPACE,
    YCbCrNorm_CSPACE,
    RGB_CSPACE,
    RGBNorm_CSPACE,
    GRAY_CSPACE,
    GRAYNorm_CSPACE,
    WEIGHTED_YCbCr_CSPACE,
    LAB_CSPACE,
    LABNorm_CSPACE,
    LABNorm2_CSPACE,
    HSV_CSPACE
} iftColorSpace;

typedef struct ift_color {
  int val[3];
  float alpha;
} iftColor;

typedef struct ift_colortable {
  iftColor *color;
  int ncolors;
} iftColorTable;

typedef struct ift_fcolor {
	float val[3];
} iftFColor;


iftColorTable *iftCreateColorTable(int ncolors);
void iftDestroyColorTable(iftColorTable **ctb);
iftColor  iftRGBtoYCbCrBT2020(iftColor cin, const int rgbBitDepth, const int yCbCrBitDepth);
iftColor iftRGBtoYCbCr(iftColor cin, int normalization_value);
iftColor iftYCbCrtoRGB(iftColor cin, int normalization_value);
iftColor  iftYCbCrBT2020toRGB(iftColor cin, const int yCbCrBitDepth, const int rgbBitDepth);
static inline iftFColor iftRGBtoLab(iftColor rgb, int normalization_value)
{
        /* with no normalization. */
    // L in [0, 99.998337]
    // a in [-86.182236, 98.258614]
    // b in [-107.867744, 94.481682]

    //RGB to XYZ
    float R = rgb.val[0]/(float)normalization_value;
    float G = rgb.val[1]/(float)normalization_value;
    float B = rgb.val[2]/(float)normalization_value;

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
iftFColor iftRGBtoLabNorm(iftColor rgb, int normalization_value);
static inline iftFColor iftRGBtoLabNorm2(iftColor rgb, int normalization_value)
{
    /* get lab values*/
    iftFColor lab=iftRGBtoLab(rgb,normalization_value);
    /*normalize each value*/
    lab.val[0]=lab.val[0]/99.998337f;
    lab.val[1]=(lab.val[1]+86.182236f)/(86.182236f+98.258614f);
    lab.val[2]=(lab.val[2]+107.867744f)/(107.867744f+94.481682f);

    return lab;
}
iftColor  iftRGBtoHSV(iftColor cin, int normalization_value);
iftColor  iftHSVtoRGB(iftColor cin, int normalization_value);

// ---------- iftColor.h end
// ---------- iftDHeap.h start

#define MINVALUE   0 /* define heap to remove node with minimum value */
#define MAXVALUE   1 /* define heap to remove node with maximum value */

#define iftDad(i) ((i - 1) / 2)
#define iftLeftSon(i) (2 * i + 1)
#define iftRightSon(i) (2 * i + 2)
#define iftSetRemovalPolicyDHeap(a,b) a->removal_policy = b

typedef struct ift_dheap {
    double *value;
    char  *color;
    int   *node;
    int   *pos;
    int    last;
    int    n;
    char removal_policy;
} iftDHeap;

iftDHeap *iftCreateDHeap(int n, double *value);
void      iftDestroyDHeap(iftDHeap **H);
char      iftFullDHeap(iftDHeap *H);
char      iftEmptyDHeap(iftDHeap *H);
char      iftInsertDHeap(iftDHeap *H, int pixel);
int       iftRemoveDHeap(iftDHeap *H);
void      iftRemoveDHeapElem(iftDHeap *H, int pixel);
void      iftGoUpDHeap(iftDHeap *H, int i);
void      iftGoDownDHeap(iftDHeap *H, int i);
void      iftResetDHeap(iftDHeap *H);

// ---------- iftDHeap.h end
// ---------- iftFIFO.h start 

typedef struct ift_fifo {
    int *FIFO;
    int n;
    int first;
    int last;
    char *color;
} iftFIFO;

iftFIFO *iftCreateFIFO(int n);
void iftDestroyFIFO(iftFIFO **F);
char iftInsertFIFO(iftFIFO *F, int elem);
int iftRemoveFIFO(iftFIFO *F);
bool iftFullFIFO(iftFIFO *F);
bool iftEmptyFIFO(iftFIFO *F);
void iftResetFIFO(iftFIFO *F);
int iftColorFIFO(iftFIFO *F, int pos);

// ---------- iftFIFO.h end
// ---------- iftFile.h start

#if defined(__WIN32) || defined(__WIN64)
#define IFT_SEP_C "\\"
#else
#define IFT_SEP_C "/"
#endif

typedef struct ift_file {
    char *path;
    int label;
    int sample;
    char *suffix;
    char status;
} iftFile;

bool iftFileExists(const char *pathname);
const char *iftFileExt(const char *pathname);
static inline bool iftPathnameExists(const char *pathname) {
    struct stat buffer;
    return (stat (pathname, &buffer) == 0);
}
char *iftJoinPathnames(long n, ...);
char *iftFilename(const char *pathname, const char *suffix);
void iftDestroyFile(iftFile **file);

// ---------- iftFile.h end
// ---------- iftBoundingBox.h start 

typedef struct ift_bounding_box {
    iftVoxel begin;
    iftVoxel end;
} iftBoundingBox;

// ---------- iftBoundingBox.h end
// ---------- iftImage.h start 

#define iftGetVoxelIndex(s, v) ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)])
#define iftImgVal(img, x, y, z) img->val[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftImgCb(img, x, y, z) img->Cb[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftImgCr(img, x, y, z) img->Cr[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftImgAlpha(img, x, y, z) img->alpha[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftCopyVoxelSize(src, dst) (dst)->dx = (src)->dx; (dst)->dy = (src)->dy; (dst)->dz = (src)->dz;
#define iftValidVoxel(img, v)((v.x >= 0) && (v.x <= ((img)->xsize - 1)) && (v.y >= 0) && (v.y <= ((img)->ysize - 1)) && (v.z >= 0) && (v.z <= ((img)->zsize - 1)))
#define iftGetXCoord(s, p) (((p) % (((s)->xsize)*((s)->ysize))) % (s)->xsize)
#define iftGetYCoord(s, p) (((p) % (((s)->xsize)*((s)->ysize))) / (s)->xsize)
#define iftGetZCoord(s, p) ((p) / (((s)->xsize)*((s)->ysize)))
#define iftImgVoxelVal(i, v) i->val[iftGetVoxelIndex(i, v)]
#define iftImgVoxelCb(i, v) i->Cb[iftGetVoxelIndex(i, v)]
#define iftImgVoxelCr(i, v) i->Cr[iftGetVoxelIndex(i, v)]

typedef struct ift_image {
    int *val;
    ushort *Cb;
    ushort *Cr;
    ushort *alpha;
    int xsize;
    int ysize;
    int zsize;
    float dx;
    float dy;
    float dz;
    int *tby, *tbz;
    int n;
} iftImage;

iftImage *iftReadImageByExt(const char *filename, ...);
iftImage *iftCreateImage(int xsize, int ysize, int zsize);
iftImage *iftCopyImage(const iftImage *img);
void iftDestroyImage(iftImage **img);
void iftWriteImageByExt(const iftImage *img, const char *filename, ...);
static inline bool iftIsColorImage(const iftImage *img) {
    return ((img->Cb != NULL) && (img->Cr != NULL));
}
int iftMaximumValue(const iftImage *img);
int iftMinimumValue(const iftImage *img);
static inline bool iftIs3DImage(const iftImage *img) {
    return (img->zsize > 1);
}
iftVoxel iftGetVoxelCoord(const iftImage *img, int p);
iftImage *iftSelectImageDomain(int xsize, int ysize, int zsize);
iftBoundingBox iftMinBoundingBox(const iftImage *img, iftVoxel *gc_out);
iftImage *iftReadImage(const char *filename, ...);
iftImage* iftReadImagePNG(const char* format, ...);
iftImage* iftReadImageJPEG(const char* format, ...);
iftImage *iftReadImageGZip(const char *filename, ...);
iftImage *iftReadImageP5(const char *filename, ...);
iftImage *iftReadImageP6(const char *filename, ...);
iftImage *iftReadImageP2(const char *filename, ...);
iftImage  *iftCreateImageFromBuffer(int xsize,int ysize,int zsize, int *val);
void iftCopyImageInplace(const iftImage *src, iftImage *dst);
void iftWriteImage(const iftImage *img, const char *filename, ...);
void iftWriteImageGZip(const iftImage *img, const char *filename, ...);
void iftWriteImageP5(const iftImage *img, const char *filename, ...);
void iftWriteImageP6(const iftImage *img, const char *filename, ...);
void iftWriteImageP2(const iftImage *img, const char *filename, ...);
void iftWriteImagePNG(const iftImage* img, const char* format, ...);
void iftWriteImageJPEG(const iftImage* img, const char* format, ...);
int iftMaximumValueInRegion(const iftImage *img, iftBoundingBox bb);
void iftSetImage(iftImage *img, int value);
void  iftSetAlpha(iftImage *img, ushort value);
void iftSetCbCr(iftImage *img, ushort value);
void iftVerifyImageDomains(const iftImage *img1, const iftImage *img2, const char *function);
uchar iftImageDepth(const iftImage* img);
static inline long iftMaxImageRange(uchar img_depth) {
    return (1L << (img_depth)) - 1; // 2^img_depth -1
}
void iftMinMaxValues(const iftImage *img, int *min, int *max);
int iftNumberOfElements(const iftImage *mask);
iftImage *iftImageGradientMagnitude(const iftImage *img, iftAdjRel *Ain);
iftImage *iftCreateImageFromImage(const iftImage *src);
iftImage *iftCreateColorImage(int xsize, int ysize, int zsize, int depth);

// ---------- iftImage.h end
// ---------- iftFImage.h start

#define iftFGetVoxelIndex(s,v) ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)])
#define iftFGetXCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) % (s)->xsize)
#define iftFGetYCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) / (s)->xsize)
#define iftFGetZCoord(s,p) ((p) / (((s)->xsize)*((s)->ysize)))

typedef struct ift_fimage {
  float *val;
  int    xsize,ysize,zsize;
  float  dx,dy,dz;
  int   *tby, *tbz;
  int    n;
} iftFImage;

iftFImage *iftCreateFImage(int xsize,int ysize,int zsize);
void iftDestroyFImage(iftFImage **img);
iftVoxel    iftFGetVoxelCoord(const iftFImage *img, int p);
char iftFValidVoxel(const iftFImage *img, iftVoxel v);

// ---------- iftFImage.h end 
// ---------- iftGraphics.h start

void iftDrawPoint(iftImage *img, iftVoxel u, iftColor YCbCr, iftAdjRel *B, int rangeValue);
void iftDrawBorders(iftImage *img, iftImage *label, iftAdjRel *A, iftColor YCbCr, iftAdjRel *B);
void iftDrawBordersSingleLabel(iftImage *img, iftImage *labelMap, const int label, iftColor YCbCr);

// ---------- iftGraphics.h end
// ---------- iftMatrix.h start 

typedef struct ift_matrix {
    int nrows;
    int ncols;
    long n;
    float *val;
    long *tbrow;
    bool allocated;
} iftMatrix;

struct ift_str_matrix {
    char **val;
    int ncols, nrows, *tbrow;
    int n;
};

iftMatrix *iftCreateMatrix(int ncols, int nrows);
iftMatrix *iftCopyMatrix(const iftMatrix *A);
void iftDestroyMatrix(iftMatrix **M);
void iftDestroyStrMatrix(iftStrMatrix **sM);
iftStrMatrix *iftCopyStrMatrix(char **val, int nrows, int ncols);
iftStrMatrix *iftCreateStrMatrix(int ncols, int nrows); 

#define iftGetMatrixCol(m,i) ((i) % (m)->ncols)
#define iftGetMatrixRow(m,i) ((i) / (m)->ncols)
#define iftGetMatrixIndex(m,c,r) ((c)+(m)->tbrow[(r)])
#define iftMatrixRowPointer(m, r) ((m)->val + iftGetMatrixIndex((m), 0, (r)))
#define iftMatrixElem(m, c, r) (m)->val[iftGetMatrixIndex((m), (c), (r))]

// ---------- iftMatrix.h end
// ---------- iftMImage.h start 

typedef struct ift_mimage {
    iftMatrix *data;
    float **val;
    int xsize,ysize,zsize;
    float dx,dy,dz; 
    int *tby, *tbz;
    int n,m;
} iftMImage;

#define iftMGetVoxelIndex(s,v) ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)])

iftMImage  *iftCreateMImage(int xsize,int ysize,int zsize, int nbands);
void iftDestroyMImage(iftMImage **img);
iftMImage *iftImageToMImage(const iftImage *img, char color_space);
iftImage *iftMImageToImage(const iftMImage *img, int Imax, int band);
iftImage *iftGridSampling(iftMImage *img, iftImage *mask, int nsamples);
iftImage *iftAltMixedSampling(iftMImage *img, iftImage *mask, int nsamples);
static inline bool iftIs3DMImage(const iftMImage *img) {
    return (img->zsize > 1);
}
iftImage   *iftBorderProbImage(iftMImage *img);
iftImage   *iftMImageBasins(const iftMImage *img, iftAdjRel *A);
iftVoxel    iftMGetVoxelCoord(const iftMImage *img, int p);
char iftMValidVoxel(const iftMImage *img, iftVoxel v);

// ---------- iftMImage.h end
// ---------- iftMemory.h start 

int *iftAllocIntArray(long n);
void iftCopyIntArray(int *array_dst, const int *array_src, int nelems);
#ifndef  __cplusplus
long long *iftAllocLongLongIntArray(long n);
void iftCopyLongLongIntArray(long long *array_dst, const long long *array_src, int nelems);
#endif
float *iftAllocFloatArray(long n);
void iftCopyFloatArray(float *array_dst, float *array_src, int nelems);
double *iftAllocDoubleArray(long n);
void iftCopyDoubleArray(double *array_dst, double *array_src, int nelems);
ushort *iftAllocUShortArray(long n);
uchar *iftAllocUCharArray(long n);
char *iftAllocCharArray(long n);
char *iftAllocString(long n);
#define iftSwap(x, y) do { __typeof__(x) _IFT_SWAP_ = x; x = y; y = _IFT_SWAP_; } while (0)

// ---------- iftMemory.h end
// ---------- iftMetrics.h start

typedef struct ift_error_classification {
	int tp;
	int fp;
	int fn;
	int tn;
} iftErrorClassification;

float iftFScoreGivenErrors(iftErrorClassification* error);
float iftPrecisionGivenErrors(iftErrorClassification* error);
float iftRecallGivenErrors(iftErrorClassification* error);

// ---------- iftMetrics.h end 
// ---------- iftSegmentation.h start 

float iftBoundaryRecall(iftImage *gt, iftImage *border, float tolerance_dist);
float iftUnderSegmentation(iftImage *gt_image, iftImage *label_image);
iftImage *iftBorderImage(const iftImage *label, bool get_margins);
iftImage* iftBorderImageToLabelImage(iftImage* border);
iftImage  *iftWaterGray(iftImage *basins, iftImage *marker, iftAdjRel *A);
iftFImage *iftSmoothWeightImage(const iftImage *basins, float beta);
iftFImage *iftWeightNormFactor(const iftFImage *weight, iftAdjRel *A);
float *iftFScoreMultiLabel (iftImage *mlabel, iftImage *mlgt, int number_of_objects);
float iftCompactness2D(iftImage *label);
float iftTopologyMeasure(iftImage *label);
float iftFScoreError(iftImage *bin, iftImage *gt);
iftErrorClassification iftSegmentationErrors(iftImage* gt_image, iftImage* cl_image);

// ---------- iftSegmentation.h end
// ---------- iftSet.h start 

typedef struct ift_set {
    int elem;
    struct ift_set *next;
} iftSet;

void iftInsertSet(iftSet **S, int elem);
int iftRemoveSet(iftSet **S);
void    iftRemoveSetElem(iftSet **S, int elem);
void    iftDestroySet(iftSet **S);
iftSet* iftSetUnion(iftSet *S1,iftSet *S2);
iftSet* iftSetConcat(iftSet *S1,iftSet *S2);
char    iftUnionSetElem(iftSet **S, int elem);
void    iftInvertSet(iftSet **S);
int 	iftSetSize(const iftSet* S);
iftSet* iftSetCopy(iftSet* S);
int     iftSetHasElement(iftSet *S, int elem);

// ---------- iftSet.h end
// ---------- iftSList.h start

typedef struct ift_snode {
    char *elem; 
    struct ift_snode *prev;
    struct ift_snode *next;
} iftSNode;

typedef struct ift_slist {
    long n;
    iftSNode *head;
    iftSNode *tail;
} iftSList;

iftSList *iftCreateSList();
void iftDestroySList(iftSList **SL);
static inline bool iftIsSListEmpty(const iftSList *SL) {
    return (SL->n == 0);
}
void iftInsertSListIntoHead(iftSList *SL, const char *elem);
void iftInsertSListIntoTail(iftSList *SL, const char *elem);
char *iftRemoveSListHead(iftSList *SL);
char *iftRemoveSListTail(iftSList *SL);

// ---------- iftSList.h end
// ---------- iftDialog.h start 

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define MSG_MEMORY_ALLOC_ERROR  "Cannot allocate memory space"
#define MSG_FILE_OPEN_ERROR  "Cannot open file %s"
#define MSG_FILE_FORMAT_ERROR  "Format error in line %d of file %s"
#define IFT_COLOR_YELLOW "\x1B[33m"
#define IFT_COLOR_BLUE "\x1B[34m"
#define IFT_COLOR_RESET "\033[0m"

#ifdef IFT_DEBUG
#define iftAssert(cond, msg, func, ...) iftAssertAux(cond, msg, func, __VA_ARGS__)
#else
#define iftAssert(cond, msg, func, ...) (void)0
#endif

void iftError(const char *msg, const char *func, ...);
void iftWarning(const char *msg, const char *func, ...);
void iftDeprecated(const char *old_function, const char *new_function, const char *message);

// ---------- iftDialog.h end
// ---------- iftStream.h start 

static inline void *iftAlloc(size_t n, size_t size) {
    return calloc(n, size);
}
static inline void iftFree(void *data) {
    if (data != NULL)
        free(data);
}
static inline void *iftRealloc(void *data, size_t size) {
    return realloc(data, size);
}
static inline void iftCopyData(void *dst, const void *src, size_t n, size_t chunk_size) {
    size_t nbytes = n * chunk_size;
    memmove(dst, src, nbytes);
}
char *iftGetLine(FILE *stream);

// ---------- iftStream.h end
// ---------- iftString.h start 

void iftRightTrim(char* s, char c);
iftSList *iftSplitString(const char *phrase, const char *delimiter);
char *iftLowerString(const char *str);
bool iftCompareStrings(const char *str1, const char *str2);
char *iftSplitStringAt(const char *phrase, const char *delimiter, long position);
char *iftCopyString(const char *str, ...);
char *iftRemoveSuffix(const char *str, const char *suffix);
bool iftEndsWith(const char *str, const char *suffix);
bool iftStartsWith(const char *str, const char *prefix);
char *iftConcatStrings(int n, ...);
char *iftRemovePrefix(const char *str, const char *prefix);

// ---------- iftString.h end
// ---------- iftIGraph.h start 

#define   COMPLETE 0 /* graph where all nodes are adjacent to each other */  
#define   EXPLICIT 1 /* graph with adjacency list of the nodes */
#define   IMPLICIT 2 /* graph with translation-invariant adjacency relation */

typedef struct ift_inode {
  int voxel; 
  float weight;
  iftSet *adj;
} iftINode;

typedef struct ift_igraph {
  iftINode *node;
  int nnodes;
  int nfeats;
  iftImage *index;
  float **feat;
  int *label, *marker, *root, *pred;
  double *pvalue;
  iftAdjRel *A;
  char type;
} iftIGraph;

iftIGraph *iftMImageToIGraph(const iftMImage *img, const iftImage *mask);
iftIGraph *iftImplicitIGraph(iftMImage *img, const iftImage *mask, iftAdjRel *A);
void iftIGraphSetWeightForRegionSmoothing(iftIGraph *igraph, const iftImage *img);
void iftIGraphSmoothRegions(iftIGraph *igraph, int num_smooth_iterations);
void iftDestroyIGraph(iftIGraph **igraph);
iftImage *iftIGraphLabel(iftIGraph *igraph);
void iftIGraphSetFWeight(iftIGraph *igraph, iftFImage *weight);
iftFImage *iftIGraphWeight(iftIGraph *igraph);

// ---------- iftIGraph.h end
// ---------- iftDataSet.h start

typedef float (*iftArcWeightFun)(float *f1, float *f2, float *alpha, int n);
float iftDistance1(float *f1, float *f2, float *alpha, int n);
float iftDistance10(float *f1, float *f2, float *alpha, int n);

// ---------- iftDataSet.h end
// ---------- iftSeeds.h start

iftIntArray * iftGridSamplingOnMask(const iftImage *bin_mask, float radius,
                                   int initial_obj_voxel_idx, long n_samples);
float iftEstimateGridOnMaskSamplingRadius(const iftImage *binMask, int initialObjVoxelIdx, int nSamples);
iftSet * iftObjectBorderSet(const iftImage *label_img, iftAdjRel *Ain);
iftImage * iftSelectKLargestRegionsAndPropagateTheirLabels(iftImage *label, iftAdjRel *A, int K);

// ---------- iftSeeds.h end
// ---------- iftDir.h start 

typedef struct ift_dir {
    char *path;
    long nfiles;
    iftFile **files;
    long nsubdirs;
    struct ift_dir **subdirs;
} iftDir;

bool iftDirExists(const char *pathname, ...);
char *iftParentDir(const char *pathname);
void iftMakeDir(const char *dir_path);
iftDir *iftLoadFilesFromDirByRegex(const char *dir_pathname, const char *regex);
iftDir *iftLoadDir(const char *dir_pathname, long hier_levels);
void iftDestroyDir(iftDir **dir);

// ---------- iftDir.h end
// ---------- iftCompression.h start 

struct iftGZipPtr {
    bool withz;
    FILE* nzfptr;
#ifdef HAVE_ZLIB
    gzFile zfptr;
#endif
};

typedef struct iftGZipPtr *iftGZipFile;

iftGZipFile iftGZipOpen(const char *path, const char *mode, bool use_compression);
char *iftGZipGets(char *str, int size, iftGZipFile file);
int iftGZipClose(iftGZipFile *file);
size_t iftGZipRead(void *buf, size_t size, size_t nmemb, iftGZipFile file);
int iftGZipPuts(const char *str, iftGZipFile file);
size_t iftGZipWrite(const void *buf, size_t size, size_t nmemb, iftGZipFile file);

// ---------- iftCompression.h end
// ---------- iftRegex.h start

bool iftRegexMatch(const char *str, const char *regex_pattern, ...);

// ---------- iftRegex.h end 
// ---------- iftImageMath.h start

iftImage *iftAddValue(const iftImage *img, int val);

// ---------- iftImageMath.h end 
// ---------- iftNumerical.h start

iftIntArray *iftIntRange(int begin, int end, int inc);

// ---------- iftNumerical.h end
// ---------- iftSort.h start

#define IFT_MAXVALUE 1

void iftFQuickSort( float *value, int *index, int i0, int i1, uchar order); 
void iftBucketSort(int *value, int *index, int nelems, uchar order);

// ---------- iftSort.h end
// ---------- iftBMap.h start

typedef struct ift_bitmap {
    char *val;
    int nbytes;
    int n;
} iftBMap;

iftBMap *iftCreateBMap(int n);
void iftDestroyBMap(iftBMap **bmap);
static inline void iftBMapSet1(iftBMap *bmap, int b) {
    bmap->val[b >> 3] |= (1 << (b & 0x07));
}
static inline bool iftBMapValue(const iftBMap *bmap, int b) {
    return ((bmap->val[b >> 3] & (1 << (b & 0x07))) != 0);
}

// ---------- iftBMap.h end
// ---------- iftList.h start 

typedef struct _ift_node {
    int elem;
    struct _ift_node *previous;
    struct _ift_node *next;
} iftNode;

typedef struct _ift_list {
    /** Number of Nodes of the Integer Doubly Linked List */
    int n;
    /** Pointer to the Head (first) Node */
    iftNode *head;
    /** Pointer to the Tail (last) Node */
    iftNode *tail;
} iftList;

iftList *iftCreateList();
void iftDestroyList(iftList **L);
bool iftIsEmptyList(const iftList *L);
void iftInsertListIntoTail(iftList *L, int elem);
int iftRemoveListTail(iftList *L);
iftIntArray *iftListToIntArray(const iftList *L);

// ---------- iftList.h end
// ---------- iftRegion.h start 

iftImage* iftRelabelRegions(iftImage* labelled, iftAdjRel* adj_rel);

// ---------- iftRegion.h end
// ---------- iftDict.h start

#define IFT_INITIAL_DICT_SIZE 103
#define IFT_HASH_FAIL -1 // "flag" to indicate that there is no available buckets in the dict
#define iftDictContainKey(KEY, DICT, KEY_INDEX) iftDictContainGKey(iftCreateGVal((KEY)), (DICT), (KEY_INDEX))
#define iftGetBoolValFromDict(KEY,DICT)     iftGetBoolVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_BOOL_TYPE))
#define iftGetCharValFromDict(KEY,DICT)     iftGetCharVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_CHAR_TYPE))
#define iftGetUCharValFromDict(KEY,DICT)    iftGetUCharVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_UCHAR_TYPE))
#define iftGetStrValFromDict(KEY,DICT)      iftGetStrVal(iftCopyGVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_STR_TYPE)))
#define iftGetConstStrValFromDict(KEY,DICT) iftGetConstStrVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_STR_TYPE))
#define iftGetLongValFromDict(KEY,DICT)     iftGetLongVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_LONG_TYPE))
#define iftGetULongValFromDict(KEY,DICT)    iftGetULongVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_ULONG_TYPE))
#define iftGetDblValFromDict(KEY,DICT)      iftGetDblVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_DBL_TYPE))
#define iftGetPtrValFromDict(KEY,DICT)      iftGetPtrVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_PTR_TYPE))
#define iftGetIntArrayFromDict(KEY, DICT)   iftGetIntArrayVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_INT_ARRAY_TYPE))
#define iftGetDblArrayFromDict(KEY, DICT)   iftGetDblArrayVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_DBL_ARRAY_TYPE))
#define iftGetStrArrayFromDict(KEY, DICT)   iftGetStrArrayVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_STR_ARRAY_TYPE))
#define iftGetIntMatrixFromDict(KEY, DICT)  iftGetIntMatrixVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_INT_MATRIX_TYPE))
#define iftGetDblMatrixFromDict(KEY, DICT)  iftGetDblMatrixVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_DBL_MATRIX_TYPE))
#define iftGetStrMatrixFromDict(KEY, DICT)  iftGetStrMatrixVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_STR_MATRIX_TYPE))
#define iftGetDictFromDict(KEY, DICT)       iftGetDictVal(iftGetGValFromDict(iftCreateGVal((KEY)), (DICT), IFT_DICT_TYPE))
#define iftInsertIntoDict(KEY,VAL,DICT) iftInsertKeyValIntoDict(iftCreateGVal(KEY), iftCreateGVal(VAL), DICT)

typedef struct ift_key_val {
    iftGVal key;
    iftGVal val;
    int prev;
    int next;
} iftKeyVal;

struct ift_dict {
    size_t capacity;
    size_t n_busy_buckets;
    iftKeyVal *table;
    bool erase_elements;
    int firstKey;
    int lastKey;
};

bool iftDictContainGKey(iftGVal key, const iftDict *dict, size_t *key_index);
iftGVal iftGetGValFromDict(iftGVal key, const iftDict *dict, iftCDataType true_val_type);
iftDict *iftCreateDict();
iftDict *iftCreateDictWithInitialSize(size_t n);
bool iftIsDictEmpty(const iftDict *dict);
iftDict *iftCopyDict(const iftDict *dict);
iftKeyVal iftDictFirst(const iftDict* dict);
iftKeyVal iftDictNext(const iftDict* dict, const iftKeyVal keyval);
bool iftDictValidKey(const iftKeyVal keyval);
bool iftInsertKeyValIntoDict(iftGVal key, iftGVal val, iftDict *dict);
bool iftIsDictFull(const iftDict *dict);
iftDict *iftCreateDictWithApproxSize(size_t n);
void iftDestroyDict(iftDict **dict);

// ---------- iftDict.h end 
// ---------- iftCommandLineParser.h start

typedef struct ift_cmd_line_opt {
    char short_name[4];
    char long_name[128];
    bool has_arg;
    iftCDataType arg_type;
    bool required;
    char help[2048];
} iftCmdLineOpt;

typedef struct ift_cmd_line_parser {
    /** Program name */
    char program_name[128];
    /** Number of options. */
    int n_opts;
    /** An array with the command line options. */
    iftCmdLineOpt *opts;
    /** Description of the program usage. It is used in the iftPrintUsage(). */
    char *description;
} iftCmdLineParser;

iftCmdLineParser *iftCreateCmdLineParser(const char *description, int n_opts, iftCmdLineOpt cmd_line_opts[]);
iftDict *iftParseCmdLine(int argc, const char *argv[], iftCmdLineParser *parser);
void iftDestroyCmdLineParser(iftCmdLineParser **parser);
void iftPrintUsage(const iftCmdLineParser *parser);
bool iftInsertKeyValIntoDict(iftGVal key, iftGVal val, iftDict *dict);

// ---------- iftCommandLineParser.h end
#endif // _IFT_H_
