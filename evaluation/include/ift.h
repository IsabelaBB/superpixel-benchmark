/*****************************************************************************\
* ift.h
*
* AUTHOR    : Felipe Belem (Org.) and Alexandre Falcao et. al.
* DATE      : 2021-02-05
* LICENSE   : MIT License
* COPYRIGHT : Alexandre Xavier Falcao 2012-2021
* EMAIL     : felipe.belem@ic.unicamp.br
*             afalcao@ic.unicamp.br
\*****************************************************************************/
#ifndef _IFT_H_
#define _IFT_H_

#ifdef __cplusplus
extern "C" {
#endif


// ---------- iftBasicDataTypes.h start 

#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <regex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define IFT_STR_DEFAULT_SIZE 4096

typedef struct timeval timer;
typedef unsigned char  uchar;
typedef unsigned short ushort;
typedef unsigned int   uint;
typedef unsigned long  ulong;
#ifndef  __cplusplus
typedef long long llong;
typedef unsigned long long ullong;
#endif


typedef struct ift_size {
    int x, y, z;
} iftSize;

// ---------- iftBasicDataTypes.h end
// ---------- iftIntArray.h start

struct ift_int_array {
    long n;
    int *val;
};

// ---------- iftFloatArray.h end
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

#define iftMax(x,y) (((x) > (y))?(x):(y))
#define iftMin(x,y) (((x) < (y))?(x):(y))
#define iftRound(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))
#define iftVoxelDistance(u,v) (sqrtf((u.x-v.x)*(u.x-v.x)+(u.y-v.y)*(u.y-v.y)+(u.z-v.z)*(u.z-v.z)))

void iftSkipComments(FILE *fp);
double iftLog(double val, double base);

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

typedef struct ift_fcolor {
	float val[3];
} iftFColor;

iftColor iftRGBColor(int R, int G, int B);
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

// ---------- iftDHeap.h end
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

// ---------- iftBMap.h end
// ---------- iftImage.h start 

#define iftImgVal(img, x, y, z) img->val[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftImgCb(img, x, y, z) img->Cb[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftImgCr(img, x, y, z) img->Cr[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]

#define iftImgR(img, x, y, z) img->R[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftImgG(img, x, y, z) img->G[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftImgB(img, x, y, z) img->B[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]

#define iftImgAlpha(img, x, y, z) img->alpha[((x)+(img)->tby[(y)]+(img)->tbz[(z)])]
#define iftCopyVoxelSize(src, dst) (dst)->dx = (src)->dx; (dst)->dy = (src)->dy; (dst)->dz = (src)->dz;
#define iftValidVoxel(img, v)((v.x >= 0) && (v.x <= ((img)->xsize - 1)) && (v.y >= 0) && (v.y <= ((img)->ysize - 1)) && (v.z >= 0) && (v.z <= ((img)->zsize - 1)))
#define iftGetXCoord(s, p) (((p) % (((s)->xsize)*((s)->ysize))) % (s)->xsize)
#define iftGetYCoord(s, p) (((p) % (((s)->xsize)*((s)->ysize))) / (s)->xsize)
#define iftGetZCoord(s, p) ((p) / (((s)->xsize)*((s)->ysize)))
#define iftImgVoxelVal(i, v) i->val[iftGetVoxelIndex(i, v)]
#define iftImgVoxelCb(i, v) i->Cb[iftGetVoxelIndex(i, v)]
#define iftImgVoxelCr(i, v) i->Cr[iftGetVoxelIndex(i, v)]
#define iftDiagonalSize(s) (iftRound(sqrtf(s->xsize*s->xsize + s->ysize*s->ysize + s->zsize*s->zsize)))

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
void iftDestroyImage(iftImage **img);
void iftWriteImageByExt(const iftImage *img, const char *filename, ...);
static inline bool iftIsColorImage(const iftImage *img) {
    return ((img->Cb != NULL) && (img->Cr != NULL));
}

iftImage *iftReadImage(const char *filename, ...);
iftImage* iftReadImagePNG(const char* format, ...);
iftImage* iftReadImageJPEG(const char* format, ...);
iftImage *iftReadImageP5(const char *filename, ...);
iftImage *iftReadImageP6(const char *filename, ...);
iftImage *iftReadImageP2(const char *filename, ...);
iftImage  *iftCreateImageFromBuffer(int xsize,int ysize,int zsize, int *val);

void  iftSetAlpha(iftImage *img, ushort value);
void iftSetCbCr(iftImage *img, ushort value);

// ---------- iftMatrix.h end
// ---------- iftMImage.h start 

int *iftAllocIntArray(long n);
ushort *iftAllocUShortArray(long n);
uchar *iftAllocUCharArray(long n);
char *iftAllocCharArray(long n);


#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define MSG_MEMORY_ALLOC_ERROR  "Cannot allocate memory space"
#define MSG_FILE_OPEN_ERROR  "Cannot open file %s"

void iftError(const char *msg, const char *func, ...);
void iftWarning(const char *msg, const char *func, ...);

// ---------- iftDialog.h end
// ---------- iftStream.h start 

static inline void *iftAlloc(size_t n, size_t size) {
    return calloc(n, size);
}
static inline void iftFree(void *data) {
    if (data != NULL)
        free(data);
}

// ---------- iftStream.h end
// ---------- iftString.h start 
char *iftLowerString(const char *str);
bool iftCompareStrings(const char *str1, const char *str2);

// ---------- iftString.h end
// ---------- iftDir.h start 

typedef struct ift_dir {
    char *path;
    long nfiles;
    iftFile **files;
    long nsubdirs;
    struct ift_dir **subdirs;
} iftDir;

bool iftDirExists(const char *pathname, ...);

// ---------- iftDir.h end
// ---------- iftRegex.h start

bool iftRegexMatch(const char *str, const char *regex_pattern, ...);

// ---------- iftRegex.h end 
// ---------- iftNumerical.h start


#ifdef __cplusplus
}
#endif

#endif // _IFT_H_
