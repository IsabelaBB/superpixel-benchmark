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

// ---------- iftBMap.c end
// ---------- iftImage.c start
#if IFT_LIBPNG
#include <png.h>
#endif
#if IFT_LIBJPEG
#include <jpeglib.h>
#endif

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

bool iftFileExists(const char *pathname) 
{
    return (iftPathnameExists(pathname) && !iftDirExists(pathname));
}

char *iftAllocCharArray(long n) 
{
    char *v = NULL;

    v = (char *) iftAlloc(n, sizeof(char));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocCharArray");
    return (v);
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

bool iftCompareStrings(const char *str1, const char *str2) 
{
    if (str1 == NULL)
        iftError("First String is NULL", "iftCompareStrings");
    if (str2 == NULL)
        iftError("Second String is NULL", "iftCompareStrings");

    return (strcmp(str1, str2) == 0);
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

ushort *iftAllocUShortArray(long n) 
{
    ushort *v = NULL;

    v = (ushort *) iftAlloc(n, sizeof(ushort));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocUShortArray");
    return(v);
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

void iftSkipComments(FILE *fp)
{
    //skip for comments
    while (fgetc(fp) == '#') {
        while (fgetc(fp) != '\n');
    }
    fseek(fp,-1,SEEK_CUR);
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

double iftLog(double val, double base) 
{
    return (log(val) / log(base));
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

uchar *iftAllocUCharArray(long n) 
{
    uchar *v = NULL;

    v = (uchar *) iftAlloc(n, sizeof(uchar));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocUCharArray");
    return (v);
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

void iftWarning(const char *msg, const char *func, ...) 
{
    va_list args;
    char    final_msg[4096];
    
    va_start(args, func);
    vsprintf(final_msg, msg, args);
    va_end(args);
    
    fprintf(stdout, "\nWarning in %s: \n%s\n", func, final_msg);
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
        case JCS_EXT_RGB:
    
            iftError("Big gamut red/green/blue color space not supported", "iftReadImageJPEG");

            break;
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

/////////////////////////////////////////////////////////////////

