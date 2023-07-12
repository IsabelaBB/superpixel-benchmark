## External Libraries

This file describes the necessary steps to enable distinct features within the library. For convenience, copy of **some** the static libraries used in this project are provided within this folder. Still, for all external libraries, their links are provided for more details on how to download and use them in your own system.

### Flags

In this project, the following compilation flags are provided:
- **IFT_DEBUG**: This flags enables the library's debug mode. It enables specific GNU GCC compiling flags for code inspection (_e.g._ unused variables, pointer arithmetic, etc.), enables printing debug and warning messages, and enables parameter checking for most functions.
- **IFT_LIBJPEG**: Enables support to **Joint Photographic Experts Group** (JPEG) files through the [LibJPEG](http://libjpeg.sourceforge.net/) library. In this project, the version 9.0 was used.
- **IFT_LIBPNG**: Enables support to **Portable Network Graphics** (PNG) files through the [LibPNG](https://libpng.sourceforge.io/index.html) library. In this project, the version 1.6.29 was used.
- **IFT_OMP**: Enables multithreading with [Open MultiProcessing](https://www.openmp.org/) (OpenMP). In this project, the version 4.5 was used.

In order to enable a supporting feature, one may simply define the respective variable in the Makefile as shown below.
```bash
IFT_DEBUG = YES
IFT_LIBJPEG = YES
IFT_LIBPNG = YES
IFT_OMP = YES
```
  Note that such value must be an **exact match** and, therefore, any mismatch will be considered as a negative option for this feature. Pay attention to whitespaces!

### Installation
  
  For convenience, in this folder, a copy of the libraries used in this project are provided. One may simply extract them within this folder and compile properly running the command below in their respective folders.
```bash
make
```
  Such steps should be enough since the Makefile expects those libraries to be located in this folder. However, if one desire, it is possible to indicate another version of one of those libraries (whether it is a static or shared one) by modifying the respective variables within the Makefile to the desired location.
```bash
INCS += -I$(HOME_DIR)/externals/libjpeg/include
LIBS_LD += -L$(HOME_DIR)/externals/libjpeg/lib
#(...)
INCS += -I$(HOME_DIR)/externals/libpng/include
LIBS_LD += -L$(HOME_DIR)/externals/libpng/lib
# (...)
INCS += -I$(HOME_DIR)/externals/libpng/zlib/include
LIBS_LD += -L$(HOME_DIR)/externals/libpng/zlib/lib
```
  It is **important** to notice that LibPNG library requires the installation of the [ZLib](https://zlib.net/) library for compressing and decompressing files. In this project, version 1.2.11 was used and it is also provided within the **libpng.zip** file. For most **Linux-based** systems, the following apt-get command should be enough for installing all the external libraries used in this project.
```bash
sudo apt-get install libjpeg9-dev libpng-dev zlib1g-dev build-essential
```