Iterative Spanning Forest - Superpixel/Supervoxel Segmentation Demo

(c) 2017 John E. Vargas-Muñoz, Ananda S. Chowdhury, Eduardo B. Alexandre, Felipe L. Galvão, Paulo A. Vechiatto Miranda, and Alexandre X. Falcão
------------------------------------------------------------------------------------------------------------------------------

This software contains superpixel segmentation programs corresponding to the paper 
"An Iterative Spanning Forest Framework for Superpixel Segmentation" submitted to IEEE TIP. 
The software consists of the GFT library and programs to generate 2D and 3D superpixel segmentations.

Please cite the version in https://arxiv.org/abs/1801.10041

HOW TO BUILD:

1) Extract "isf_demo" files to some folder, e.g. /home/pmiranda/prog/lib/isf_demo

2) Create a environment variable GFT_DIR and add to .bashrc:
export GFT_DIR=/home/pmiranda/prog/lib/isf_demo

3) In "isf_demo" root folder execute:
make clean
make

4) If you get the error message "fatal error: zlib.h: No such file or directory", then install zlib package: zlib1g-dev
Repeat step 3.

HOW TO USE:

After a successful build, there will be two binaries for 2D images and 3D images, respectively:
${GFT_DIR}/bin/superpixels
${GFT_DIR}/bin/supervoxels

The 3D images:

We are using a file format (SCN) very similar to pgm (P5) for
grayscale images. The header is

<type>
<xsize> <ysize> <zsize>
<dx> <dy> <dz>
<depth>
<binary array of voxel values> 

See dat/image2.scn in some text editor.

type is SCN.
xsize, ysize, and zsize are the dimensions of the image.
dx, dy, and dz are the dimensions of the voxels (it should be
    	       isotropic for better results).
depth is the number of bits per voxel (it can be 8 or 16 bits)
binary array of voxel values is written in the increasing order of x,
y, and then z coordinates.

The resulting .scn file contains the supervoxel labels from 1 to k.


Run either binary without arguments to check their usage.

Examples: 

./bin/superpixels 1 dat/image1.png 100 0.12 0 0 out/image1
Will generate the ISF-GRID-MEAN segmentation label map "out/image1.pgm" and a visualization "out/image1.ppm".

./bin/superpixels 1 dat/image1.png 100 0.12 1 0 out/image1
Will generate the ISF-MIX-MEAN segmentation label map "out/image1.pgm" and a visualization "out/image1.ppm".

./bin/superpixels 1 dat/image1.png 100 0.12 0 1 out/image1
Will generate the ISF-GRID-ROOT segmentation label map "out/image1.pgm" and a visualization "out/image1.ppm".

./bin/superpixels 1 dat/image1.png 100 0.12 1 1 out/image1
Will generate the ISF-MIX-ROOT segmentation label map "out/image1.pgm" and a visualization "out/image1.ppm".

Similarly,

./bin/supervoxels dat/image2.scn 1000 0.08 0 0 out/image2 Will generate
the ISF-GRID-MEAN segmentation label map "out/image2.scn".

This is a MR brain image. From the same link you have gotten ISF, you
can get IVS -- a software tool to visualize the 3D rendition of the
label image by typing

ivs ./out/image2.scn

ivs can also render the 3D supervoxels as objects if you load the
original image and then load the label image. However, it is
unfortunately limited to manipulate 10 objects (supervoxels) only. If
you use more, you can display the rendition but not manipulate the
labels correctly.

Another example involves a CT image of the knee. Try 

./bin/supervoxels dat/image3.scn 10 0.5 0 0 out/image3. It will
generate the ISF-GRID-MEAN segmentation label map "out/image3.scn"
where the patella bone corresponds to one of the supervoxels.
