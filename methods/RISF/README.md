Image segmentation using dense and sparse hierarchies of superpixels
(c) 2020 Felipe L. Galvao, Silvio J. F. Guimaraes, and Alexandre X. Falcao
------------------------------------------------------------------------------

This software contains superpixel segmentation programs corresponding to the
paper "Image segmentation using dense and sparse hierarchies of superpixels"
submitted to Elsevier's Pattern Recognition (PR).

This software includes:
  - A program to run Recursive Iterative Spanning Forest methods;
  - A program to run the Mid-level region merging method;
  - A program to compute the image segmentation measures used in the paper;
  - Some general superpixel related utility programs.

Please cite the corresponding paper if you use any of this software.

For any inquiries about this software contact felipelemes@outlook.com

------------------------------------------------------------------------------
 1 - How to build
------------------------------------------------------------------------------

This demo is self-contained and includes the libjpeg, libpng, zlib, and a
  subset of the IFT library.

The only requirement is the gcc compiler (tested with version 5.4.0).

This demo was only tested on Linux Ubuntu 18.04, but it should also work
  on other Linux distributions, OSX, and Windows. On Windows you will
  need an environment that supports running a Makefile; for that we
  recommend using the Windows Subsystem for Linux (WSL).

Here are the steps to actually build the demo:

1) Extract "RISFSuperpixels" to some folder, e.g. 
     /home/alice/dev/RISFSuperpixels

2) Move to the root folder 'RISFSuperpixels'

3) Execute:
     make

4) On a successful execution, the folder 'RISFSuperpixels/bin' should
   contain the following programs:
     iftComputeSuperpixelMetrics
     iftForceLabelMapConnectivity
     iftMidLevelRegionMerge
     iftOverlayBorders
     iftRISF_segmentation
     iftShowConnectivityProblems

------------------------------------------------------------------------------
 2 - Usage
------------------------------------------------------------------------------

In all programs, what we refer to as LABEL MAP is an image where pixels
  sharing an intensity value represent a single REGION (superpixel).

Note that 3D images are supported using the SCN format (.scn), which is very
  similar to the Netpbm P5 (.pgm) format.

For the examples provided in each program, it is assumed that the user is
  at the root folder 'RISFSuperpixels'.

Next, we describe what each of the compiled programs does and how to use them.

------------------------------------------------------------------------------
 2.1 - iftRISF_segmentation
------------------------------------------------------------------------------

This program runs specific Recursive Iterative Spanning Forest (RISF) methods
  according to the supplied parameters. More specifically, it performs a
  single recursive superpixel segmentation step and should be executed
  multiple times to build a hierarchy.

Note that this implementation is geared towards the recursive operation over
  Region Adjancecy Graphs (RAGs) and requires an existing segmentation to work.
  For the optimized pixel-level operation, we use the Iterative Spanning 
  Forest (ISF) implementation available at:
  'https://www.ic.unicamp.br/~afalcao/downloads.html'.

Usage: 

    iftRISF_segmentation [OPTIONS]

The following options are required:

    -i, --input-image     [HAS ARG] Input image path.
    -o, --output-image    [HAS ARG] Output label map image.
    -n, --superpixel-num  [HAS ARG] Target number of superpixels.
    -l, --prev-seg        [HAS ARG] Previous segmentation label map.

The following options are optional and customize how the method operates: 
  
    -s, --sampling        [HAS ARG] Seed sampling strategy:
        1 = Grid (DEFAULT) - Only strategy used in the corresponding paper
        2 = Mixed - The entropy based strategy from (Vargas et al., 2019)
        3 = Random
        4 = Geodesic - The maximum geodesic cost strategy from (Galvao et al., 2018)
      
    -c, --path-cost       [HAS ARG] Path-cost function for non-trivial paths:
        1 = Root (DEFAULT)
        2 = Mean
      
    -r, --seed-recomp     [HAS ARG] Seed recomputation method:
        1 = Medoid - Node with feature vector closes to mean
        2 = Centroid
        DEFAULT is based on path cost, medoid for root and centroid for mean.
      
    -a, --alpha           [HAS ARG] Linear weight for the parametric distance.
        Smaller is more compact. Recommended range [0.12,0.5] for natural
        images. DEFAULT is 0.2.
      
    -b, --beta            [HAS ARG] Exponential weight for the parametric distance.
        DEFAULT is 12.
      
    -t, --iters-num       [HAS ARG] Number of iterations. 
        DEFAULT is 10.
      
    -p, --print-opt       [HAS ARG] Printing options:
        0 = Nothing is printed
        1 = Print human-readable stats (DEFAULT)
        2 = Print CSV formatted stats
      
    -x, --output-format   [HAS ARG] Output label map format:
        0 = Deduce by file extension (DEFAULT)
        1 = Netpbm PGM P2 (ASCII)
        2 = Netpbm PGM P5 (binary)

Help Options

    -h, --help            Show help options

Example of how to use this program:

    ./bin/iftRISF_segmentation -i dat/berkley_horse.jpg -o dat/horse_risf.pgm -n 200 -l dat/horse_isf.pgm

    This will compute a superpixel segmentation 'dat/horse_risf.pgm' 
      over the image 'dat/berkley_horse.jpg' with approximately 
      200 regions. Following from the default parameters, this defines
      a RISF method with grid sampling, root path-cost function (alpha = 0.2
      and beta = 12), seed recomputation selecting node closest to mean,
      and 10 iterations.
      
    See iftOverlayBorders to generate a visualization of the label map over
      the original image.

------------------------------------------------------------------------------
 2.2 - iftMidLevelRegionMerge
------------------------------------------------------------------------------

This program runs the proposed mid-level region merge algorithm. Instead of
  storing the entire hierarchy, this program only outputs the specified scale
  (in number of regions). To store multiple observation scales, execute the
  program recursively while lowering the target number of regions.

Note that its interface is analogous to the iftRISF_segmentation program,
  minus the parameters defining an RISF method. The mid-level region merge
  algorithm in itself is parameterless.

Usage:

    iftMidLevelRegionMerge [OPTIONS]

The following options are required:

    -i, --input-image     [HAS ARG] Input image path.
    -o, --output-image    [HAS ARG] Output label map image.
    -n, --superpixel-num  [HAS ARG] Target number of superpixels.
    -l, --prev-seg        [HAS ARG] Previous segmentation label map.

The following options are optional and customize how the method operates:

    -p, --print-opt       [HAS ARG] Printing options:
        0 = Nothing is printed
        1 = Print human-readable stats (DEFAULT)
        2 = Print CSV formatted stats

    -x, --output-format   [HAS ARG] Output label map format:
        0 = Deduce by file extension (DEFAULT)
        1 = Netpbm PGM P2 (ASCII)
        2 = Netpbm PGM P5 (binary)

Help Options

    -h, --help            Show help options

Example of how to use this program:

    ./bin/iftMidLevelRegionMerge -i dat/berkeley_sphinx.jpg -o dat/sphinx_merge.pgm -n 100 -l dat/sphinx_risf.pgm
    
    This will compute a superpixel segmentation 'dat/sphinx_merge.pgm' 
      over the image 'dat/berkley_sphinx.jpg' with exactly 100 regions.

------------------------------------------------------------------------------
 2.3 - iftComputeSegmentationMetrics
------------------------------------------------------------------------------

This program computes the segmentation metrics for the superpixel and large
  segment regimes presented in our work. Given an image segmentation and the
  corresponding set of ground truth segmentations, it outputs a selection of
  metrics according to the supplied options.

The segmentations are interpreted as label maps by default, but binary border
  maps can also be used (see options -ib and -gb). 

Usage:

    iftComputeSegmentationMetrics [OPTIONS]

The following options are required:

    -i, --input-segmentation  [HAS ARG] Input image segmentation path.
    -g, --input-groundtruth   [HAS ARG] Input ground-truth prefix.
        All label maps matching the supplied prefix will be considered.

The following options are optional and customize how the method operates:

    -ib, --is-border-seg      Interpret segmentation as a border map.
    -gb, --is-border-gt       Interpret ground-truth as a border map.
    -t, --boundary-tolerance  [HAS ARG] Maximum distance to count as valid border pixel (default = 2.0).
    -br, --boundary-recall    Compute boundary recall (BR).
    -bp, --boundary-precision Compute boundary precision (BP).
    -bf, --boundary-fscore    Compute boundary f-score (BF).
    -u, --border-union        Compute border statistics for union of all gt images (UBR, UBP, UBF).
    -ue, --undersegmentation  Compute undersegmentation error (UE).
    -as, --achievable-seg     Compute achievable segmentation accuracy (ASA).
    -c, --compactness         Compute compactness.
    -sp, --std-superpixels    Compute standard superpixel metrics (equivalent to -br -ue -as -c).
    -d, --dice                [HAS ARG] Compute dice for [param] objects.
    -to, --topology           Compute topology metric.
    -cv, --covering           Compute covering.
    -vi, --variation-of-info  Compute variation of information.
    -ri, --prob-rand-index    Compute probabilistic rand index.
    -ls, --std-large-segment  Compute standard large segment metrics (equivalent to -cv -vi -ri).
    -p, --print-opt           [HAS ARG] Printing options:
        1 = Print human-readable stats (DEFAULT)
        2 = Print CSV formatted stats
        3 = Print CSV header (dummy values for -i and -g can be used)
        4 = Print CSV header with standard deviation (std)

Help Options

    -h, --help            Show help options

Example of how to use this program:

    ./bin/iftComputeSegmentationMetrics -i dat/sphinx_risf.pgm -g dat/sphinx_gt -sp -gb
    
    This will compute and print the boundary recall, undersegmentation error,
      achievable segmentation accuracy and compactness of the segmentation
      'dat/sphinx_risf' against the border map ground truth segmentations
      'dat/sphinx_gt_1.pgm' and 'dat/sphinx_gt_2.pgm'. 

------------------------------------------------------------------------------
 2.4 - iftOverlayBorders
------------------------------------------------------------------------------

This program creates a visualization of the original image overlaid by the
  borders from one or more segmentations (either label or border map).

Usage:

    iftOverlayBorders <original image> <result path (with extension)> <segmentation 1> [<segmentation 2>...]

Example:

    ./bin/iftOverlayBorders dat/berkeley_horse.jpg dat/vis_horse.png dat/horse_isf.pgm
    
      This will create the visualization 'dat/vis_horse.png' with the borders
        from the label map 'dat/horse_isf.pgm' over the image 'dat/berkeley_horse'.

------------------------------------------------------------------------------
 2.5 - iftShowConnectivityProblems
------------------------------------------------------------------------------

This program examines a given label map to check if it is 4-connected,
  8-connected or not connected. It optionally creates a visualization
  of identified connectivity problems.

Usage:

    iftShowConnectivityProblems <input label map> [<original image> <output visualization prefix>]

Example:

    ./bin/iftShowConnectivityProblems dat/horse_isf.pgm dat/berkeley_horse.jpg dat/vis_connectivity
    
      This will tell if the label map 'dat/horse_isf.pgm' is connected and
        what type of connectivity it is (4-connected in this case). If any
        disconnected regions were found, each one would create a visualization
        with prefix 'dat/vis_connectivity' overlaying the disjoint region borders
        over the original image 'dat/berkeley_horse'.

------------------------------------------------------------------------------
 2.6 - iftForceLabelMapConnectivity
------------------------------------------------------------------------------

This program uses a connected components algorithm to relabel the regions
  of a label map so that disjoint regions have different labels. It also
  optionally filters all regions below a given size in pixels, merging
  them with an arbitrary neighboring region.

Usage:

    iftForceLabelMapConnectivity <input label map> <output label map> [minimum region size to avoid filter]

Example:

    ./bin/iftForceLabelMapConnectivity dat/horse_isf.pgm dat/horse_isf_relabel.pgm 50
  
    This will relabel each connected component of the label map 
      'dat/horse_isf.pgm', filtering those with size below 50 pixels
      and saving the result on 'dat/horse_isf_relabel.pgm'.


