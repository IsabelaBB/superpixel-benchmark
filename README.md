# Superpixel Segmentation Review

## Intrduction
Superpixel segmentation consists of partitioning images into regions composed of similar and connected
pixels. Its methods have been widely used in many computer vision applications since it allows for reducing
the workload, removing redundant information, and preserving regions with meaningful features. Due to
the rapid progress in this area, the literature fails to catch up on more recent works among the compared
ones and to categorize the methods according to all existing strategies. This repository allows to evaluate 20 strategies based on nine criteria: connectivity, compactness, delineation, control over the number of superpixels, color homogeneity, robustness, running time, stability, and visual quality. In [1] we provide a comprehensive review with new taxonomy for superpixel segmentation and revisit the recent and popular literature according to our taxonomy. Our experiments show the trends of each approach in pixel clustering and discuss individual trade-offs. 

## Superpixel methods
We evaluated 23 superpixel methods (and a grid segmentation baseline) and provided the code for each in this repository in [methods](methods/). Be aware that we have adapted these codes to match our system setup, and the benchmark inputs and outputs. Some of these methods require a saliency map (ODISF and SICLE) or a contour prior map (SCALP) to compute superpixels. In [others](others/) folder, we provide the code of [U-2-Net](https://github.com/xuebinqin/U-2-Net) [2] and [Structured Edge Detection Toolbox](https://github.com/pdollar/edges) [3,4,5] to compute a saliency map and a contour prior map, repectively. Also, DAL-HERS require a pre-trained model available [here]([https://drive.google.com/file/d/14-uaeMAihLdMepfZAth19T1pfZIoMcaE/view?usp=sharing](https://github.com/hankuipeng/DAL-HERS/tree/master)). The following Table presents the reference paper and code of each superpixel method. 

| Method   | Ref. | Code link |
|----------|------|-----------|
| AINET      |[Paper](https://openaccess.thecvf.com/content/ICCV2021/html/Wang_AINet_Association_Implantation_for_Superpixel_Segmentation_ICCV_2021_paper.html)|[AINET](https://github.com/davidstutz/superpixel-benchmark)|[AINET]([https://github.com/davidstutz/superpixel-benchmark](https://github.com/YanFangCS/AINET))|
| CRS      |[Paper](https://doi.org/10.1007/978-3-642-40395-8_21)|[CRS](https://github.com/davidstutz/superpixel-benchmark)|[CRS](https://github.com/davidstutz/superpixel-benchmark)|
| DAL-HERS |[Paper](http://openaccess.thecvf.com/content/WACV2022/html/Peng_HERS_Superpixels_Deep_Affinity_Learning_for_Hierarchical_Entropy_Rate_Segmentation_WACV_2022_paper.html)|[DAL-HERS](https://github.com/hankuipeng/DAL-HERS)|
| DISF     |[Paper](https://doi.org/10.1109/LSP.2020.3015433)|[DISF](https://github.com/LIDS-UNICAMP/ODISF)|
| DRW      |[Paper](https://doi.org/10.1109/TIP.2020.2967583)|[DRW](https://github.com/zh460045050/DRW)|
| ERGC     |[Paper](https://doi.org/10.1016/j.irbm.2013.12.007)|[ERGC](https://github.com/davidstutz/superpixel-benchmark)|
| ERS      |[Paper](https://doi.org/10.1109/CVPR.2011.5995323)|[ERS](https://github.com/akanazawa/collective-classification/tree/master/segmentation)|
| ETPS     |[Paper](http://openaccess.thecvf.com/content_cvpr_2015/html/Yao_Real-Time_Coarse-to-Fine_Topologically_2015_CVPR_paper.html)|[ETPS](https://bitbucket.org/mboben/spixel/src/master/)|
| GMMSP    |[Paper](https://doi.org/10.1109/TIP.2018.2836306)|[GMMSP](https://github.com/ahban/GMMSP-superpixel)|
| IBIS     |[Paper](https://doi.org/10.1109/ACCESS.2021.3081919)|[IBIS](https://github.com/xapha/IBIS)|
| ISF      |[Paper](https://doi.org/10.1109/TIP.2019.2897941)|[ISF](https://www.ic.unicamp.br/afalcao/downloads.html)|
| LNSNet   |[Paper](http://openaccess.thecvf.com/content/CVPR2021/html/Zhu_Learning_the_Superpixel_in_a_Non-Iterative_and_Lifelong_Manner_CVPR_2021_paper.html)|[LNSNet](https://github.com/zh460045050/LNSNet)|
| LSC      |[Paper](http://openaccess.thecvf.com/content_cvpr_2015/html/Li_Superpixel_Segmentation_Using_2015_CVPR_paper.html)|[LSC](https://jschenthu.weebly.com/projects.html)|
| ODISF    |[Paper](https://doi.org/10.1109/SIBGRAPI54419.2021.00054)|[OISF](https://github.com/LIDS-UNICAMP/ODISF)|
| RSS      |[Paper](https://doi.org/10.1007/s11263-020-01352-9)|[RSS](https://github.com/dfchai/Rooted-Spanning-Superpixels)|
| SCALP    |[Paper](https://doi.org/10.1016/j.cviu.2018.01.006)|[SCALP](https://github.com/rgiraud/scalp)|
| SEEDS    |[Paper](https://doi.org/10.1007/978-3-642-33786-4_2),[Paper](https://doi.org/10.1007/s11263-014-0744-2)|[SEEDS](https://github.com/davidstutz/superpixel-benchmark)
| SH       |[Paper](https://doi.org/10.1109/TIP.2018.2836300)|[SH](https://github.com/semiquark1/boruvka-superpixel)|
| SICLE    |[Paper](https://doi.org/10.1007/978-3-031-19897-7_21),[Paper](https://doi.org/10.48550/arXiv.2204.03533)|[SICLE](https://github.com/LIDS-UNICAMP/SICLE)|
| SIN      |[Paper](https://doi.org/10.1007/978-3-030-89370-5_22)|[SIN](https://github.com/yuanqqq/SIN)|
| SLIC     |[Paper](https://doi.org/10.1109/TPAMI.2012.120)|[SLIC](https://www.epfl.ch/labs/ivrl/research/slic-superpixels/)|
| SNIC     |[Paper](https://openaccess.thecvf.com/content_cvpr_2017/html/Achanta_Superpixels_and_Polygons_CVPR_2017_paper.html)|[SNIC](https://github.com/achanta/SNIC)|
| SSFCN    |[Paper](https://openaccess.thecvf.com/content_CVPR_2020/html/Yang_Superpixel_Segmentation_With_Fully_Convolutional_Networks_CVPR_2020_paper.html)|[SSFCN](https://github.com/fuy34/superpixel_fcn)|

## Datasets
We selected five datasets which impose different challenges for superpixel segmentation:
[Birds](https://doi.org/10.1109/SIBGRAPI.2016.047) [6]; [Insects](https://doi.org/10.1109/SIBGRAPI.2016.047) [6]; [Sky](https://doi.org/10.1109/SIBGRAPI.2015.20) [7]; and [ECSSD](https://doi.org/10.1109/TPAMI.2015.2465960) [8]; [NYUV2](https://doi.org/10.1007/978-3-642-33715-4_54) [9]. Birds has 150 natural images of birds with thin and elongated objects’ parts. Similarly, Insects has 130 images of invertebrates with less texture on background regions. Sky has 60 images for sky segmentation with large homogeneous regions with subtle luminosity variations. The Extended Complex Scene Saliency Dataset (ECSSD) is composed of 1000 images with objects and backgrounds whose textures are complex. Finally, the NYUV2 dataset is composed of 1449 video sequences from several indoor scenes recorded by Microsoft Kinect. 

## Evaluation measures
We assess connectivity, compactness, delineation, color homogeneity, robustness, running time, stability, control over the number of superpixels, and visual quality in superpixel segmentation methods. This repository contais our evaluation code ([evaluation](evaluation/README.md)) with five superpixel evaluation measures: Similarity between Image and Reconstruction from Superpixels (SIRS) [1], Boundary Recall (BR) [10], Explained Variation (EV) [11], Undersegmentation Error (UE) [12], and Compactness (CO) [13]. In addition, we provide code to assess running time, control over the number of superpixels, connectivity, and robustness. 

## Compiling and Running
- To **compile** all files: `bash make.sh`
  - Every method that needs an executable file contains a _Makefile_ in its folder. `make.sh` just call each one and call the Makefile in evaluation folder.
- [Scripts](Scripts/) has bash and python scripts to run saliency/contour maps ([others](Scripts/others/)), segmentation ([Segmentation](Scripts/Segmentation/)), and evaluation ([Evaluation](Scripts/Evaluation/)).
- To run **segmentation**, execute `bash run_segmentation.sh` ([run_segmentation.sh](Scripts/Segmentation/run_segmentation.sh)).
    - The folder [Methods](Scripts/Segmentation/Methods) has scripts for each method.
    - One may specify the methods, datasets, directories and parameters in `run_segmentation.sh`.
    - Similarly, one can generate segmentation for images with average blur or salt and peper noise. 
- To generate **segmentation** for images with average blur or salt and peper **noise**:
  - Run `python3 add_noise.py` ([add_noise.py](Scripts/Evaluation/add_noise.py)) to generate noised images
  - Then, generate the segmentations by running `run_segmentation_robustness.sh` ([run_segmentation_robustness.sh](Scripts/Segmentation/run_segmentation_robustness.sh)).
  - One may specify the methods, datasets, directories and parameters in `run_segmentation_robustness.sh`.
- To **evaluate**:
  - Run `bash run_eval.sh` ([run_eval.sh](Scripts/Evaluation/run_eval.sh)) to assess BR, UE, SIRS, EV, CO, superpiixel connectivity, and the control over the superpixel number.
  - Run `bash run_eval_robustness.sh` ([run_eval.sh](Scripts/Evaluation/run_eval_robustness.sh)) to assess BR, UE, SIRS, EV, CO, superpiixel connectivity, and the control over the superpixel number in the noised images.
  - Run `python3 eval_stability.py` ([eval_stability.py](Scripts/Evaluation/eval_stability.py)) to assess the stability of superpixel evaluation.


## References
[1] Isabela B Barcelos, Felipe De C Belém, Leonardo De M João, Alexandre X Falcão, and Guimarães Silvio JF. 2022. Improving color homogeneity measure in superpixel segmentation assessment. In 2022 35th SIBGRAPI Conference on Graphics, Patterns and Images (SIBGRAPI), Vol. 1. 79–84. https://doi.org/10.1109/SIBGRAPI55357.2022.9991772.

[2] QIN, Xuebin et al. U2-Net: Going deeper with nested U-structure for salient object detection. Pattern recognition, v. 106, p. 107404, 2020.

[3] DOLLÁR, Piotr; ZITNICK, C. Lawrence. Structured forests for fast edge detection. In: Proceedings of the IEEE international conference on computer vision. 2013. p. 1841-1848.

[4] DOLLÁR, Piotr; ZITNICK, C. Lawrence. Fast edge detection using structured forests. IEEE transactions on pattern analysis and machine intelligence, v. 37, n. 8, p. 1558-1570, 2014.

[5] ZITNICK, C. Lawrence; DOLLÁR, Piotr. Edge boxes: Locating object proposals from edges. In: Computer Vision–ECCV 2014: 13th European Conference, Zurich, Switzerland, September 6-12, 2014, Proceedings, Part V 13. Springer International Publishing, 2014. p. 391-405.

[6] Lucy A. C. Mansilla and Paulo A. V. Miranda. 2016. Oriented Image Foresting Transform Segmentation: Connectivity Constraints with Adjustable Width. In 2016 29th SIBGRAPI Conference on Graphics, Patterns and Images (SIBGRAPI). 289–296. https://doi.org/10.1109/SIBGRAPI.2016.047

[7] Eduardo Barreto Alexandre, Ananda Shankar Chowdhury, Alexandre Xavier Falcao, and Paulo A. Vechiatto Miranda. 2015. IFT-SLIC: A General Framework for Superpixel Generation Based on Simple Linear Iterative Clustering and Image Foresting Transform. In 2015 28th SIBGRAPI Conference on Graphics, Patterns and Images. 337–344. https://doi.org/10.1109/SIBGRAPI.2015.20

[8] Jianping Shi, Qiong Yan, Li Xu, and Jiaya Jia. 2015. Hierarchical image saliency detection on extended CSSD. IEEE transactions on pattern analysis and machine intelligence 38, 4 (2015), 717–729.

[9] SILBERMAN, Nathan et al. Indoor segmentation and support inference from rgbd images. In: Computer Vision–ECCV 2012: 12th European Conference on Computer Vision, Florence, Italy, October 7-13, 2012, Proceedings, Part V 12. Springer Berlin Heidelberg, 2012. p. 746-760.

[10] David R Martin, Charless C Fowlkes, and Jitendra Malik. 2004. Learning to detect natural image boundaries using local brightness, color, and texture cues. IEEE transactions on pattern analysis and machine intelligence 26, 5 (2004), 530–549. https://doi.org/10.1109/TPAMI.2004.1273918

[11] Alastair P Moore, Simon JD Prince, Jonathan Warrell, Umar Mohammed, and Graham Jones. 2008. Superpixel lattices. In 2008 IEEE conference on computer vision and pattern recognition. IEEE, 1–8.

[12] Peer Neubert and Peter Protzel. 2012. Superpixel benchmark and comparison. In Proc. Forum Bildverarbeitung, Vol. 6. KIT Scientific Publishing, 1–12.

[13] Alexander Schick, Mika Fischer, and Rainer Stiefelhagen. 2012. Measuring and evaluating the compactness of superpixels. In Proceedings of the 21st international conference on pattern recognition (ICPR2012). IEEE, 930–934.

