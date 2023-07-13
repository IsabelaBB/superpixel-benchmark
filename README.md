# Superpixel Segmentation Review

## Intrduction
Superpixel segmentation consists of partitioning images into regions composed of similar and connected
pixels. Its methods have been widely used in many computer vision applications since it allows for reducing
the workload, removing redundant information, and preserving regions with meaningful features. Due to
the rapid progress in this area, the literature fails to catch up on more recent works among the compared
ones and to categorize the methods according to all existing strategies. This repository allows to evaluate 20 strategies based on nine criteria: connectivity, compactness, delineation, control over the number of superpixels, color homogeneity, robustness, running time, stability, and visual quality. In [1] we provide a comprehensive review with new taxonomy for superpixel segmentation and revisit the recent and popular literature according to our taxonomy. Our experiments show the trends of each approach in pixel clustering and discuss individual trade-offs. 

## Superpixel methods
We evaluated 20 superpixel methods and provided the code for each in this repository in [methods](methods/). Be aware that we have adapted these codes to match our system setup, and the benchmark inputs and outputs. The following Table presents the reference paper and code of each method. 

| Method   | Ref. | Code link |
|----------|------|-----------|
| CRS      |[Paper](https://doi.org/10.1007/978-3-642-40395-8_21)|[CRS](https://github.com/davidstutz/superpixel-benchmark)|[CRS](https://github.com/davidstutz/superpixel-benchmark)
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
| SLIC     |[Paper](https://doi.org/10.1109/TPAMI.2012.120)|[SLIC](https://www.epfl.ch/labs/ivrl/research/slic-superpixels/)|
| SNIC     |[Paper](https://openaccess.thecvf.com/content_cvpr_2017/html/Achanta_Superpixels_and_Polygons_CVPR_2017_paper.html)|[SNIC](https://github.com/achanta/SNIC)|

## Datasets
We selected four datasets which impose different challenges for superpixel segmentation:
[Birds](https://doi.org/10.1109/SIBGRAPI.2016.047); [Insects](https://doi.org/10.1109/SIBGRAPI.2016.047); [Sky](https://doi.org/10.1109/SIBGRAPI.2015.20); and [ECSSD](https://doi.org/10.1109/TPAMI.2015.2465960). Birds has 150 natural images of birds with thin and elongated objects’ parts. Similarly, Insects has 130 images of invertebrates with less
texture on background regions. Sky has 60 images for sky segmentation with large homogeneous
regions with subtle luminosity variations. Finally, the Extended Complex Scene Saliency Dataset
(ECSSD) is composed of 1000 images with objects and backgrounds whose textures are
complex. 

## Evaluation measures
We assess connectivity, compactness, delineation, color homogeneity, robustness, running time, stability, control over the number of superpixels, and visual quality in superpixel segmentation methods. The results provide valuable insights into the pros and cons of the methods, supporting the choice of the most suitable one for a given application. This code contais five superpixel evaluation measures: Similarity between Image and Reconstruction from Superpixels (SIRS)[1], Boundary Recall (BR)[2], Explained Variation (EV)[3], Undersegmentation Error (UE)[4], and Compactness (CO)[5] ([evaluation](evaluation/README.md)). In addition, we provide code to assess running time, control over the number of superpixels, connectivity, and robustness ([evaluation](evaluation/README.md)). 


