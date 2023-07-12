# Superpixel Segmentation Review

## Intrduction

Superpixel segmentation consists of partitioning images into regions composed of similar and connected
pixels. Its methods have been widely used in many computer vision applications since it allows for reducing
the workload, removing redundant information, and preserving regions with meaningful features. Due to
the rapid progress in this area, the literature fails to catch up on more recent works among the compared
ones and to categorize the methods according to all existing strategies. This repository allows to evaluate 20 strategies based on nine criteria: connectivity, compactness, delineation, control over the number of superpixels, color homogeneity, robustness, running time, stability, and visual quality. In [1] we provide a comprehensive review with new taxonomy for superpixel segmentation and revisit the recent and popular literature according to our taxonomy. Our experiments show the trends of each approach in pixel clustering and discuss individual trade-offs. 

## Superpixel methods

| Method   | Ref. | Code link |
|----------|------|-----------|
| CRS      |[Paper]()|[CRS](https://github.com/davidstutz/superpixel-benchmark)|
| DAL-HERS |[Paper](http://openaccess.thecvf.com/content/WACV2022/html/Peng_HERS_Superpixels_Deep_Affinity_Learning_for_Hierarchical_Entropy_Rate_Segmentation_WACV_2022_paper.html)|[DAL-HERS](https://github.com/hankuipeng/DAL-HERS)|
| DSIF     |[Paper]()|[DISF](https://github.com/LIDS-UNICAMP/ODISF)|
| DRW      |[Paper]()|[DRW](https://github.com/zh460045050/DRW)|
| ERGC     |[Paper]()|[ERGC](https://github.com/davidstutz/superpixel-benchmark)|
| ERS      |[Paper](https://ieeexplore.ieee.org/abstract/document/5995323/)|[ERS](https://github.com/akanazawa/collective-classification/tree/master/segmentation)|
| ETPS     |[Paper](http://openaccess.thecvf.com/content_cvpr_2015/html/Yao_Real-Time_Coarse-to-Fine_Topologically_2015_CVPR_paper.html)|[ETPS](https://bitbucket.org/mboben/spixel/src/master/)|
| GMMSP    |[Paper]()|[GMMSP](https://github.com/ahban/GMMSP-superpixel)|
| IBIS     |[Paper]()|[IBIS](https://github.com/xapha/IBIS)|
| ISF      |[Paper](https://ieeexplore.ieee.org/abstract/document/8636546/)|[ISF](https://www.ic.unicamp.br/afalcao/downloads.html)|
| LNSNet   |[Paper](http://openaccess.thecvf.com/content/CVPR2021/html/Zhu_Learning_the_Superpixel_in_a_Non-Iterative_and_Lifelong_Manner_CVPR_2021_paper.html)|[LNSNet](https://github.com/zh460045050/LNSNet)|
| LSC      |[Paper](http://openaccess.thecvf.com/content_cvpr_2015/html/Li_Superpixel_Segmentation_Using_2015_CVPR_paper.html)|[LSC](https://jschenthu.weebly.com/projects.html)|
| ODISF    |[Paper]()|[OISF](https://github.com/LIDS-UNICAMP/ODISF)|
| RSS      |[Paper]()|[RSS](https://github.com/dfchai/Rooted-Spanning-Superpixels)|
| SCALP    |[Paper]()|[SCALP](https://github.com/rgiraud/scalp)|
| SEEDS    |[Paper](https://link.springer.com/chapter/10.1007/978-3-642-33786-4_2)|[SEEDS](https://github.com/davidstutz/superpixel-benchmark),[Paper](https://link.springer.com/article/10.1007/s11263-014-0744-2)|
| SH       |[Paper](https://ieeexplore.ieee.org/abstract/document/8360136/)|[SH](https://github.com/semiquark1/boruvka-superpixel)|
| SICLE    |[Paper]()|[SICLE](https://github.com/LIDS-UNICAMP/SICLE)|
| SLIC     |[Paper]()|[SLIC](https://www.epfl.ch/labs/ivrl/research/slic-superpixels/)|
| SNIC     |[Paper]()|[SNIC](https://github.com/achanta/SNIC)|

