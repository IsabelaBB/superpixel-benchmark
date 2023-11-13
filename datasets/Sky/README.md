# Datasets

We selected five datasets which impose different challenges for superpixel segmentation: _Birds_[1]; _Insects_[1]; _Sky_[2]; _ECSSD_[3]; and _NYUV2_[4]. The Table bellow summarizes their main characteristics. 

|                    	|  Birds  	| Insects 	|   Sky   	|       ECSSD       	|  NYUV2  	|
|:------------------:	|:-------:	|:-------:	|:-------:	|:-----------------:	|:-------:	|
|      Image content 	| Natural 	| Natural 	| Natural 	| Natural and urban 	|  Indoor 	|
|   Number of images 	|   150   	|   130   	|    60   	|        1000       	|   1449  	|
| Minimum image size 	| 300x300 	| 640x359 	| 599x399 	|      400x139      	| 608x448 	|
| Maximum image size 	| 640x640 	| 640x640 	| 825x600 	|      400x400      	| 608x448 	|

Birds, Insects, and Sky have natural images. The former contains images of birds, which have thin and elongated parts, hard to delineate with compact superpixels. When there are more birds, they often overlap, making it difficult to delineate them separately. In the Birds dataset, the background does not have a specific pattern and may (or may not) be blurred, colored, and textured. Similarly, the Insects dataset has images containing one or more insects with thin and elongated parts. Compared to the Birds dataset, it has more blurred and less textured backgrounds, and their objects (the insects) have thinner parts. Therefore, it has more challenging objects but less difficult backgrounds than the Birds dataset. In contrast, the Sky dataset has images with one plane each. Most images in the Sky dataset have large regions with low color and subtle luminosity variations. The ground truth in Birds, Insects, and Sky datasets are binary masks with only one connected object per image. The objects in Birds, Insects, and Sky, are a bird, an insect, and the sky, respectively. 

In contrast with the aforementioned datasets, ECSSD and NYUV2 datasets have urban scenes. Specifically, the ECSSD dataset has images in both natural and urban environments, while NYUV2 is composed of video sequences from several indoor scenes recorded by Microsoft Kinect. The images in the ECSSD dataset have complex scenes, most with non-uniform regions and backgrounds composed of several parts. In the ECSSD dataset, the images may have more objects, many without well-defined boundaries. Furthermore, some objects have transparency, which makes them difficult to identify. Conversely, the RGBD images in the NYUV2 dataset have rich geometric structures with large planar surfaces, such as the floor, walls, and table tops. Its images also have small objects and occlusion, accentuated by the mess and disorder common in inhabited environments. In the ECSSD dataset the ground truth images are binary masks, each one with at least one connected object. In the NYUV2 dataset, the ground truth images have dense multi-class labels. However, for those in the NYUV2 dataset, we remove unlabeled pixels similar to [5]. 

## References

[1] Lucy A. C. Mansilla and Paulo A. V. Miranda. 2016. Oriented Image Foresting Transform Segmentation: Connectivity Constraints with Adjustable Width. In 2016 29th SIBGRAPI Conference on Graphics, Patterns and Images (SIBGRAPI).
289–296. https://doi.org/10.1109/SIBGRAPI.2016.047

[2] Eduardo Barreto Alexandre, Ananda Shankar Chowdhury, Alexandre Xavier Falcao, and Paulo A. Vechiatto Miranda.
2015. IFT-SLIC: A General Framework for Superpixel Generation Based on Simple Linear Iterative Clustering
and Image Foresting Transform. In 2015 28th SIBGRAPI Conference on Graphics, Patterns and Images. 337–344.
https://doi.org/10.1109/SIBGRAPI.2015.20

[3] Jianping Shi, Qiong Yan, Li Xu, and Jiaya Jia. 2015. Hierarchical image saliency detection on extended CSSD. IEEE
transactions on pattern analysis and machine intelligence 38, 4 (2015), 717–729.

[4] SILBERMAN, Nathan et al. Indoor segmentation and support inference from rgbd images. In: Computer Vision–ECCV 2012: 12th European Conference on Computer Vision, Florence, Italy, October 7-13, 2012, Proceedings, Part V 12. Springer Berlin Heidelberg, 2012. p. 746-760.

[5] STUTZ, David; HERMANS, Alexander; LEIBE, Bastian. Superpixels: An evaluation of the state-of-the-art. Computer Vision and Image Understanding, v. 166, p. 1-27, 2018.
