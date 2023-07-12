## SCALP: Superpixels with Contour Adherence using Linear Path

### Overview

Implementation of paper:   [PDF](https://hal.archives-ouvertes.fr/hal-01510063/document)
```
@article{giraud2018scalp,
  title={Robust superpixels using color and contour features along linear path},
  author={Giraud, R{\'e}mi and Ta, Vinh-Thong and Papadakis, Nicolas},
  journal={Computer Vision and Image Understanding},
  volume={170},
  pages={1--13},
  year={2018}
}
```

- SCALP superpixels generated using color and contour features on the linear path between the pixel and the superpixel barycenter: 

![image](./Figures/scalp_method.png)


### Requirements

- Linux

- For C++ version:  [CImg library](http://cimg.eu/)  (unique .h already included)

- A contour prior map can be provided to our method (for C++: an image with 255 for highest contour intensity)  
The contour detection method used in this work is available [here](https://github.com/pdollar/edges)  
Other contour detection methods can be found [here](https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html)



### Execution

#### MATLAB / C-Mex
```
run main.m    %call SCALP_mex
```


#### C++

- Compilation:
```
make
```

- Execution prototype:
```
./SCALP -i img_name [-k superpixel_nbr(450)] [-m compactness(0.075)]  [-outm output_map_name(./res/labelmap.png)] [-outb output_border_name(./res/borders.png)]  [-c contour(NULL)]
```
- Execution with contour map:  (make test)
``` 
./SCALP -i ./data/test_img.jpg -k 300 -m 0.075 -outm test_img_map.png -outb test_img_border.png -c ./data/test_img_contour.png
```
- Execution on an image list:  (make test_list)
```
./scripts/test_list.sh ./data/list_file.txt ./data/ 450 0.075
```


### Data

The Berkeley Segmentation Dataset (BSD) containing 500 images of size 321x481 pixels with segmentation and contour ground truths is available 
[here](https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html)

### License

(C) Rémi Giraud, 2020  
remi.giraud@u-bordeaux.fr  
[https://rgiraud.vvv.enseirb-matmeca.fr](https://rgiraud.vvv.enseirb-matmeca.fr)  
ENSEIRB-MATMECA (Bordeaux INP), Laboratory IMS

This code is free to use, share and modify for any non-commercial purposes.
Any commercial use is strictly prohibited without the authors' consent, also including authors from (chen2017) since SCALP uses some part of their code:
```
@InProceedings{chen2017,
    author   = {Jiansheng Chen and Zhengqin Li and Bo Huang},
    title    = {Linear Spectral Clustering Superpixel},
    booktitle  = {IEEE Transactions on Image Processing},
    YEAR   = {2017},
}
```


