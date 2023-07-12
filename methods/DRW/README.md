# DRW
## Overview
The official code of the paper: [Dynamic Random Walk for Superpixel Segmentation][arxiv] (ACCV2018/IEEE TIP)

## Getting Started
To run our code:

1. Run ``MakeiFile.m" ro generate the matlab function based on our "C++ codes"

2. Setting hyperparameters of ``exampleDRW.m":

2.1. Set "k" as the superpixel number

2.2. Set "type_seed" to select our seed initialize strategy

2.3. Set "type_grad" to select seed prior (if use pb, please firstly download the source code of pb edge detection)

3. Run the code.


## Citation

If you find our work useful in your research, please cite:

@inproceedings\{DRW_ACCV,</br>
  title=\{Dynamic Random Walk for Superpixel Segmentation\},</br>
  author=\{Zhu, Lei and Kang, Xuejing and Ming, Anlong and Zhang, Xuesong\},</br>
  booktitle=\{Asian Conference on Computer Vision\},</br>
  pages=\{540--554\},</br>
  year=\{2018\},</br>
  organization=\{Springer\}</br>
\}</br>

@article\{DRW_TIP,</br>
  title=\{Dynamic random walk for superpixel segmentation\},</br>
  author=\{Kang, Xuejing and Zhu, Lei and Ming, Anlong\},</br>
  journal=\{IEEE Transactions on Image Processing\},</br>
  volume=\{29\},</br>
  pages=\{3871--3884\},</br>
  year=\{2020\},</br>
  publisher=\{IEEE\}</br>
\}</br>

[arxiv]: https://ieeexplore.ieee.org/abstract/document/8967213/
