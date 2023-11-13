import shutil
import os
import argparse

parser = argparse.ArgumentParser(description='PyTorch StereoSpixel inference on a folder of img pairs',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--src', default='../output/NYU', help='path to spixel test dir')
parser.add_argument('--dst', default= './pretrain_ckpt/NYU_only_test_padding_allLabloss_v1.tar', help='path to collect all the evaluation results',)

args = parser.parse_args()

src = args.src
dst = args.dst

# list = ["54" ,"96", "150" ,"216" ,"294", "384" ,"486", "600", "726" ,"864", "1014", "1176", "1350", "1536", "1944" ]
list = ["300", "432", "588", "768", "972", "1200", "1452", "1728", "2028", "2352"]
for l in list:
    src_pth = src + '/SPixelNet_nSpixel_' + l +'/map_csv/results.csv'
    dst_pth = dst + '/' + l
    if not os.path.isdir(dst_pth):
        os.makedirs(dst_pth)
    dst_path =dst_pth + '/results.csv'
    shutil.copy(src_pth, dst_path)
