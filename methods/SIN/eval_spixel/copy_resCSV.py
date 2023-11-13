import shutil
import os
import argparse

parser = argparse.ArgumentParser(description='PyTorch StereoSpixel inference on a folder of img pairs',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--src', default='/media/yuanqing/ssd/code/visual_tracking/superpixel_deconv/output/NYU/finetune', help='path to spixel test dir')
parser.add_argument('--dst', default= './result_set/deconv_finetune/nyu', help='path to collect all the evaluation results',)

args = parser.parse_args()

src = args.src
dst = args.dst

# list = ["54", "96", "150", "216" ,"294", "384" ,"486", "600", "726" ,"864", "1014", "1176", "1350", "1536", "1944" ]
list = ["300", "432", "588", "768", "972", "1200", "1452", "1728", "2028", "2352"]
# list = ["200", "300", "400", "500", "600", "700", "800", "900", "1000", "1100", "1200"]
# list = ["300", "500", "700", "900", "1100", "1300", "1500", "1700", "1900", "2100", "2300"]
# list = ["600"]
for l in list:
    src_pth = src + '/SPixelNet_nSpixel_' + l +'/map_csv/results.csv'
    dst_pth = dst + '/' + l
    if not os.path.isdir(dst_pth):
        os.makedirs(dst_pth)
    dst_path =dst_pth + '/results.csv'
    shutil.copy(src_pth, dst_path)
