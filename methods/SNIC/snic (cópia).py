from PIL import Image
import numpy as np
import skimage.color
from pysnic.algorithms.snic import snic

import argparse
import os
import time
from statistics import mean
import sys

## input arguments 
def parse_arguments(argv):
    description = ('calculate superpixels, '
            'output orig image with color averaged within superpixels')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--img', type=str, help='input image file')
    parser.add_argument('--label', type=str, default=None, help='output image file')
    parser.add_argument('--k', type=int, help='number of superpixels')
    parser.add_argument('--dtime', type=str, default=None)
    parser.add_argument('--time', type=str, default=None)
    parser.add_argument('--cmp', type=float, default=10.00, help='compactness')
    
    args = parser.parse_args(argv)
    return args


def runDirectory(img_path, num_spx, label_path, dtime_file, time_file, compactness):
    
    if os.path.isdir(img_path):
        img_timeList = []
        for image_file in os.listdir(img_path):
            
            if(dtime_file is not None):
                f = open(dtime_file, 'a')

            if os.path.isfile(os.path.join(img_path, image_file)):
                label_file = label_path
                if(label_path is not None):
                    label_file = label_path+'/'+image_file.split('.')[0]+'.png'
                img_time = runSNIC(img_path+'/'+image_file, num_spx, label_file, compactness)
                
                if(dtime_file is not None):
                    f.write("{} {:5f}\n".format(image_file, img_time))
                img_timeList.append(img_time)

            if(dtime_file is not None):
                f.close()

        if(time_file is not None):
            with open(time_file, 'a') as f:
                f.write("{} {:5f}\n".format(num_spx, mean(img_timeList)))

    else:
        img_time = runSNIC(img_path, num_spx, label_path)


def runSNIC(img_file, num_spx, label_file, compactness):

    # load image
    rgb_image = np.array(Image.open(img_file))
    if(len(rgb_image.shape) == 2):
        rgb_image = np.stack((rgb_image,)*3, axis=-1)

    # convert the image from RGB to CIELAB
    lab_image = skimage.color.rgb2lab(rgb_image)

    # SNIC
    # lab_image should be a list for 'snic'
    lab_image = lab_image.tolist()

    start_time = time.time()
    segmentation, _, centroids = snic(lab_image, num_spx, compactness)
    final_time = (time.time() - start_time)
    
    if(label_file is not None):
        segmentation = np.array(segmentation).astype(np.uint16)
        image = Image.fromarray(segmentation)
        image.save(label_file)

    return final_time

def main():
    args = parse_arguments(sys.argv[1:])
    runDirectory(args.img, args.k, args.label, args.dtime, args.time, args.cmp)

if __name__ == '__main__':
    sys.exit(main())

