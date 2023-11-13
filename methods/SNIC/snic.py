from PIL import Image
import numpy as np
import skimage.color
from pysnic.algorithms.snic import snic

import argparse
import os
import time
import sys
from glob import glob 

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


def runDirectory(args):
    
    if(os.path.isdir(args.img)):
        # read all the image files in the folder 
        img_files = glob(args.img+'/*')
        img_files.sort()
    else:
        img_files = [args.img]

    overall_time = 0.0

    # for each image:
    for n in range(len(img_files)):
        img_file = img_files[n]
        img_name = os.path.basename(img_file).split('.')[0]
        
        # load image
        rgb_image = np.array(Image.open(img_file))
        if(len(rgb_image.shape) == 2):
            rgb_image = np.stack((rgb_image,)*3, axis=-1)

        # convert the image from RGB to CIELAB
        lab_image = skimage.color.rgb2lab(rgb_image)
        # lab_image should be a list for 'snic'
        lab_image = lab_image.tolist()

        # SNIC
        start_time = time.time()
        segmentation, _, centroids = snic(lab_image, args.k, args.cmp)
        final_time = (time.time() - start_time)
        overall_time += final_time

        # save segmentation labels
        if(args.label is not None):
            segmentation = np.array(segmentation).astype(np.uint16)
            image = Image.fromarray(segmentation)
            spixl_save_name = args.label 
            if(os.path.isdir(spixl_save_name)):
                spixl_save_name = spixl_save_name + '/' + img_name + '.png'
            image.save(spixl_save_name)

        # save the run times 
        if(args.dtime is not None):
            with open(args.dtime, "a+") as fileTime:
                fileTime.write(img_name + " " + "%.5f" % final_time + "\n")
    
        if(len(img_files) == 0):
            print("ERROR: no image founded")
            return
        
    if(args.time is not None):
        with open(args.time, "a+") as fileTime:
            fileTime.write(str(args.k) + " " + "%.5f" % (overall_time/len(img_files)) + "\n") 
     
def main():
    args = parse_arguments(sys.argv[1:])
    runDirectory(args)

if __name__ == '__main__':
    sys.exit(main())

