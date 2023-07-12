#!/usr/bin/env python3

# standard lib
import sys
import argparse

# numpy family
import numpy as np

# 3rd party
from skimage import io
from PIL import Image

import os
import time
from glob import glob 

# local
sys.path.insert(0, "../pybuild")
sys.path.insert(0, "pybuild")
import boruvka_superpixel

description = ('calculate superpixels, '
            'output orig image with color averaged within superpixels')
parser = argparse.ArgumentParser(description=description)
parser.add_argument('--img', type=str, help='input image file')
parser.add_argument('--label', type=str, default=None, help='output image file')
parser.add_argument('--k', type=int, help='number of superpixels')
parser.add_argument('--dtime', type=str, default=None)
parser.add_argument('--time', type=str, default=None)

args = parser.parse_args()


def runDirectory():
    
    if(os.path.isdir(args.img)):
        # read all the image files in the folder 
        img_files = glob(args.img+'/*')
        img_files.sort()
    else:
        img_files = [args.img]

    if(args.k is None):
        nC_seq = [25,50,75,100,200,300,400,500,600,700,800,900,1000]
    else:
        nC_seq = [int(args.k)]
    
    overall_time = {}
    for nC in nC_seq:
        overall_time[nC] = 0

    # for each image:
    for n in range(len(img_files)):
        
        ## input image 
        img_file = img_files[n]
        imgId = os.path.basename(img_file).split('.')[0]
        img_in = io.imread(img_file)
        img_edge = np.zeros((img_in.shape[:2]), dtype=img_in.dtype)

        start_time = time.time()
        bosupix = boruvka_superpixel.BoruvkaSuperpixel()
        bosupix.build_2d(img_in, img_edge)
        elapsed_time = (time.time() - start_time)

        for nC in nC_seq:
            image = Image.fromarray(bosupix.label(nC))
            
            # Labels
            if(args.label is not None):
                spixl_save_name = args.label 
                
                if(os.path.isdir(spixl_save_name)):
                    if len(nC_seq) > 1:
                        spixl_save_name = spixl_save_name + '/' + str(nC)
                    spixl_save_name = spixl_save_name + '/' + imgId + '.pgm'
            
                image.save(spixl_save_name)


            # save the run times 
            if(args.dtime is not None):
                time_save_name = args.dtime
                if os.path.isdir(time_save_name):
                    time_save_name = time_save_name + "-" + str(nC) + ".txt"
                with open(time_save_name, "a+") as fileTime:
                    fileTime.write(imgId + " " + "%.5f" % elapsed_time + "\n")
            
            overall_time[nC] += elapsed_time

    if(len(img_files) == 0):
        print("ERROR: no image founded")
        return

    if(args.time is not None):
        with open(args.time, "a+") as fileTime:
            for nC in nC_seq:
                fileTime.write(str(nC) + " " + "%.5f" % (overall_time[nC]/len(img_files)) + "\n") 



if __name__ == '__main__':
    runDirectory()

