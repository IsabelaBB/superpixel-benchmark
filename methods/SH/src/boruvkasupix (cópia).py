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
from statistics import mean

# local
sys.path.insert(0, "../pybuild")
sys.path.insert(0, "pybuild")
import boruvka_superpixel

def boruvkasupix(infile, outfile, superpixels):
    img_in = io.imread(infile)
    img_edge = np.zeros((img_in.shape[:2]), dtype=img_in.dtype)

    start_time = time.time()
    bosupix = boruvka_superpixel.BoruvkaSuperpixel()
    bosupix.build_2d(img_in, img_edge)
    final_time = (time.time() - start_time)
    
    if outfile is not None:
        for k in superpixels:
            out = bosupix.label(k)
            image = Image.fromarray(out)
            if len(superpixels) > 1:
                image.save(outfile + '/'+ str(k) + '/'+ (infile.split('/')[-1]).split('.')[0] + '.pgm')
            else:
                image.save(outfile + '/'+ (infile.split('/')[-1]).split('.')[0] + '.pgm')

    return final_time

        

def runDirectory(img_path, label_path, dtime_path, time_file, superpixels=[25,50,75,100,200,300,400,500,600,700,800,900,1000]):
    
    if os.path.isdir(img_path):
        img_timeList = []
        out_dtime_path = []
        for image_file in os.listdir(img_path):
            
            if os.path.isfile(os.path.join(img_path, image_file)):
                img_time = boruvkasupix(img_path+'/'+image_file, label_path, superpixels)
                
                img_timeList.append(img_time)
                if(dtime_path is not None):
                    out_dtime_path.append("{} {:5f}\n".format(image_file, img_time))

        for k in superpixels:
            if(dtime_path is not None):
                if(len(superpixels) > 1):
                    with open(dtime_path+'/SH-'+str(k)+'.txt', 'a') as f:
                        f.writelines(out_dtime_path)
                else:
                    with open(dtime_path, 'a') as f:
                        f.writelines(out_dtime_path)

        if(time_file is not None):
            with open(time_file, 'a') as f:
                for k in superpixels:
                    f.write("{} {:6f}\n".format(k, mean(img_timeList)))

    else:
        img_time = boruvkasupix(img_path, label_path)


def parse_arguments(argv):
    description = ('calculate superpixels, '
            'output orig image with color averaged within superpixels')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--img', type=str, help='input image file')
    parser.add_argument('--label', type=str, default=None, help='output image file')
    parser.add_argument('--k', type=int, help='number of superpixels')
    parser.add_argument('--dtime', type=str, default=None)
    parser.add_argument('--time', type=str, default=None)
    
    args = parser.parse_args(argv)
    return args

def main():
    args = parse_arguments(sys.argv[1:])
    runDirectory(args.img, args.label, args.dtime, args.time, [args.k])

if __name__ == '__main__':
    sys.exit(main())

