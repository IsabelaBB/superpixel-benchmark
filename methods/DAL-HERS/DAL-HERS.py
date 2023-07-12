#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@author: Hankui Peng

"""

## import necessary modules 
# system 
from skimage.segmentation import mark_boundaries
from glob import glob 
import numpy as np
import argparse 
import random 
import torch
from time import perf_counter
import sys 
import cv2
import os 
from PIL import Image

# local
sys.path.insert(0, "../pybuild")
sys.path.insert(0, "./pybuild")
import hers_superpixel

#from utils.analysis_util import *
#from utils.hed_edges import *
#from model.network import *

sys.path.insert(0, "./utils")
import analysis_util
from hed_edges import *

sys.path.insert(0, "./model")
from network import *


## input arguments 
parser = argparse.ArgumentParser(description='DAL-HERS Superpixel Segmentation on a folder of images',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--pretrained', default='./pretrained/DAL_loss=bce-rgb_date=23Feb2021.tar', help='path to the pretrained model')
parser.add_argument('--img', default=None, help='path to images folder')
parser.add_argument('--label', default=None, help='The path of the labeled image')
parser.add_argument('--border', default=None, help='The path of the borders image')
parser.add_argument('--k', default=None, help='Number of superpixels')
parser.add_argument('--dtime', default=None, help='Time log')
parser.add_argument('--time', default=None, help='Time log')
parser.add_argument('--edge', default=False, help='whether to incorporate edge information')
parser.add_argument('--device', default=torch.device("cuda" if torch.cuda.is_available() else "cpu"), help='default device (CPU / GPU)')
args = parser.parse_args()

random.seed(100)

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
        image = cv2.imread(img_file)
        image = image.astype(np.float32)

        # load the model 
        network_data = torch.load(args.pretrained, map_location=args.device)
        model = DAL(nr_channel=8, conv1_size=7)
        model.load_state_dict(network_data['state_dict'])
        model.to(args.device)
        model.eval()

        start = perf_counter()
        bosupix = DALHERS(image, args.edge, img_file, model)
        end = perf_counter()
        elapsed_time = (end - start)

        for nC in nC_seq:
            # obtain the segmentation labels 
            sp_label = bosupix.label(nC)
           
            # Borders
            if args.border is not None:
                border_img = np.max(image) * mark_boundaries(image, sp_label.astype(int), color = (0,0,255)) # candidate color choice: (220,20,60)
                spixl_save_name = args.border
                
                if os.path.isdir(spixl_save_name):
                    if len(nC_seq) > 1:
                        spixl_save_name = spixl_save_name + '/' + str(nC)
                    spixl_save_name = spixl_save_name + '/' + imgId + '.png'
                cv2.imwrite(spixl_save_name, border_img)
            
            # Labels
            if(args.label is not None):
                spixl_save_name = args.label 
                
                if(os.path.isdir(spixl_save_name)):
                    if len(nC_seq) > 1:
                        spixl_save_name = spixl_save_name + '/' + str(nC)
                    spixl_save_name = spixl_save_name + '/' + imgId + '.pgm'
            
                image = Image.fromarray(sp_label)
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
     

def DALHERS(img, edge, img_file, model):
    data_type = np.float32
    h, w, ch = img.shape

    ## input affinities 
    affinities = analysis_util.ProduceAffMap(img_file, model)
        
    ## HED edge information 
    if edge:
        Input = torch.FloatTensor(np.ascontiguousarray(np.array(Image.open(img_file))[:, :, ::-1].transpose(2, 0, 1).astype(np.float32)))*(1.0 / 255.0)
        edge_prob = estimate(Input) 
        input_edge = np.array(edge_prob.squeeze(0), dtype=data_type)
    else:
        input_edge = np.ones((h, w), dtype=data_type) # Provide no external edge information by default 
        
    ## build the hierarchical segmentation tree 
    bosupix = hers_superpixel.BoruvkaSuperpixel()
    bosupix.build_2d(img, affinities, input_edge)
    return bosupix

            
if __name__ == '__main__':
    runDirectory()
