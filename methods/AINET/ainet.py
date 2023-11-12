import argparse
import os

import torch.backends.cudnn as cudnn

import models
import torchvision.transforms as transforms
import flow_transforms
from PIL import Image
import numpy as np
from glob import glob
from math import sqrt, ceil

from loss import *
import time
import random

os.environ['CUDA_VISIBLE_DEVICES'] = '0'

model_names = sorted(name for name in models.__dict__
                     if name.islower() and not name.startswith("__"))

parser = argparse.ArgumentParser(description='PyTorch SPixelNet inference on a folder of imgs',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--img', metavar='DIR', default='/home/name/superpixel/dataset/NYUV2/', help='path to images folder')
parser.add_argument('--pretrained', metavar='PTH', help='path to pre-trained model', default= './pretrain_ckpt/SpixelNet_bsd_ckpt.tar')
parser.add_argument('--label', metavar='DIR', default= '' , help='path to output folder')

parser.add_argument('--k', type=int, default=100, help='Number of superpixels')

parser.add_argument('--dtime', default=None, help='The path of the time log')
parser.add_argument('--time', default=None, help='The path of the time log')

parser.add_argument('--same_size', type=bool, default=False, help='Wether the images have the same size')

parser.add_argument('--downsize', default=16, type=float,help='superpixel grid cell, must be same as training setting')
parser.add_argument('-nw', '--num_threads', default=1, type=int,  help='num_threads')
parser.add_argument('-b', '--batch-size', default=1, type=int, metavar='N', help='mini-batch size')

# 
parser.add_argument('--input_img_height', '-v_imgH', default=480,  type=int, help='img height_must be 16*n')
parser.add_argument('--input_img_width', '-v_imgW', default=640,   type=int, help='img width must be 16*n')

args = parser.parse_args()
#args.test_list = args.data_dir + '/nyuv2_test_subset.txt'


random.seed(100)
@torch.no_grad()

def getScale(H_, W_, spx):
    return float(sqrt((spx*16.0*16.0)/(H_*W_)))
    

def main():
    global args, save_path
    print("=> fetching img pairs in '{}'".format(args.img))

    if(os.path.isdir(args.img)):
      # read all the image files in the folder 
      img_files = glob(args.img+'/*')
      img_files.sort()
    else:
      img_files = [args.img]
    print('{} samples found'.format(len(img_files)))

    overall_time = 0.0

    # create model
    network_data = torch.load(args.pretrained)
    print("=> using pre-trained model '{}'".format(args.pretrained))
    model = models.__dict__[network_data['arch']]( data = network_data).cuda()
    model.eval()
    args.arch = network_data['arch']
    cudnn.benchmark = True

    if args.label is not None:
      if not os.path.isdir(args.label):
        os.makedirs(args.label)
      save_path = args.label #+ '/{0}'.format(args.k)
      print('=> will save everything to {}'.format(save_path))
    
    for n in range(len(img_files)):
      img_file = img_files[n]
      img_name = os.path.basename(img_file).split('.')[0]
      
      # Data loading code
      input_transform = transforms.Compose([
        flow_transforms.ArrayToTensor(),
        transforms.Normalize(mean=[0,0,0], std=[255,255,255]),
        transforms.Normalize(mean=[0.411,0.432,0.45], std=[1,1,1])
      ])

      
      img_ = np.asarray(Image.open(img_file))
      
      if args.label is not None and not os.path.isdir(save_path):
        os.makedirs(save_path)

      H_ = img_.shape[0]
      W_ = img_.shape[1]

      if(len(img_.shape) < 3):
        img_ = np.stack((img_,)*3, axis=-1)

      scale = getScale(H_, W_, args.k)
      args.input_img_height, args.input_img_width = ceil(((float(H_) * scale)/16.0))*16.0, ceil(((float(W_) * scale)/16.0))*16.0
      spixlId, XY_feat_stack = init_spixel_grid(args, b_train=False)
      
      img = cv2.resize(img_, (int(args.input_img_width), int(args.input_img_height)), interpolation=cv2.INTER_CUBIC)
      img1 = input_transform(img)
      ori_img = input_transform(img_)

      # compute output
      tic = time.time()
      output,_ = model(img1.cuda().unsqueeze(0))
      
      # assign the spixel map
      curr_spixl_map = update_spixl_map(spixlId, output)
      ori_sz_spixel_map = F.interpolate(curr_spixl_map.type(torch.float), size=(H_, W_), mode='nearest').type(torch.int)
      
      mean_values = torch.tensor([0.411, 0.432, 0.45], dtype=img1.cuda().unsqueeze(0).dtype).view(3, 1, 1)
      spixel_viz, spixel_label_map = get_spixel_image((ori_img + mean_values).clamp(0, 1), ori_sz_spixel_map.squeeze(), n_spixels= int((H_ * scale * W_ * scale) / (16 * 16)),  b_enforce_connect=True)

      torch.cuda.synchronize()
      toc = time.time() - tic
      overall_time += toc

      # save spixel viz
      if(args.label is not None):
        spixl_save_name = os.path.join(save_path, img_name + '.png')
        image = Image.fromarray((spixel_label_map + 1).astype(np.uint16))
        image.save(spixl_save_name)

      if (args.dtime is not None):
        with open(args.dtime, "a+") as fileTime:
          fileTime.write(img_name + " " + "%.5f" % toc + "\n")

    if (args.time is not None):
      with open(args.time, "a+") as fileTime:
        fileTime.write(str(args.k) + " " + "%.5f" % (overall_time/len(img_files)) + "\n")


if __name__ == '__main__':
    main()
