import argparse
import os
import torch.backends.cudnn as cudnn
import models
import torchvision.transforms as transforms
import flow_transforms
#from scipy.misc import imread
#from scipy.misc import imsave
from loss import *
import time
import random
from glob import glob
from PIL import Image
from math import sqrt, ceil
import numpy as np

#import matplotlib.pyplot as plt

# import sys
# sys.path.append('../cython')
# from connectivity import enforce_connectivity

'''
Infer from custom dataset:
author:Fengting Yang 
last modification: Mar.5th 2020

usage:
1. set the ckpt path (--pretrained) and output
2. comment the output if do not need

results will be saved at the args.label

'''

#os.environ['CUDA_VISIBLE_DEVICES'] = '1'

model_names = sorted(name for name in models.__dict__
                     if name.islower() and not name.startswith("__"))


parser = argparse.ArgumentParser(description='PyTorch SPixelNet inference on a folder of imgs',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--img', metavar='DIR', default='./demo/inputs', help='path to images folder')
parser.add_argument('--pretrained', metavar='PTH', help='path to pre-trained model',
                                    default= './pretrain_ckpt/SpixelNet_bsd_ckpt.tar')
parser.add_argument('--label', metavar='DIR', default= './demo' , help='path to output folder')

parser.add_argument('--k', default=25, type=int,help='number of superpixels')
parser.add_argument('--downsize', default=16, type=float,help='superpixel grid cell, must be same as training setting')

parser.add_argument('--dtime', default=None, help='The path of the time log')
parser.add_argument('--time', default=None, help='The path of the time log')


args = parser.parse_args()

def getScale(H_, W_, spx):
    return float(sqrt((spx*args.downsize*args.downsize)/(H_*W_)))

random.seed(100)
@torch.no_grad()
def test(args, model, img_paths, save_path, idx):
    # Data loading code
    input_transform = transforms.Compose([
        flow_transforms.ArrayToTensor(),
        transforms.Normalize(mean=[0,0,0], std=[255,255,255]),
        transforms.Normalize(mean=[0.411,0.432,0.45], std=[1,1,1])
    ])

    img_file = img_paths[idx]
    load_path = img_file
    imgId = os.path.basename(img_file)[:-4]

    # may get 4 channel (alpha channel) for some format
    img_ = cv2.imread(load_path)[:, :, :3]

    H_ = img_.shape[0]
    W_ = img_.shape[1]
    
    scale = getScale(H_, W_, args.num_spx)
    new_height, new_width = ceil(((float(H_) * scale)/args.downsize))*args.downsize, ceil(((float(W_) * scale)/args.downsize))*args.downsize
    
    img = cv2.resize(img_, (int(new_width), int(new_height)), interpolation=cv2.INTER_CUBIC)

    # get spixel id
    n_spixl_h = int(np.floor(new_height / args.downsize))
    n_spixl_w = int(np.floor(new_width / args.downsize))

    spix_values = np.int32(np.arange(0, n_spixl_w * n_spixl_h).reshape((n_spixl_h, n_spixl_w)))
    spix_idx_tensor_ = shift9pos(spix_values)

    spix_idx_tensor = np.repeat(
      np.repeat(spix_idx_tensor_, args.downsize, axis=1), args.downsize, axis=2)

    spixeIds = torch.from_numpy(np.tile(spix_idx_tensor, (1, 1, 1, 1))).type(torch.float).cuda()

    n_spixel =  int(n_spixl_h * n_spixl_w)

    img1 = input_transform(img)
    ori_img = input_transform(img_)

    # compute output
    tic = time.time()
    output = model(img1.cuda().unsqueeze(0))
    
    # assign the spixel map
    curr_spixl_map = update_spixl_map(spixeIds, output)
    ori_sz_spixel_map = F.interpolate(curr_spixl_map.type(torch.float), size=(H_, W_), mode='nearest').type(torch.int)

    mean_values = torch.tensor([0.411, 0.432, 0.45], dtype=img1.cuda().unsqueeze(0).dtype).view(3, 1, 1)
    spixel_viz, spixel_label_map = get_spixel_image((ori_img + mean_values).clamp(0, 1), ori_sz_spixel_map.squeeze(), n_spixels= int((H_ * scale * W_ * scale) / (args.downsize * args.downsize)),  b_enforce_connect=True)
    
    torch.cuda.synchronize()
    toc = time.time() - tic

    # ************************ Save all result********************************************
    # save img, uncomment it if needed
    # if not os.path.isdir(os.path.join(save_path, 'img')):
    #     os.makedirs(os.path.join(save_path, 'img'))
    # spixl_save_name = os.path.join(save_path, 'img', imgId + '.jpg')
    # img_save = (ori_img + mean_values).clamp(0, 1)
    # imsave(spixl_save_name, img_save.detach().cpu().numpy().transpose(1, 2, 0))

    spixel_viz = spixel_viz.transpose(1, 2, 0)
    spixel_label_map = spixel_label_map.astype(np.uint16)
    # save spixel viz
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    spixl_save_name = os.path.join(save_path, imgId + '.png')
    
    image = Image.fromarray(spixel_label_map)
    image.save(spixl_save_name)

    if idx % 10 == 0:
        print("processing %d"%idx)

    return toc

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

    save_path = args.label
    print('=> will save everything to {}'.format(save_path))
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    print('{} samples found'.format(len(img_files)))

    overall_time = 0.0

    # create model
    network_data = torch.load(args.pretrained)
    print("=> using pre-trained model '{}'".format(network_data['arch']))
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
      
        H_ = img_.shape[0]
        W_ = img_.shape[1]
        if(len(img_.shape) < 3):
            img_ = np.stack((img_,)*3, axis=-1)

        scale = getScale(H_, W_, args.k)
        new_height, new_width = ceil(((float(H_) * scale)/args.downsize))*args.downsize, ceil(((float(W_) * scale)/args.downsize))*args.downsize
        
        img = cv2.resize(img_, (int(new_width), int(new_height)), interpolation=cv2.INTER_CUBIC)
        
        # get spixel id
        n_spixl_h = int(np.floor(new_height / args.downsize))
        n_spixl_w = int(np.floor(new_width / args.downsize))
        
        spix_values = np.int32(np.arange(0, n_spixl_w * n_spixl_h).reshape((n_spixl_h, n_spixl_w)))
        spix_idx_tensor_ = shift9pos(spix_values)
        spix_idx_tensor = np.repeat(
            np.repeat(spix_idx_tensor_, args.downsize, axis=1), args.downsize, axis=2)
        spixeIds = torch.from_numpy(np.tile(spix_idx_tensor, (1, 1, 1, 1))).type(torch.float).cuda()

        n_spixel =  int(n_spixl_h * n_spixl_w)

        img1 = input_transform(img)
        ori_img = input_transform(img_)

        # compute output
        tic = time.time()
        output = model(img1.cuda().unsqueeze(0))

        # assign the spixel map
        curr_spixl_map = update_spixl_map(spixeIds, output)
        ori_sz_spixel_map = F.interpolate(curr_spixl_map.type(torch.float), size=(H_, W_), mode='nearest').type(torch.int)

        mean_values = torch.tensor([0.411, 0.432, 0.45], dtype=img1.cuda().unsqueeze(0).dtype).view(3, 1, 1)
        spixel_viz, spixel_label_map = get_spixel_image((ori_img + mean_values).clamp(0, 1), ori_sz_spixel_map.squeeze(), n_spixels= int((H_ * scale * W_ * scale) / (args.downsize * args.downsize)),  b_enforce_connect=True)

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
