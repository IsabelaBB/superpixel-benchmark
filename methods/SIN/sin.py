import argparse
import os
import torch.backends.cudnn as cudnn
import models
import torchvision.transforms as transforms
import flow_transforms
from imageio import imread, imsave
import numpy as np
from PIL import Image
from loss import *
import time
import random
from glob import glob
from models.model_util import update_spixel_map

from math import sqrt, ceil, floor #, max, min

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

results will be saved at the args.output

'''

# os.environ['CUDA_VISIBLE_DEVICES'] = '1'

model_names = sorted(name for name in models.__dict__
                     if name.islower() and not name.startswith("__"))

parser = argparse.ArgumentParser(description='PyTorch SPixelNet inference on a folder of imgs',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--img', metavar='DIR', default='./demo/inputs', help='path to images folder')
parser.add_argument('--pretrained', metavar='PTH', help='path to pre-trained model',
                                    default= './pretrain_ckpt/model_best.tar')
parser.add_argument('--label', metavar='DIR', default= './demo' , help='path to output folder')

parser.add_argument('--k', default=25, type=int,help='number of superpixels')
parser.add_argument('--downsize', default=16, type=float,help='superpixel grid cell, must be same as training setting')

parser.add_argument('--dtime', default=None, help='The path of the time log')
parser.add_argument('--time', default=None, help='The path of the time log')

args = parser.parse_args()

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
    imgId = os.path.basename(img_file).split('.')[0]

    # may get 4 channel (alpha channel) for some format
    img_ = imread(load_path)[:, :, :3]
    H, W, _ = img_.shape
    H_, W_  = int(np.ceil(H/16.)*16-15), int(np.ceil(W/16.)*16-15)

    img1 = cv2.resize(img_, (W_, H_), interpolation=cv2.INTER_CUBIC)
    img1 = input_transform(img1)
    ori_img = input_transform(img_)

    # compute output
    tic = time.time()
    prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h = model(img1.cuda().unsqueeze(0))
    toc = time.time() - tic

    # assign the spixel map
    curr_spixl_map = update_spixel_map(img1.cuda().unsqueeze(0), prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h)
    # curr_spixl_map = map0
    ori_sz_spixel_map = F.interpolate(curr_spixl_map.type(torch.float), size=( H_,W_), mode='nearest').type(torch.int)

    mean_values = torch.tensor([0.411, 0.432, 0.45], dtype=img1.cuda().unsqueeze(0).dtype).view(3, 1, 1)
    spixel_viz, spixel_label_map = get_spixel_image((ori_img + mean_values).clamp(0, 1), ori_sz_spixel_map.squeeze(), n_spixels= 0,  b_enforce_connect=False)

    # ************************ Save all result********************************************
    # save img, uncomment it if needed
    # if not os.path.isdir(os.path.join(save_path, 'img')):
    #     os.makedirs(os.path.join(save_path, 'img'))
    # spixl_save_name = os.path.join(save_path, 'img', imgId + '.jpg')
    # img_save = (ori_img + mean_values).clamp(0, 1)
    # imsave(spixl_save_name, img_save.detach().cpu().numpy().transpose(1, 2, 0))


    # save spixel viz
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    spixl_save_name = os.path.join(save_path, '/', imgId + '.png')
    imsave(spixl_save_name, spixel_viz.transpose(1, 2, 0))

    # save the unique maps as csv, uncomment it if needed
    #if not os.path.isdir(os.path.join(save_path, 'map_csv')):
    #    os.makedirs(os.path.join(save_path, 'map_csv'))
    #output_path = os.path.join(save_path, 'map_csv', imgId + '.csv')
      # plus 1 to make it consistent with the toolkit format
    #np.savetxt(output_path, spixel_label_map.astype(int), fmt='%i',delimiter=",")


    if idx % 10 == 0:
        print("processing %d"%idx)

    return toc


def getSpxVals_bkp(H_, W_, spx):
  #divisors = 0
  ref_scale = min(H_, W_) / max(H_,W_)
  min_diff = spx
  best_val = sqrt(spx)
  for i in range(2, int(floor(sqrt(spx)))):
    if(spx % i == 0):
      #divisors += 1
      if(abs(ref_scale - (i / (spx / i))) < min_diff):
        min_diff = abs(ref_scale - (i / (spx / i)))
        best_val = spx/i
  
  return 12, 9
  if(H_ > W_):
     return best_val, spx/best_val
  else:
      return spx/best_val, best_val


def nextStep(step, x, y):
  if(step == 1):
    return 2, x+1, y
  if(step == 2 or step == 4):
    return step+1, x-1, y-1
  if(step == 3 or step == 5):
    return step+1, x+2, y+1
  if(step == 6):
    return 1, x-1, y-1

def getSpxVals(H_, W_, spx):
  x,y = int(ceil(sqrt(spx))), int(floor(sqrt(spx)))
  rel = max(H_,W_)/min(H_,W_)
  min_diff_rel = rel - x/y
  step, next_x, next_y = nextStep(1, x, y)

  while(abs(rel - (next_x/next_y)) > min_diff_rel):
    x, y = next_x, next_y
    min_diff_rel = abs(rel - (x/y))
    step, next_x, next_y = nextStep(step, x, y)

  if(H_ > W_): return x,y
  else: return y,x


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

        #H_, W_  = int(np.ceil(H/16.)*16-15), int(np.ceil(W/16.)*16-15)
        #scale = getScale(H_, W_, args.k)
        #new_height, new_width = ceil(((float(H_) * scale)/args.downsize))*args.downsize, ceil(((float(W_) * scale)/args.downsize))*args.downsize
        
        h, w = getSpxVals(H_, W_, args.k)
        new_height, new_width  = int(np.ceil(h)*16-15), int(np.ceil(w)*16-15)

        img = cv2.resize(img_, (int(new_width), int(new_height)), interpolation=cv2.INTER_CUBIC)
        
        img1 = input_transform(img)
        ori_img = input_transform(img_)

        # compute output
        tic = time.time()
        prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h = model(img1.cuda().unsqueeze(0))

        # assign the spixel map
        curr_spixl_map = update_spixel_map(img1.cuda().unsqueeze(0), prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h)
        ori_sz_spixel_map = F.interpolate(curr_spixl_map.type(torch.float), size=(H_, W_), mode='nearest').type(torch.int)

        mean_values = torch.tensor([0.411, 0.432, 0.45], dtype=img1.cuda().unsqueeze(0).dtype).view(3, 1, 1)
        #spixel_viz, spixel_label_map = get_spixel_image((ori_img + mean_values).clamp(0, 1), ori_sz_spixel_map.squeeze(), n_spixels= int((H_ * scale * W_ * scale) / (args.downsize * args.downsize)))
        spixel_viz, spixel_label_map = get_spixel_image((ori_img + mean_values).clamp(0, 1), ori_sz_spixel_map.squeeze(), n_spixels=int(np.ceil(h)*np.ceil(w)))

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
