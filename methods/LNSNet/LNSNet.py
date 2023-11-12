import torch

# from tensorboardX import SummaryWriter
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from skimage.segmentation import mark_boundaries
from glob import glob

from libs.model import *
from libs.test import *
from PIL import Image
from time import perf_counter
import random

random.seed(1)
torch.manual_seed(1)

parser = argparse.ArgumentParser()

###### Data Setting#######
parser.add_argument('--img', default='demo_imgs/012.jpg',
                    help='The path of the source image')
parser.add_argument('--label', default='result_012.jpg',
                    help='The path of the labeled image')
parser.add_argument('--border', default=None,
                    help='The path of the borders image')
parser.add_argument('--dtime', default=None, help='The path of the time log')
parser.add_argument('--time', default=None, help='The path of the time log')

###### Model Setting#######
parser.add_argument('--device', default="cuda" if torch.cuda.is_available() else "cpu", help='use cuda or cpu')
parser.add_argument('--k', type=int, default=100, help='Number of superpixels')
parser.add_argument('--kn', type=int, default=16,
                    help='The number of dis limit')
parser.add_argument('--seed_strategy', type=str,
                    default='network', help='network/grid')
##### Optimizing Setting####
parser.add_argument('--lr', type=float, default=3e-4, help='learning rate')
parser.add_argument('--use_gal', default=True, help='Using sowm weight update')
parser.add_argument('--use_gbl', default=True, help='Using sowm weight update')
parser.add_argument('--is_dilation', default=True,
                    help='Using dilation convolution')
parser.add_argument('--check_path', type=str,
                    default='./lnsnet_BSDS_checkpoint.pth')

################

args = parser.parse_args()

def runDirectory():
  if(os.path.isdir(args.img)):
    # read all the image files in the folder 
    img_files = glob(args.img+'/*')
    img_files.sort()
  else:
    img_files = [args.img]

  overall_time = 0.0
  model = LNSN(args.k, args)
  model.load_state_dict(torch.load(args.check_path))

  # for each image:
  for n in range(len(img_files)):
    img_file = img_files[n]
    img_name = os.path.basename(img_file).split('.')[0]
    image = plt.imread(img_file)
    
    if (len(image.shape) == 2):
      image = np.stack((image,)*3, axis=-1)
    
    input = preprocess(image, args.device)

    with torch.no_grad():
      b, _, h, w = input.size()
      if args.device == "cuda":
        model.cuda()
        input.cuda()

      start = perf_counter()

      recons, cx, cy, f, probs = model.forward(input, torch.zeros(h, w))
      labels = assignment_test(f, input, cx, cy)
      labels = labels.permute(0, 2, 1).contiguous().view(b, -1, h, w)
      labels = labels.argmax(1).squeeze().to("cpu").detach().numpy()

      end = perf_counter()
      t = end - start
      
      if image.shape[:2] != labels.shape[-2:]:
        labels = labels.transpose(1, 0)

    if(args.label is not None):
      label_image = Image.fromarray(labels.astype(np.uint16))
      spixl_save_name = args.label 
      if(os.path.isdir(spixl_save_name)):
        spixl_save_name = spixl_save_name + "/" + img_name + '.png'
      label_image.save(spixl_save_name)

    if (args.border is not None):
      boundaries = mark_boundaries(image, labels, color=(1, 0, 0))
      spixl_save_name = args.border
      if(os.path.isdir(spixl_save_name)):
        spixl_save_name = spixl_save_name + "/" + img_name + '.png'
      plt.imsave(spixl_save_name, boundaries)

    if (args.dtime is not None):
      with open(args.dtime, "a+") as fileTime:
        fileTime.write(img_name + " " + "%.5f" % t + "\n")
  
    overall_time += t

  if (args.time is not None):
    with open(args.time, "a+") as fileTime:
      fileTime.write(str(args.k) + " " + "%.5f" % (overall_time/len(img_files)) + "\n")


if __name__ == '__main__':
    runDirectory()
