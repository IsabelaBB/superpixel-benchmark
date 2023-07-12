import skimage
import argparse
import sys
import numpy as np
import cv2 as cv
import os
from glob import glob

## input arguments 
def parse_arguments(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help='input image file/path')
    parser.add_argument('--output', type=str, default=None, help='output image file/path')
    parser.add_argument('--noise', type=str, help='noise type (avgblur or salt_pepper)')
    parser.add_argument('--param', type=float, default=None, help='Noise probability or blur kernel')
    
    args = parser.parse_args(argv)
    return args

def addNoise(input, output, noise, param):
  if noise == 'avgblur':
    img = cv.imread(input)
    gimg = cv.blur(img,(int(param),int(param)))
    cv.imwrite(output, gimg)
  else:
    img = skimage.io.imread(input)/255.0
    gimg = skimage.util.random_noise(img, mode=noise, seed=42, amount=param)
    gimg = np.array(255*gimg, dtype = 'uint8')
    skimage.io.imsave(output, gimg)

def main():
  args = parse_arguments(sys.argv[1:])

  if(os.path.isdir(args.input)):
    if(os.path.isdir(args.output)== False):
      print("ERROR: When input is a path, output must also be a path.")
      return
    
    img_files = glob(args.input+'/*')
    img_files.sort()

    for n in range(len(img_files)):
      img_file = img_files[n]
      img_name = os.path.basename(img_file).split('.')[0]
      output_file = args.output + '/' + img_name + '.png'
      addNoise(img_file, output_file, args.noise, args.param)
  
  else:
    if(os.path.isdir(args.output)):
      print("ERROR: When input is a file, output must also be a file.")
      return

    addNoise(args.input, args.output, args.noise, args.param)


if __name__ == '__main__':
  sys.exit(main())



