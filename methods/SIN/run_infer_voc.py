import argparse
import os
from glob import glob
import torch.backends.cudnn as cudnn
import models
import torchvision.transforms as transforms
import flow_transforms
# from scipy.ndimage import imread
# from scipy.misc import imsave, imresize
from imageio import imread, imsave
from skimage import img_as_ubyte
# from skimage.transform import resize as imresize
from loss import *
import time
import random
from models.model_util import update_spixel_map

# import sys
# sys.path.append('./third_party/cython')
# from connectivity import enforce_connectivity

'''
Infer from bsds500 dataset:
author:Fengting Yang 
last modification:  Mar.14th 2019

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

parser.add_argument('--data_dir', metavar='DIR', default='/media/yuanqing/ssd/code/visual_tracking/bilateralinceptions-master/data/voc_test/VOCdevkit/VOC2012/JPEGImages/',help='path to images folder')
# parser.add_argument('--data_dir', metavar='DIR', default='/media/yuanqing/ssd/code/visual_tracking/bilateralinceptions-master/data/VOCdevkit/VOC2012/JPEGImages/',help='path to images folder')
parser.add_argument('--pretrained', metavar='PTH', help='path to pre-trained model', default='./pretrain_ckpt/voc_model/model_best_20m.tar')
parser.add_argument('--output', metavar='DIR', default= './output/voc_test/finetune',help='path to output folder')

parser.add_argument('--downsize', default=16, type=float, help='superpixel grid cell, must be same as training setting')
parser.add_argument('-b', '--batch-size', default=1, type=int,  metavar='N', help='mini-batch size')

# the BSDS500 has two types of image, horizontal and veritical one, here I use train_img and input_img to presents them respectively
parser.add_argument('--train_img_height', '-tcimgH', default=320 ,  type=int, help='img height must be 16*n')
parser.add_argument('--train_img_width', '-t_imgW', default=480,  type=int, help='img width must be 16*n')
parser.add_argument('--input_img_height', '-v_imgH', default=480,   type=int, help='img height_must be 16*n')  #
parser.add_argument('--input_img_width', '-v_imgW', default=320,    type=int, help='img width must be 16*n')

args = parser.parse_args()
args.test_list = '/media/yuanqing/ssd/code/visual_tracking/bilateralinceptions-master/data/voc_test/VOCdevkit/VOC2012/' + 'test.txt'
# args.test_list = '/media/yuanqing/ssd/code/visual_tracking/bilateralinceptions-master/data/' + 'reducedval.txt'

random.seed(100)
@torch.no_grad()
def test(model, img_paths, save_path, idx):
    # Data loading code
    input_transform = transforms.Compose([
        flow_transforms.ArrayToTensor(),
        transforms.Normalize(mean=[0,0,0], std=[255,255,255]),
        transforms.Normalize(mean=[0.411,0.432,0.45], std=[1,1,1])
    ])

    imgId = img_paths[idx]
    img_paths = args.data_dir + imgId + '.jpg'
    img_file = img_paths
    load_path = img_file
    # imgId = os.path.basename(img_file)[:-4]

    # origin size 481*321 or 321*481
    img_ = imread(load_path)
    H, W = img_.shape[0:2]
    if img_.size == (H*W):
        img_ = cv2.cvtColor(img_, cv2.COLOR_GRAY2RGB)

    H_, W_  = int(np.ceil(H*1./16.)*16-15), int(np.ceil(W*1./16.)*16-15)

    img = cv2.resize(img_, (W_, H_), interpolation=cv2.INTER_CUBIC)
    # img1 = input_transform(img1)
    # ori_img = input_transform(img_)
    #
    # # choose the right spixelIndx
    # if H_ == 321 and W_==481:
    #     # spixl_map_idx_tensor = spixeIds[0]
    #     img = cv2.resize(img_, (int(480 * scale)-15, int(320 * scale)-15), interpolation=cv2.INTER_CUBIC)
    # elif H_ == 481 and W_ == 321:
    #     # spixl_map_idx_tensor = spixeIds[1]
    #     img = cv2.resize(img_, (int(320 * scale)-15, int(480 * scale)-15), interpolation=cv2.INTER_CUBIC)
    # else:
    #     print('The image size is wrong!')
    #     return

    img1 = input_transform(img)
    ori_img = input_transform(img_)
    mean_values = torch.tensor([0.411, 0.432, 0.45], dtype=img1.cuda().unsqueeze(0).dtype).view(3, 1, 1)

    # compute output
    tic = time.time()
    prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h = model(img1.cuda().unsqueeze(0))
    # assign the spixel map and  resize to the original size
    curr_spixl_map = update_spixel_map(img1.cuda().unsqueeze(0), prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h)
    # curr_spixl_map = map0
    ori_sz_spixel_map = F.interpolate(curr_spixl_map.type(torch.float), size=(H, W), mode='nearest').type(torch.int)

    spix_index_np = ori_sz_spixel_map.squeeze().detach().cpu().numpy().transpose(0, 1)
    spix_index_np = spix_index_np.astype(np.int64)
    # segment_size = (spix_index_np.shape[0] * spix_index_np.shape[1]) / (int( 600*scale*scale) * 1.0)
    # min_size = int(0.06 * segment_size)
    # max_size = int(3 * segment_size)
    # spixel_label_map = enforce_connectivity(spix_index_np[None, :, :], min_size, max_size)[0]
    spixel_label_map = spix_index_np

    torch.cuda.synchronize()
    toc = time.time() - tic

    n_spixel = len(np.unique(spixel_label_map))
    given_img_np = (ori_img + mean_values).clamp(0, 1).detach().cpu().numpy().transpose(1, 2, 0)
    spixel_bd_image = mark_boundaries(given_img_np / np.max(given_img_np), spixel_label_map.astype(int), color=(1, 1, 1))
    spixel_viz = spixel_bd_image.astype(np.float32).transpose(2, 0, 1)

    # ************************ Save all result********************************************
    # save img, uncomment it if needed
    # if not os.path.isdir(os.path.join(save_path, 'img')):
    #     os.makedirs(os.path.join(save_path, 'img'))
    # spixl_save_name = os.path.join(save_path, 'img', imgId + '.jpg')
    # img_save = (ori_img + mean_values).clamp(0, 1)
    # imsave(spixl_save_name, img_save.detach().cpu().numpy().transpose(1, 2, 0))

    # save spixel viz
    if not os.path.isdir(os.path.join(save_path, 'spixel_viz')):
        os.makedirs(os.path.join(save_path, 'spixel_viz'))
    spixl_save_name = os.path.join(save_path, 'spixel_viz', imgId + '_sPixel.png')
    imsave(spixl_save_name, img_as_ubyte(spixel_viz.transpose(1, 2, 0)))

    # save the unique maps as csv for eval
    if not os.path.isdir(os.path.join(save_path, 'map_csv')):
        os.makedirs(os.path.join(save_path, 'map_csv'))
    output_path = os.path.join(save_path, 'map_csv', imgId + '.csv')
    # plus 1 to make it consistent with the toolkit format
    np.savetxt(output_path, spixel_label_map .astype(int), fmt='%i', delimiter=",")


    if idx % 10 == 0:
        print("processing %d"%idx)

    return toc, n_spixel

def main():
    global args, save_path
    data_dir = args.data_dir
    print("=> fetching img pairs in '{}'".format(data_dir))

    train_img_height = args.train_img_height
    train_img_width = args.train_img_width
    input_img_height = args.input_img_height
    input_img_width = args.input_img_width


    save_path = args.output

    print('=> will save everything to {}'.format(save_path))
    if not os.path.isdir(save_path):
        os.makedirs(save_path)

    tst_lst = []
    with open(args.test_list, 'r') as tf:
        img_path = tf.readlines()
        for path in img_path:
            tst_lst.append(path[:-1])
    print('{} samples found'.format(len(tst_lst)))

    # create model
    network_data = torch.load(args.pretrained)
    print("=> using pre-trained model '{}'".format(network_data['arch']))
    model = models.__dict__[network_data['arch']]( data = network_data).cuda()
    model.eval()
    args.arch = network_data['arch']
    cudnn.benchmark = True

    # the following code is for debug
    for n in range(len(tst_lst)):
      time, n_spixel = test(model, tst_lst, save_path, n)

    print("generate {} spixels".format(n_spixel))


if __name__ == '__main__':
    main()
