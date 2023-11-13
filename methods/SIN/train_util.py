import torch
import torch.nn.functional as F
import torch.nn as nn
import numpy as np
from skimage.segmentation import mark_boundaries
import cv2

import sys
sys.path.append('./third_party/cython')
# from connectivity import enforce_connectivity


def init_spixel_grid(args,  b_train=True):
    if b_train:
        img_height, img_width = args.train_img_height, args.train_img_width
    else:
        img_height, img_width = args.input_img_height, args.input_img_width

    # get spixel id for the final assignment
    n_spixl_h = int(np.floor(img_height/args.downsize))
    n_spixl_w = int(np.floor(img_width/args.downsize))

    spixel_height = int(img_height / (1. * n_spixl_h))
    spixel_width = int(img_width / (1. * n_spixl_w))

    spix_values = np.int32(np.arange(0, n_spixl_w * n_spixl_h).reshape((n_spixl_h, n_spixl_w)))
    spix_idx_tensor_ = shift9pos(spix_values)

    spix_idx_tensor =  np.repeat(
        np.repeat(spix_idx_tensor_, spixel_height,axis=1), spixel_width, axis=2)

    torch_spix_idx_tensor = torch.from_numpy(
                np.tile(spix_idx_tensor, (args.batch_size, 1, 1, 1))).type(torch.float).cuda()


    curr_img_height = int(np.floor(img_height))
    curr_img_width = int(np.floor(img_width))

    # pixel coord
    all_h_coords = np.arange(0, curr_img_height, 1)
    all_w_coords = np.arange(0, curr_img_width, 1)
    curr_pxl_coord = np.array(np.meshgrid(all_h_coords, all_w_coords, indexing='ij'))

    coord_tensor = np.concatenate([curr_pxl_coord[1:2, :, :], curr_pxl_coord[:1, :, :]])

    all_XY_feat = (torch.from_numpy(
        np.tile(coord_tensor, (args.batch_size, 1, 1, 1)).astype(np.float32)).cuda())

    return  torch_spix_idx_tensor, all_XY_feat

#===================== pooling and upsampling feature ==========================================

def shift9pos(input, h_shift_unit=1,  w_shift_unit=1):
    # input should be padding as (c, 1+ height+1, 1+width+1)
    input_pd = np.pad(input, ((h_shift_unit, h_shift_unit), (w_shift_unit, w_shift_unit)), mode='edge')
    input_pd = np.expand_dims(input_pd, axis=0)

    # assign to ...
    top     = input_pd[:, :-2 * h_shift_unit,          w_shift_unit:-w_shift_unit]
    bottom  = input_pd[:, 2 * h_shift_unit:,           w_shift_unit:-w_shift_unit]
    left    = input_pd[:, h_shift_unit:-h_shift_unit,  :-2 * w_shift_unit]
    right   = input_pd[:, h_shift_unit:-h_shift_unit,  2 * w_shift_unit:]

    center = input_pd[:,h_shift_unit:-h_shift_unit,w_shift_unit:-w_shift_unit]

    bottom_right    = input_pd[:, 2 * h_shift_unit:,   2 * w_shift_unit:]
    bottom_left     = input_pd[:, 2 * h_shift_unit:,   :-2 * w_shift_unit]
    top_right       = input_pd[:, :-2 * h_shift_unit,  2 * w_shift_unit:]
    top_left        = input_pd[:, :-2 * h_shift_unit,  :-2 * w_shift_unit]

    shift_tensor = np.concatenate([     top_left,    top,      top_right,
                                        left,        center,      right,
                                        bottom_left, bottom,    bottom_right], axis=0)
    return shift_tensor


def poolfeat(input, prob, sp_h=2, sp_w=2):

    def feat_prob_sum(feat_sum, prob_sum, shift_feat):
        feat_sum += shift_feat[:, :-1, :, :]
        prob_sum += shift_feat[:, -1:, :, :]
        return feat_sum, prob_sum

    b, _, h, w = input.shape

    h_shift_unit = 1
    w_shift_unit = 1
    p2d = (w_shift_unit, w_shift_unit, h_shift_unit, h_shift_unit)
    feat_ = torch.cat([input, torch.ones([b, 1, h, w]).cuda()], dim=1)  # b* (n+1) *h*w
    prob_feat = F.avg_pool2d(feat_ * prob.narrow(1, 0, 1), kernel_size=(sp_h, sp_w), stride=(sp_h, sp_w)) # b * (n+1) * h* w
    send_to_top_left = F.pad(prob_feat, p2d, mode='constant', value=0)[:, :, 2 * h_shift_unit:, 2 * w_shift_unit:]
    feat_sum = send_to_top_left[:, :-1, :, :].clone()
    prob_sum = send_to_top_left[:, -1:, :, :].clone()

    prob_feat = F.avg_pool2d(feat_ * prob.narrow(1, 1, 1), kernel_size=(sp_h, sp_w), stride=(sp_h, sp_w))  # b * (n+1) * h* w
    top = F.pad(prob_feat, p2d, mode='constant', value=0)[:, :, 2 * h_shift_unit:, w_shift_unit:-w_shift_unit]
    feat_sum, prob_sum = feat_prob_sum(feat_sum,prob_sum,top )

    prob_feat = F.avg_pool2d(feat_ * prob.narrow(1, 2, 1), kernel_size=(sp_h, sp_w), stride=(sp_h, sp_w))  # b * (n+1) * h* w
    top_right = F.pad(prob_feat, p2d, mode='constant', value=0)[:, :, 2 * h_shift_unit:, :-2 * w_shift_unit]
    feat_sum, prob_sum = feat_prob_sum(feat_sum, prob_sum, top_right)

    prob_feat = F.avg_pool2d(feat_ * prob.narrow(1, 3, 1), kernel_size=(sp_h, sp_w), stride=(sp_h, sp_w))  # b * (n+1) * h* w
    left = F.pad(prob_feat, p2d, mode='constant', value=0)[:, :, h_shift_unit:-h_shift_unit, 2 * w_shift_unit:]
    feat_sum, prob_sum = feat_prob_sum(feat_sum, prob_sum, left)

    prob_feat = F.avg_pool2d(feat_ * prob.narrow(1, 4, 1), kernel_size=(sp_h, sp_w), stride=(sp_h, sp_w))  # b * (n+1) * h* w
    center = F.pad(prob_feat, p2d, mode='constant', value=0)[:, :, h_shift_unit:-h_shift_unit, w_shift_unit:-w_shift_unit]
    feat_sum, prob_sum = feat_prob_sum(feat_sum, prob_sum, center)

    prob_feat = F.avg_pool2d(feat_ * prob.narrow(1, 5, 1), kernel_size=(sp_h, sp_w), stride=(sp_h, sp_w))  # b * (n+1) * h* w
    right = F.pad(prob_feat, p2d, mode='constant', value=0)[:, :, h_shift_unit:-h_shift_unit, :-2 * w_shift_unit]
    feat_sum, prob_sum = feat_prob_sum(feat_sum, prob_sum, right)

    prob_feat = F.avg_pool2d(feat_ * prob.narrow(1, 6, 1), kernel_size=(sp_h, sp_w), stride=(sp_h, sp_w))  # b * (n+1) * h* w
    bottom_left = F.pad(prob_feat, p2d, mode='constant', value=0)[:, :, :-2 * h_shift_unit, 2 * w_shift_unit:]
    feat_sum, prob_sum = feat_prob_sum(feat_sum, prob_sum, bottom_left)

    prob_feat = F.avg_pool2d(feat_ * prob.narrow(1, 7, 1), kernel_size=(sp_h, sp_w), stride=(sp_h, sp_w))  # b * (n+1) * h* w
    bottom = F.pad(prob_feat, p2d, mode='constant', value=0)[:, :, :-2 * h_shift_unit, w_shift_unit:-w_shift_unit]
    feat_sum, prob_sum = feat_prob_sum(feat_sum, prob_sum, bottom)

    prob_feat = F.avg_pool2d(feat_ * prob.narrow(1, 8, 1), kernel_size=(sp_h, sp_w), stride=(sp_h, sp_w))  # b * (n+1) * h* w
    bottom_right = F.pad(prob_feat, p2d, mode='constant', value=0)[:, :, :-2 * h_shift_unit, :-2 * w_shift_unit]
    feat_sum, prob_sum = feat_prob_sum(feat_sum, prob_sum, bottom_right)

    pooled_feat = feat_sum / (prob_sum + 1e-8)

    return pooled_feat


def  upfeat(input, prob, up_h=2, up_w=2):
    # input b*n*H*W  downsampled
    # prob b*9*h*w
    b, c, h, w = input.shape

    h_shift = 1
    w_shift = 1

    p2d = (w_shift, w_shift, h_shift, h_shift)
    feat_pd = F.pad(input, p2d, mode='constant', value=0)

    gt_frm_top_left = F.interpolate(feat_pd[:, :, :-2 * h_shift, :-2 * w_shift], size=(h * up_h, w * up_w),mode='nearest')
    feat_sum = gt_frm_top_left * prob.narrow(1,0,1)

    top = F.interpolate(feat_pd[:, :, :-2 * h_shift, w_shift:-w_shift], size=(h * up_h, w * up_w), mode='nearest')
    feat_sum += top * prob.narrow(1, 1, 1)

    top_right = F.interpolate(feat_pd[:, :, :-2 * h_shift, 2 * w_shift:], size=(h * up_h, w * up_w), mode='nearest')
    feat_sum += top_right * prob.narrow(1,2,1)

    left = F.interpolate(feat_pd[:, :, h_shift:-w_shift, :-2 * w_shift], size=(h * up_h, w * up_w), mode='nearest')
    feat_sum += left * prob.narrow(1, 3, 1)

    center = F.interpolate(input, (h * up_h, w * up_w), mode='nearest')
    feat_sum += center * prob.narrow(1, 4, 1)

    right = F.interpolate(feat_pd[:, :, h_shift:-w_shift, 2 * w_shift:], size=(h * up_h, w * up_w), mode='nearest')
    feat_sum += right * prob.narrow(1, 5, 1)

    bottom_left = F.interpolate(feat_pd[:, :, 2 * h_shift:, :-2 * w_shift], size=(h * up_h, w * up_w), mode='nearest')
    feat_sum += bottom_left * prob.narrow(1, 6, 1)

    bottom = F.interpolate(feat_pd[:, :, 2 * h_shift:, w_shift:-w_shift], size=(h * up_h, w * up_w), mode='nearest')
    feat_sum += bottom * prob.narrow(1, 7, 1)

    bottom_right =  F.interpolate(feat_pd[:, :, 2 * h_shift:, 2 * w_shift:], size=(h * up_h, w * up_w), mode='nearest')
    feat_sum += bottom_right * prob.narrow(1, 8, 1)

    return feat_sum


# ================= - spixel related -=============
def assign2uint8(assign):
    #red up, green mid, blue down, for debug only
    b,c,h,w = assign.shape

    red = torch.cat([torch.ones(size=assign.shape),  torch.zeros(size=[b,2,h,w])],dim=1).cuda()

    green = torch.cat([ torch.zeros(size=[b,1,h,w]),
                      torch.ones(size=assign.shape),
                      torch.zeros(size=[b,1,h,w])],dim=1).cuda()

    blue  = torch.cat([torch.zeros(size=[b,2,h,w]),
                       torch.ones(size=assign.shape)],dim=1).cuda()

    black = torch.zeros(size=[b,3,h,w]).cuda()
    white = torch.ones(size=[b,3,h,w]).cuda()
    # up probablity
    mat_vis = torch.where(assign.type(torch.float) < 0. , white, black)
    mat_vis = torch.where(assign.type(torch.float) >= 0. , red* (assign.type(torch.float)+1)/3, mat_vis)
    mat_vis = torch.where(assign.type(torch.float) >= 3., green*(assign.type(torch.float)-2)/3, mat_vis)
    mat_vis = torch.where(assign.type(torch.float) >= 6., blue * (assign.type(torch.float) - 5.) / 3, mat_vis)

    return (mat_vis * 255.).type(torch.uint8)

def val2uint8(mat,maxVal):
    maxVal_mat = torch.ones(mat.shape).cuda() * maxVal
    mat_vis = torch.where(mat > maxVal_mat, maxVal_mat, mat)
    return (mat_vis * 255. / maxVal).type(torch.uint8)


def update_spixl_map (spixl_map_idx_in, assig_map_in):
    assig_map = assig_map_in.clone()

    b,_,h,w = assig_map.shape
    _, _, id_h, id_w = spixl_map_idx_in.shape

    if (id_h == h) and (id_w == w):
        spixl_map_idx = spixl_map_idx_in
    else:
        spixl_map_idx = F.interpolate(spixl_map_idx_in, size=(h,w), mode='nearest')

    assig_max,_ = torch.max(assig_map, dim=1, keepdim= True)
    assignment_ = torch.where(assig_map == assig_max, torch.ones(assig_map.shape).cuda(),torch.zeros(assig_map.shape).cuda())
    new_spixl_map_ = spixl_map_idx * assignment_ # winner take all
    new_spixl_map = torch.sum(new_spixl_map_,dim=1,keepdim=True).type(torch.int)

    return new_spixl_map


def get_spixel_image(given_img, spix_index, n_spixels = 600, b_enforce_connect = False):

    if not isinstance(given_img, np.ndarray):
        given_img_np_ = given_img.detach().cpu().numpy().transpose(1,2,0)
    else: # for cvt lab to rgb case
        given_img_np_ = given_img

    if not isinstance(spix_index, np.ndarray):
        spix_index_np = spix_index.detach().cpu().numpy().transpose(0,1)
    else:
        spix_index_np = spix_index


    h, w = spix_index_np.shape
    given_img_np = cv2.resize(given_img_np_, dsize=(w, h), interpolation=cv2.INTER_CUBIC)

    if b_enforce_connect:
        spix_index_np = spix_index_np.astype(np.int64)
        segment_size = (given_img_np_.shape[0] * given_img_np_.shape[1]) / (int(n_spixels) * 1.0)
        min_size = int(0.06 * segment_size)
        max_size =  int(3 * segment_size)
        spix_index_np = enforce_connectivity(spix_index_np[None, :, :], min_size, max_size)[0]
    cur_max = np.max(given_img_np)
    spixel_bd_image = mark_boundaries(given_img_np/cur_max, spix_index_np.astype(int), color = (0,1,1)) #cyna
    return (cur_max*spixel_bd_image).astype(np.float32).transpose(2,0,1), spix_index_np #

# ============ accumulate Q =============================
def spixlIdx(args, b_train = False):
    # code modified from ssn
    if b_train:
        n_spixl_h = int(np.floor(args.train_img_height / args.downsize))
        n_spixl_w = int(np.floor(args.train_img_width / args.downsize))
    else:
        n_spixl_h = int(np.floor(args.input_img_height / args.downsize))
        n_spixl_w = int(np.floor(args.input_img_width / args.downsize))

    spix_values = np.int32(np.arange(0, n_spixl_w * n_spixl_h).reshape((n_spixl_h, n_spixl_w)))
    spix_idx_tensor = shift9pos(spix_values)

    torch_spix_idx_tensor = torch.from_numpy(
        np.tile(spix_idx_tensor, (args.batch_size, 1, 1, 1))).type(torch.float).cuda()

    return torch_spix_idx_tensor

class AverageMeter(object):
    """Computes and stores the average and current value"""

    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count

    def __repr__(self):
        return '{:.3f} ({:.3f})'.format(self.val, self.avg)

def batch2img(img):
    b,_,h,w = img.shape
    tmp = img.permute(0,2,3,1)
    for i in range(b):
        if i ==0:
            tmp_stack = tmp[i,:,:,:]
        else:
            tmp_stack = torch.cat([tmp_stack,tmp[i,:,:,:]],dim=-2)
    return  tmp_stack


def build_LABXY_feat(label_in, XY_feat):

    img_lab = label_in.clone().type(torch.float)

    b, _, curr_img_height, curr_img_width = XY_feat.shape
    scale_img =  F.interpolate(img_lab, size=(curr_img_height,curr_img_width), mode='nearest')
    LABXY_feat = torch.cat([scale_img, XY_feat],dim=1)

    return LABXY_feat


def rgb2Lab_torch(img_in, mean_values = None):
    # self implemented function that convert RGB image to LAB
    # inpu img intense should be [0,1] float b*3*h*w

    img= (img_in.clone() + mean_values.cuda()).clamp(0, 1)

    assert img.min() >= 0 and img.max() <= 1

    mask = img > 0.04045
    img[mask] = torch.pow((img[mask] + 0.055) / 1.055, 2.4)
    img[~mask] /= 12.92

    xyz_from_rgb = torch.tensor([[0.412453, 0.357580, 0.180423],
                             [0.212671, 0.715160, 0.072169],
                             [0.019334, 0.119193, 0.950227]]).cuda()
    rgb = img.permute(0,2,3,1)

    xyz_img = torch.matmul(rgb, xyz_from_rgb.transpose_(0,1))


    xyz_ref_white = torch.tensor([0.95047, 1., 1.08883]).cuda()

    # scale by CIE XYZ tristimulus values of the reference white point
    lab = xyz_img / xyz_ref_white

    # Nonlinear distortion and linear transformation
    mask = lab > 0.008856
    lab[mask] = torch.pow(lab[mask], 1. / 3.)
    lab[~mask] = 7.787 * lab[~mask] + 16. / 116.

    x, y, z = lab[..., 0:1], lab[..., 1:2], lab[..., 2:3]

    # Vector scaling
    L = (116. * y) - 16.
    a = 500.0 * (x - y)
    b = 200.0 * (y - z)

    return torch.cat([L, a, b], dim=-1).permute(0,3,1,2)


def label2one_hot_torch(labels, C=14):
    # w.r.t http://jacobkimmel.github.io/pytorch_onehot/
    '''
        Converts an integer label torch.autograd.Variable to a one-hot Variable.

        Parameters
        ----------
        labels : torch.autograd.Variable of torch.cuda.LongTensor
            N x 1 x H x W, where N is batch size.
            Each value is an integer representing correct classification.
        C : integer.
            number of classes in labels.

        Returns
        -------
        target : torch.cuda.FloatTensor
            N x C x H x W, where C is class number. One-hot encoded.
        '''
    b,_, h, w = labels.shape
    one_hot = torch.zeros(b, C, h, w, dtype=torch.long).cuda()
    target = one_hot.scatter_(1, labels.type(torch.long).data, 1) #require long type

    return target.type(torch.float32)


# ===============compute labxy loss(unsupervised)=======================
def get_lab_loss(img, img_grad, map_bound, map_non_bound):
    # _, _, h, w = img_grad.shape
    # img_grad = img_grad*img_grad
    # max_grad = img_grad.max(dim=-1)[0].max(dim=-1)[0]
    # mean_grad = torch.sum(torch.sum(img_grad, dim=-1), dim=-1)/(h*w)
    # img_grad = F.relu(img_grad-mean_grad.unsqueeze(-1).unsqueeze(-1))
    # overlap = torch.sum(torch.sum(map_bound*img_grad, dim=-1), dim=-1)
    # gt_bound_all = torch.sum(torch.sum(img_grad, dim=-1), dim=-1)
    # pred_bound_all = torch.sum(torch.sum(map_bound, dim=-1), dim=-1)
    # lab_loss = torch.sum(gt_bound_all-overlap, dim=-1)

    test_loss = torch.sum(torch.sum(img*img*map_non_bound, dim=-1), dim=-1)
    return test_loss

def gaussian_kernel(img):
    b, c, h, w = img.shape
    kernel = [[2, 4, 5, 4, 2],
              [4, 9, 12, 9, 4],
              [5, 12, 15, 12, 5],
              [4, 9, 12, 9, 4],
              [2, 4, 5, 4, 2]]
    kernel = torch.tensor(kernel).float().cuda().repeat(c, 1, 1).unsqueeze(1)*(1/139.)
    out = F.conv2d(img, kernel, padding=0, groups=c)
    out = F.pad(out, (2, 2, 2, 2), mode='replicate')
    return out



def labxy_loss(map, img):
    device = map.device
    bz, c, h, w = img.shape
    c = 1
    y_feat = torch.arange(0, h).repeat(bz, h, 1).unsqueeze(1).float().to(device)
    x_feat = y_feat.transpose(2, 3)
    kernel_img = [[[-1, 0, 0], [0, 1, 0], [0, 0, 0]],
                  [[0, -1, 0], [0, 1, 0], [0, 0, 0]],
                  [[0, 0, -1], [0, 1, 0], [0, 0, 0]],
                  [[0, 0, 0], [-1, 1, 0], [0, 0, 0]],
                  [[0, 0, 0], [0, 1, -1], [0, 0, 0]],
                  [[0, 0, 0], [0, 1, 0], [-1, 0, 0]],
                  [[0, 0, 0], [0, 1, 0], [0, -1, 0]],
                  [[0, 0, 0], [0, 1, 0], [0, 0, -1]]]
    # kernel_xy = [[[-1, 0, 0], [0, 1, 0], [0, 0, 0]],
    #               [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
    #               [[0, 0, -1], [0, 1, 0], [0, 0, 0]],
    #               [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
    #               [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
    #               [[0, 0, 0], [0, 1, 0], [-1, 0, 0]],
    #               [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
    #               [[0, 0, 0], [0, 1, 0], [0, 0, -1]]]
    # kernel_img = [[[0, -1, 0], [0, 1, 0], [0, 0, 0]],
    #               [[0, 0, 0], [-1, 1, 0], [0, 0, 0]],
    #               [[0, 0, 0], [0, 1, -1], [0, 0, 0]],
    #               [[0, 0, 0], [0, 1, 0], [0, -1, 0]]]
    # kernel_map = [[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]]
    kernel_img = torch.tensor(kernel_img).float().to(device).repeat(c+3, 1, 1).unsqueeze(1)
    # kernel_xy = torch.tensor(kernel_xy).float().to(device).repeat(2, 1, 1).unsqueeze(1)
    # kernel_map = torch.tensor(kernel_map).float().to(device).repeat(1, 1, 1, 1)
    # lab_img = rgb2Lab_torch(img, torch.tensor([0.411, 0.432, 0.45]).unsqueeze(-1).unsqueeze(-1))
    img = gaussian_kernel(img)
    gray_img = (0.2989*img[:, 0, :, :] + 0.5870*img[:, 1, :, :] + 0.1140*img[:, 2, :, :]).unsqueeze(1)
    cat_feat = torch.cat((gray_img, x_feat, y_feat, map.float()), dim=1)
    cat_feat = F.pad(cat_feat, (1, 1, 1, 1), mode='replicate')
    feat = F.conv2d(cat_feat, kernel_img, groups=c+3)
    cat_xy = F.pad(torch.cat((x_feat, y_feat), dim=1), (1, 1, 1, 1), mode='replicate')
    # xy_feat = F.conv2d(cat_xy, kernel_xy, groups=2)
    img_grad = feat[:, :8*c, :, :]
    xy_grad = feat[:, 8*c:-8, :, :]
    # xy_grad = xy_feat
    map_bound = feat[:, -8:, :, :]
    map_bound_8dim = -F.relu(-map_bound * map_bound+1)+1
    map_bound_1dim = -F.relu(-torch.sum(map_bound*map_bound, dim=1)+1)+1
    map_non_bound_1dim = -map_bound_1dim + 1
    map_non_bound_8dim = -map_bound_8dim + 1
    # img_grad_non_bound = torch.sum(img_grad*img_grad*map_non_bound_8dim.repeat(1, c, 1, 1), dim=1)
    # img_grad_all = torch.sum(img_grad*img_grad, dim=1)
    # img_grad *= map_non_bound_1dim
    # lab_loss = torch.sum(torch.sum(img_grad_non_bound, dim=-1), dim=-1)/torch.sum(torch.sum(img_grad_all, dim=-1), dim=-1)
    # lab_loss = torch.sum(torch.sum(img_grad_non_bound, dim=-1), dim=-1)
    # lab_loss = lab_loss/torch.sum(torch.sum(img_grad_all, dim=-1), dim=-1)
    lab_loss = get_lab_loss(gray_img.squeeze(1), img_grad, map_bound_1dim, map_non_bound_1dim)
    xy_grad = torch.sum(xy_grad * xy_grad * map_bound_8dim.repeat(1, 2, 1, 1), dim=1)
    xy_loss1 = torch.sum(torch.sum(xy_grad, dim=-1), dim=-1)
    xy_loss2 = torch.sum(torch.sum(map_bound_1dim, dim=-1), dim=-1)
    xy_loss = xy_loss1
    # map_bound = F.conv2d(map, kernel_map)
    return lab_loss, xy_loss


def labxy_v_loss(prob, label):
    _, color_c, _, _ = label.shape
    kernel = [[[0, -1, 0], [0, 1, 0], [0, 0, 0]],
              [[0, 0, 0], [0, 1, 0], [0, -1, 0]]]
    kernel = torch.tensor(kernel).float().cuda().repeat(color_c, 1, 1).unsqueeze(1)
    # cat_feat = img
    label = F.pad(label, (1, 1, 0, 0), mode='replicate')
    cat_feat = F.conv2d(label, kernel, stride=(2, 1), padding=(0, 0), groups=color_c)
    cat_feat = cat_feat*cat_feat

    cat_feat = F.relu(-F.relu(cat_feat)+1)
    b, c, h, w = cat_feat.shape
    _, gt_id = cat_feat.permute(0, 2, 3, 1).reshape(-1, 2).max(1, keepdim=False)

    # color_prob = cat_feat
    # if color_prob.shape[1] > 2:
    #     color_prob = color_prob[:, 0:2, :, :] + color_prob[:, 2:4, :, :] + color_prob[:, 4:, :, :]
    # color_prob = color_prob.permute(0, 2, 3, 1).reshape(-1, 2)
    # _, color_id = color_prob.min(1, keepdim=False)

    # b, _, h, w = gt.shape
    # gt_id = torch.where(gt > 0, torch.ones(gt.shape).cuda(), torch.zeros(gt.shape).cuda())
    # gt_id = gt_id.reshape(-1, ).long()
    # gt = F.softmax(gt, dim=1)
    # _, gt_id = gt.permute(0, 2, 3, 1).reshape(-1, 2).min(1, keepdim=False)

    cross_loss = nn.CrossEntropyLoss(reduction='none')
    color_loss = cross_loss(prob[:, :, 1:-1:2, :].permute(0, 2, 3, 1).reshape(-1, 2), gt_id)
    color_loss = color_loss.view(b, h, w)

    # gt_prob = F.softmax(cat_feat, dim=1)
    # pred_prob = prob[:, :, 1:-1:2, :]
    # loss = (gt_prob[:, 0, :, :]-pred_prob[:, 0, :, :])*(gt_prob[:, 0, :, :]-pred_prob[:, 0, :, :]) +\
    #        (gt_prob[:, 1, :, :] - pred_prob[:, 1, :, :]) * (gt_prob[:, 1, :, :] - pred_prob[:, 1, :, :])
    # loss = torch.sum(torch.sum(loss, dim=-1), dim=-1)

    gt = cat_feat[:, 0, :, :] - cat_feat[:, 1, :, :]
    weight = gt * gt
    color_loss = weight*color_loss
    color_loss = torch.sum(torch.sum(color_loss, dim=-1), dim=-1)

    # regular_weight = 1-weight
    # regular_loss = regular_weight * (prob[:, 0, 1:-1:2, :]-prob[:, 1, 1:-1:2, :])
    # regular_loss = torch.sum(torch.sum(regular_loss, dim=-1), dim=-1)
    # regular_loss = regular_loss * regular_loss
    # color_loss += regular_loss

    return color_loss


def labxy_h_loss(prob, label):
    _, color_c, _, _ = label.shape
    kernel = [[[0, 0, 0], [-1, 1, 0], [0, 0, 0]],
              [[0, 0, 0], [0, 1, -1], [0, 0, 0]]]
    kernel = torch.tensor(kernel).float().cuda().repeat(color_c, 1, 1).unsqueeze(1)
    # cat_feat = img
    label = F.pad(label, (0, 0, 1, 1), mode='replicate')
    cat_feat = F.conv2d(label, kernel, stride=(1, 2), padding=(0, 0), groups=color_c)
    cat_feat = cat_feat * cat_feat

    # b, c, h, w = cat_feat.shape
    # gt_prob = F.softmax(cat_feat.view(-1, 2, h, w), dim=1).view(b, -1, h, w)
    #
    # color_prob = cat_feat
    # if color_prob.shape[1] > 2:
    #     color_prob = color_prob[:, 0:2, :, :] + color_prob[:, 2:4, :, :] + color_prob[:, 4:, :, :]
    # color_prob = color_prob.permute(0, 2, 3, 1).reshape(-1, 2)
    # _, color_id = color_prob.min(1, keepdim=False)

    # b, _, h, w = gt.shape
    # gt_id = torch.where(gt > 0, torch.ones(gt.shape).cuda(), torch.zeros(gt.shape).cuda())
    # gt_id = gt_id.reshape(-1, ).long()
    # gt = F.softmax(gt, dim=1)
    # _, gt_id = gt.permute(0, 2, 3, 1).reshape(-1, 2).min(1, keepdim=False)

    cat_feat = F.relu(-F.relu(cat_feat) + 1)
    b, c, h, w = cat_feat.shape
    _, gt_id = cat_feat.permute(0, 2, 3, 1).reshape(-1, 2).max(1, keepdim=False)

    cross_loss = nn.CrossEntropyLoss(reduction='none')
    color_loss = cross_loss(prob[:, :, :, 1:-1:2].permute(0, 2, 3, 1).reshape(-1, 2), gt_id)
    color_loss = color_loss.view(b, h, w)

    gt = cat_feat[:, 0, :, :] - cat_feat[:, 1, :, :]
    weight = gt * gt
    color_loss = weight * color_loss
    color_loss = torch.sum(torch.sum(color_loss, dim=-1), dim=-1)

    # regular_weight = 1 - weight
    # regular_loss = regular_weight * (prob[:, 0, :, 1:-1:2] - prob[:, 1, :, 1:-1:2])
    # regular_loss = torch.sum(torch.sum(regular_loss, dim=-1), dim=-1)
    # regular_loss = regular_loss*regular_loss
    # color_loss += regular_loss

    # gt_prob = F.softmax(cat_feat, dim=1)
    # pred_prob = prob[:, :, :, 1:-1:2]
    # loss = (gt_prob[:, 0, :, :] - pred_prob[:, 0, :, :]) * (gt_prob[:, 0, :, :] - pred_prob[:, 0, :, :]) + \
    #        (gt_prob[:, 1, :, :] - pred_prob[:, 1, :, :]) * (gt_prob[:, 1, :, :] - pred_prob[:, 1, :, :])
    # loss = torch.sum(torch.sum(loss, dim=-1), dim=-1)
    return color_loss


def compute_lr_grad(img):
    img_expand = F.pad(img, (1, 1, 0, 0), mode='replicate')
    img_l = img_expand[:, :, :, :-2]
    img_r = img_expand[:, :, :, 2:]
    l_grad = img - img_l
    r_grad = img_r - img
    lr_grad = l_grad*l_grad + r_grad*r_grad
    lr_grad = torch.sum(lr_grad, dim=1, keepdim=True)
    return lr_grad


def compute_tb_grad(img):
    img_expand = F.pad(img, (0, 0, 1, 1), mode='replicate')
    img_t = img_expand[:, :, :-2, :]
    img_b = img_expand[:, :, 2:, :]
    t_grad = img - img_t
    b_grad = img_b - img
    tb_grad = t_grad*t_grad + b_grad*b_grad
    tb_grad = torch.sum(tb_grad, dim=1, keepdim=True)
    return tb_grad


def compute_gt(img):
    _, c, _, _ = img.shape
    kernel_17 = [[1, 1, 1, 1, 1, 1, 1, 1, -8, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, -8, 1, 1, 1, 1, 1, 1, 1, 1]]
    kernel_9 = [[1, 1, 1, 1, -4, 0, 0, 0, 0],
                [0, 0, 0, 0, -4, 1, 1, 1, 1]]
    kernel_5 = [[1, 1, -2, 0, 0],
                [0, 0, -2, 1, 1]]
    kernel_3 = [[1, -1, 0],
                [0, -1, 1]]
    kernel_17 = torch.tensor(kernel_17).float().cuda().repeat(c, 1).unsqueeze(1).unsqueeze(1)
    kernel_9 = torch.tensor(kernel_9).float().cuda().repeat(c, 1).unsqueeze(1).unsqueeze(1)
    kernel_5 = torch.tensor(kernel_5).float().cuda().repeat(c, 1).unsqueeze(1).unsqueeze(1)
    kernel_3 = torch.tensor(kernel_3).float().cuda().repeat(c, 1).unsqueeze(1).unsqueeze(1)
    gt_h_17 = F.conv2d(img, kernel_17, stride=16, groups=c)
    gt_h_9 = F.conv2d(img, kernel_9, stride=8, groups=c)
    gt_h_5 = F.conv2d(img, kernel_5, stride=4, groups=c)
    gt_h_3 = F.conv2d(img, kernel_3, stride=2, groups=c)
    gt_v_17 = F.conv2d(img, kernel_17.transpose(2, 3), stride=(16, 8), groups=c)
    gt_v_9 = F.conv2d(img, kernel_9.transpose(2, 3), stride=(8, 4), groups=c)
    gt_v_5 = F.conv2d(img, kernel_5.transpose(2, 3), stride=(4, 2), groups=c)
    gt_v_3 = F.conv2d(img, kernel_3.transpose(2, 3), stride=(2, 1), groups=c)

    gt_h_17 = gt_h_17 * gt_h_17
    gt_h_17_l = torch.sum(gt_h_17[:, 0::2, :, :], dim=1, keepdim=True)
    gt_h_17_r = torch.sum(gt_h_17[:, 1::2, :, :], dim=1, keepdim=True)
    gt_h_17 = torch.cat((gt_h_17_l, gt_h_17_r), dim=1)

    gt_h_9 = gt_h_9 * gt_h_9
    gt_h_9_l = torch.sum(gt_h_9[:, 0::2, :, :], dim=1, keepdim=True)
    gt_h_9_r = torch.sum(gt_h_9[:, 1::2, :, :], dim=1, keepdim=True)
    gt_h_9 = torch.cat((gt_h_9_l, gt_h_9_r), dim=1)

    gt_h_5 = gt_h_5 * gt_h_5
    gt_h_5_l = torch.sum(gt_h_5[:, 0::2, :, :], dim=1, keepdim=True)
    gt_h_5_r = torch.sum(gt_h_5[:, 1::2, :, :], dim=1, keepdim=True)
    gt_h_5 = torch.cat((gt_h_5_l, gt_h_5_r), dim=1)

    gt_h_3 = gt_h_3 * gt_h_3
    gt_h_3_l = torch.sum(gt_h_3[:, 0::2, :, :], dim=1, keepdim=True)
    gt_h_3_r = torch.sum(gt_h_3[:, 1::2, :, :], dim=1, keepdim=True)
    gt_h_3 = torch.cat((gt_h_3_l, gt_h_3_r), dim=1)

    gt_v_17 = gt_v_17 * gt_v_17
    gt_v_17_l = torch.sum(gt_v_17[:, 0::2, :, :], dim=1, keepdim=True)
    gt_v_17_r = torch.sum(gt_v_17[:, 1::2, :, :], dim=1, keepdim=True)
    gt_v_17 = torch.cat((gt_v_17_l, gt_v_17_r), dim=1)

    gt_v_9 = gt_v_9 * gt_v_9
    gt_v_9_l = torch.sum(gt_v_9[:, 0::2, :, :], dim=1, keepdim=True)
    gt_v_9_r = torch.sum(gt_v_9[:, 1::2, :, :], dim=1, keepdim=True)
    gt_v_9 = torch.cat((gt_v_9_l, gt_v_9_r), dim=1)

    gt_v_5 = gt_v_5 * gt_v_5
    gt_v_5_l = torch.sum(gt_v_5[:, 0::2, :, :], dim=1, keepdim=True)
    gt_v_5_r = torch.sum(gt_v_5[:, 1::2, :, :], dim=1, keepdim=True)
    gt_v_5 = torch.cat((gt_v_5_l, gt_v_5_r), dim=1)

    gt_v_3 = gt_v_3 * gt_v_3
    gt_v_3_l = torch.sum(gt_v_3[:, 0::2, :, :], dim=1, keepdim=True)
    gt_v_3_r = torch.sum(gt_v_3[:, 1::2, :, :], dim=1, keepdim=True)
    gt_v_3 = torch.cat((gt_v_3_l, gt_v_3_r), dim=1)

    return gt_h_17, gt_h_9, gt_h_5, gt_h_3, gt_v_17, gt_v_9, gt_v_5, gt_v_3





if __name__ == '__main__':
    import os
    from torchvision import transforms
    w = 255
    s_w = 15
    map_in = torch.arange(1, s_w*s_w+1).reshape(s_w, s_w).repeat(1, 1, 1, 1).cuda()
    map_in = F.interpolate(map_in.float(), (w, w), mode='nearest')
    output = './demo/lab_loss'
    if not os.path.isdir(output):
        os.makedirs(output)
    img = cv2.imread('./demo/inputs/birds.jpg')[:225, :225, :]
    img = torch.tensor(img).float().cuda().permute(2, 0, 1).contiguous()
    norm = transforms.Normalize([0, 0, 0], [255, 255, 255])
    img = norm(img).unsqueeze(0)
    # img = img.unsqueeze(0)
    img0 = img[0]
    cv2.imwrite(os.path.join(output, 'norm.jpg'), img[0].detach().cpu().numpy().transpose(1,2,0))
    img = gaussian_kernel(img)
    cv2.imwrite(os.path.join(output, 'gaussian_smooth.jpg'), img[0].detach().cpu().numpy().transpose(1,2,0))
    lab_img = rgb2Lab_torch(img, torch.tensor([0.411,0.432,0.45]).unsqueeze(-1).unsqueeze(-1))
    img_lr_grad = compute_lr_grad(lab_img)
    img_tb_grad = compute_tb_grad(lab_img)
    gt_h_17, gt_h_9, gt_h_5, gt_h_3, gt_v_17, gt_v_9, gt_v_5, gt_v_3 = compute_gt(lab_img)
    p0v = torch.randn(1, 2, 225, 225).cuda()
    lab_v_loss = labxy_v_loss(p0v, gt_v_3)

