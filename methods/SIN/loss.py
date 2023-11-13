import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from train_util import *

'''
Loss function
author:Fengting Yang 
Mar.1st 2019

We only use "compute_semantic_pos_loss" func. in our final version, best result achieved with weight = 3e-3
'''


def compute_semantic_pos_loss(prob_in, labxy_feat,  pos_weight = 0.003,  kernel_size=16):
    # this wrt the slic paper who used sqrt of (mse)

    # rgbxy1_feat: B*50+2*H*W
    # output : B*9*H*w
    # NOTE: this loss is only designed for one level structure

    # todo: currently we assume the downsize scale in x,y direction are always same
    S = kernel_size
    m = pos_weight
    prob = prob_in.clone()

    b, c, h, w = labxy_feat.shape
    pooled_labxy = poolfeat(labxy_feat, prob, kernel_size, kernel_size)
    reconstr_feat = upfeat(pooled_labxy, prob, kernel_size, kernel_size)

    loss_map = reconstr_feat[:,-2:,:,:] - labxy_feat[:,-2:,:,:]

    # self def cross entropy  -- the official one combined softmax
    logit = torch.log(reconstr_feat[:, :-2, :, :] + 1e-8)
    loss_sem = - torch.sum(logit * labxy_feat[:, :-2, :, :]) / b
    loss_pos = torch.norm(loss_map, p=2, dim=1).sum() / b * m / S

    # empirically we find timing 0.005 tend to better performance
    loss_sum = 0.005 * (loss_sem + loss_pos)
    loss_sem_sum = 0.005 * loss_sem
    loss_pos_sum = 0.005 * loss_pos

    return loss_sum, loss_sem_sum,  loss_pos_sum


def compute_labxy_loss(prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h, label):
    p0v = prob0_v.clone()
    p0h = prob0_h.clone()
    p1v = prob1_v.clone()
    p1h = prob1_h.clone()
    p2v = prob2_v.clone()
    p2h = prob2_h.clone()
    p3v = prob3_v.clone()
    p3h = prob3_h.clone()
    b, c, h, w = label.shape
    # weight = [0.25, 0.5, 1., 2., 4., 8., 16., 32.]
    weight = [1, 2.5, 2, 5, 4, 10, 8, 20]
    # weight = [1, 1, 1, 1, 1, 1, 1, 1]
    weight = torch.tensor(weight).cuda().reshape(8, ).float()

    # img = rgb2Lab_torch(img, torch.tensor([0.411, 0.432, 0.45]).unsqueeze(-1).unsqueeze(-1))
    # img = gaussian_kernel(img)
    # img_lr_grad = compute_lr_grad(img)
    # img_tb_grad = compute_tb_grad(img)
    # gt_h_17, gt_h_9, gt_h_5, gt_h_3, gt_v_17, gt_v_9, gt_v_5, gt_v_3 = compute_gt(img)
    # img = (0.2989 * img[:, 0, :, :] + 0.5870 * img[:, 1, :, :] + 0.1140 * img[:, 2, :, :]).unsqueeze(1)

    lab_loss = torch.zeros(8, b).cuda()

    # todo: complete labxy_loss
    label_0v = label
    lab_loss[0] = labxy_v_loss(p0v, label_0v)
    label_0h = label_0v[:, :, 0::2, :]
    lab_loss[1] = labxy_h_loss(p0h, label_0h)

    label_1v = label_0h[:, :, :, 0::2]
    lab_loss[2] = labxy_v_loss(p1v, label_1v)
    label_1h = label_1v[:, :, 0::2, :]
    lab_loss[3] = labxy_h_loss(p1h, label_1h)

    label_2v = label_1h[:, :, :, 0::2]
    lab_loss[4] = labxy_v_loss(p2v, label_2v)
    label_2h = label_2v[:, :, 0::2, :]
    lab_loss[5] = labxy_h_loss(p2h, label_2h)

    label_3v = label_2h[:, :, :, 0::2]
    lab_loss[6] = labxy_v_loss(p3v, label_3v)
    label_3h = label_3v[:, :, 0::2, :]
    lab_loss[7] = labxy_h_loss(p3h, label_3h)

    lab_loss = torch.sum(lab_loss, dim=-1) / b

    lab_loss = torch.sum(lab_loss * weight, dim=0)

    lab_loss = 0.005 * lab_loss

    return lab_loss
