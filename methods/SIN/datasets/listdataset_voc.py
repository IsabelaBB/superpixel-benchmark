import torch.utils.data as data
import os
import cv2
import numpy as np


class ListDataset(data.Dataset):
    def __init__(self, root, dataset, path_list, transform=None, target_transform=None,
                 co_transform=None, loader=None, datatype=None):

        self.root = root
        self.dataset = dataset
        self.img_path_list =path_list
        self.transform = transform
        self.target_transform = target_transform
        self.co_transform = co_transform
        self.loader = loader
        self.datatype = datatype

    def __getitem__(self, index):
        img_root = '/media/yuanqing/ssd/code/visual_tracking/bilateralinceptions-master/data/VOCdevkit/VOC2012/JPEGImages/'
        img_id = self.img_path_list[index][:-1]
        img_path = img_root + img_id + '.jpg'
        label_root = '/media/yuanqing/ssd/code/visual_tracking/bilateralinceptions-master/data/VOCdevkit/VOC2012/SegmentationClass/'
        label_path = label_root + img_id + '.png'
        # print('img_id:', img_id)
        # img_path = self.img_path_list[0][:-1]
        # We do not consider other datsets in this work
        # assert self.dataset == 'bsd500'
        assert (self.transform is not None) and (self.target_transform is not None)

        inputs, label = self.loader(img_path, label_path)
        # print('input1 shape:', inputs.shape)
        # if inputs.shape[0]<=193 or inputs.shape[1]<=193:
        #     inputs = cv2.resize(inputs, (255, 255), interpolation=cv2.INTER_CUBIC)
        #     label = cv2.resize(label, (255, 255), interpolation=cv2.INTER_CUBIC)
        #     label = label[:, :, np.newaxis]
        inputs = cv2.resize(inputs, (193, 193), interpolation=cv2.INTER_CUBIC)
        label = cv2.resize(label, (193, 193), interpolation=cv2.INTER_CUBIC)
        label = label[:, :, np.newaxis]

        if self.co_transform is not None:
            inputs, label = self.co_transform([inputs], label)

        # print('inputs2 shape:', inputs[0].shape)
        if self.transform is not None:
            image = self.transform(inputs[0])

        # print('label shape:', label.shape )
        if self.target_transform is not None:
            label = self.target_transform(label)

        return image, label

    def __len__(self):
        return len(self.img_path_list)
