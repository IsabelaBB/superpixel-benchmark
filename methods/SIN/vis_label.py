import os
import torchvision.transforms as transforms
import flow_transforms
from imageio import imread, imsave
from skimage import img_as_ubyte
from loss import *


def vis(csv_path, save_path):
    spixel_label_map = np.loadtxt(csv_path, delimiter=",")
    max_id = np.max(spixel_label_map)
    h, w = spixel_label_map.shape
    rand_c = torch.rand((int(max_id), 3))

    spixel_viz = torch.zeros((3, h, w))
    for i in range(h):
        for j in range(w):
            spixel_viz[:, i, j] = rand_c[int(spixel_label_map[i, j])-1]
    spixel_viz = np.array(spixel_viz)

    n_spixel = len(np.unique(spixel_label_map))

    imgId = os.path.basename(csv_path)[:-4]

    # ************************ Save all result*******************************************
    # save spixel viz
    if not os.path.isdir(os.path.join(save_path, 'gt_viz')):
        os.makedirs(os.path.join(save_path, 'gt_viz'))
    spixl_save_name = os.path.join(save_path, 'gt_viz', imgId + '_sPixel.png')
    # imsave(spixl_save_name, img_as_ubyte(spixel_viz))
    imsave(spixl_save_name, img_as_ubyte(spixel_viz.transpose(1, 2, 0)))

    return


def main():
    csv_path = './output/test_multiscale/SPixelNet_nSpixel_600/map_csv/87015_img.csv'
    # csv_path = './nyu_test_set/nyu_preprocess_tst/label_csv/00044.csv'
    output_path = './output/vis/pred'
    vis(csv_path, output_path)


if __name__ == '__main__':
    main()
