import os
import torchvision.transforms as transforms
import flow_transforms
from imageio import imread, imsave
from skimage import img_as_ubyte
from loss import *


def vis(img_path, csv_path, save_path):
    input_transform = transforms.Compose([
        flow_transforms.ArrayToTensor(),
        transforms.Normalize(mean=[0, 0, 0], std=[255, 255, 255]),
        transforms.Normalize(mean=[0.411, 0.432, 0.45], std=[1, 1, 1])
    ])
    img_ = imread(img_path)
    ori_img = input_transform(img_)
    mean_values = torch.tensor([0.411, 0.432, 0.45], dtype=ori_img.cuda().unsqueeze(0).dtype).view(3, 1, 1)

    spixel_label_map = np.loadtxt(csv_path, delimiter=",")

    n_spixel = len(np.unique(spixel_label_map))
    given_img_np = (ori_img + mean_values).clamp(0, 1).detach().cpu().numpy().transpose(1, 2, 0)
    spixel_bd_image = mark_boundaries(given_img_np / np.max(given_img_np), spixel_label_map.astype(int), color=(1, 1, 1))
    spixel_viz = spixel_bd_image.astype(np.float32).transpose(2, 0, 1)

    imgId = os.path.basename(img_path)[:-4]

    # ************************ Save all result*******************************************
    # save spixel viz
    if not os.path.isdir(os.path.join(save_path, 'spixel_viz')):
        os.makedirs(os.path.join(save_path, 'spixel_viz'))
    spixl_save_name = os.path.join(save_path, 'spixel_viz', imgId + '_sPixel.png')
    imsave(spixl_save_name, img_as_ubyte(spixel_viz.transpose(1, 2, 0)))

    return


def main():
    img_path = './BSD500/ori_sz/test/100099_img.jpg'
    # img_path = './nyu_test_set/nyu_preprocess_tst/img/00044.jpg'
    csv_path = '/media/yuanqing/ssd/code/visual_tracking/SNIC_mex/output/BSD/600/100099_img.csv'
    output_path = './output/vis/bsd/snic'
    # if not os.path.isdir(output_path):
    #     os.makedirs(output_path)
    vis(img_path, csv_path, output_path)


if __name__ == '__main__':
    main()
