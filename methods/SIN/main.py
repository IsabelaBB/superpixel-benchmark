import argparse
import os
import shutil
import time

import torch.nn.parallel
import torch.backends.cudnn as cudnn
import torch.optim
import torch.utils.data
import torchvision.transforms as transforms
from torchvision.utils import make_grid

import flow_transforms
import models
import datasets
from loss import compute_labxy_loss
import datetime
from tensorboardX import SummaryWriter
from train_util import *
from models.model_util import update_spixel_map

'''
Main code for training 

author: Fengting 
Last modification: March 8th, 2019
'''


model_names = sorted(name for name in models.__dict__
                     if not name.startswith("__"))
dataset_names = sorted(name for name in datasets.__all__)


parser = argparse.ArgumentParser(description='PyTorch SpixelFCN Training on BSDS500',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# ================ training setting ====================
parser.add_argument('--dataset', metavar='DATASET', default='VOC',  choices=dataset_names,
                    help='dataset type : ' +  ' | '.join(dataset_names))
parser.add_argument('--arch', '-a', metavar='ARCH', default='SpixelNet1l_bn',  help='model architecture')
parser.add_argument('--data', metavar='DIR',default='./BSD500/pre_pro', help='path to input dataset')
parser.add_argument('--savepath',default='./snapshot', help='path to save ckpt')


parser.add_argument('--train_img_height', '-t_imgH', default=193,  type=int, help='img height')
parser.add_argument('--train_img_width', '-t_imgW', default=193, type=int, help='img width')
parser.add_argument('--input_img_height', '-v_imgH', default=319, type=int, help='img height_must be 16*n')  #
parser.add_argument('--input_img_width', '-v_imgW', default=319,  type=int, help='img width must be 16*n')

# ======== learning schedule ================
parser.add_argument('-j', '--workers', default=4, type=int, metavar='N', help='number of data loading workers')
parser.add_argument('--epochs', default=3000000, type=int, metavar='N', help='number of total epoches, make it big enough to follow the iteration maxmium')
parser.add_argument('--start_epoch', default=0, type=int, metavar='N',  help='manual epoch number (useful on restarts)')
parser.add_argument('--epoch_size', default= 6000,  help='choose any value > 408 to use all the train and val data')
parser.add_argument('-b', '--batch-size', default=4, type=int,   metavar='N', help='mini-batch size')

parser.add_argument('--solver', default='adam',choices=['adam', 'sgd'], help='solver algorithms, we use adam')
parser.add_argument('--lr', '--learning-rate', default=5e-05, type=float,metavar='LR', help='initial learning rate')
parser.add_argument('--momentum', default=0.9, type=float, metavar='M',  help='momentum for sgd, alpha parameter for adam')
parser.add_argument('--beta', default=0.999, type=float, metavar='M',   help='beta parameter for adam')
parser.add_argument('--weight_decay', '--wd', default=4e-4, type=float, metavar='W', help='weight decay')
parser.add_argument('--bias_decay', default=0, type=float, metavar='B', help='bias decay, we never use it')
parser.add_argument('--milestones', default=[200000], metavar='N', nargs='*', help='epochs at which learning rate is divided by 2')
parser.add_argument('--additional_step', default= 100000, help='the additional iteration, after lr decay')

# ================= other setting ===================
parser.add_argument('--gpu', default= '0', type=str, help='gpu id')
parser.add_argument('--print_freq', '-p', default=10, type=int,  help='print frequency (step)')
parser.add_argument('--record_freq', '-rf', default=5, type=int,  help='record frequency (epoch)')
parser.add_argument('--label_factor', default=5, type=int, help='constant multiplied to label index for viz.')
parser.add_argument('--pretrained', dest='pretrained', default='./snapshot/BSD500/full_model_12_23/300000_step.tar', help='path to pre-trained model')
parser.add_argument('--no-date', action='store_true',  help='don\'t append date timestamp to folder' )


best_EPE = -1
n_iter = 0
args = parser.parse_args()

# !----- NOTE the current code does not support cpu training -----!
os.environ['CUDA_VISIBLE_DEVICES'] = args.gpu
if torch.cuda.is_available():
    device = torch.device("cuda")
else:
    print('Current code does not support CPU training! Sorry about that.')
    exit(1)


def main():
    global args, best_EPE, save_path, intrinsic

    # ============= savor setting ===================
    save_path = '{}_{}_{}epochs{}_b{}_lr{}'.format(
        args.arch,
        args.solver,
        args.epochs,
        '_epochSize'+str(args.epoch_size) if args.epoch_size > 0 else '',
        args.batch_size,
        args.lr,
    )
    if not args.no_date:
        timestamp = datetime.datetime.now().strftime("%y_%m_%d_%H_%M")
    else:
        timestamp = ''
    save_path = os.path.abspath(args.savepath) + '/' + os.path.join(args.dataset, save_path  +  '_' + timestamp )

    # ==========  Data loading code ==============
    input_transform = transforms.Compose([
        flow_transforms.ArrayToTensor(),
        transforms.Normalize(mean=[0,0,0], std=[255,255,255]),
        transforms.Normalize(mean=[0.411,0.432,0.45], std=[1,1,1])
    ])

    val_input_transform = transforms.Compose([
        flow_transforms.ArrayToTensor(),
        transforms.Normalize(mean=[0, 0, 0], std=[255, 255, 255]),
        transforms.Normalize(mean=[0.411, 0.432, 0.45], std=[1, 1, 1])
    ])

    target_transform = transforms.Compose([
        flow_transforms.ArrayToTensor(),
    ])

    co_transform = flow_transforms.Compose([
            flow_transforms.RandomCrop((args.train_img_height, args.train_img_width)),
            flow_transforms.RandomVerticalFlip(),
            flow_transforms.RandomHorizontalFlip()
        ])

    print("=> loading img pairs from '{}'".format(args.data))
    train_set, val_set = datasets.__dict__[args.dataset](
        args.data,
        transform=input_transform,
        val_transform=val_input_transform,
        target_transform=target_transform,
        co_transform=co_transform
    )
    print('{} samples found, {} train samples and {} val samples '.format(len(val_set)+len(train_set),
                                                                           len(train_set),
                                                                           len(val_set)))
    train_loader = torch.utils.data.DataLoader(
        train_set, batch_size=args.batch_size,
        num_workers=args.workers, pin_memory=True, shuffle=True, drop_last=True)

    val_loader = torch.utils.data.DataLoader(
        val_set, batch_size=args.batch_size,
        num_workers=args.workers, pin_memory=True, shuffle=False, drop_last=True)

    # ============== create model ====================
    if args.pretrained:
        network_data = torch.load(args.pretrained)
        args.arch = network_data['arch']
        print("=> using pre-trained model '{}'".format(args.arch))
    else:
        network_data = None
        print("=> creating model '{}'".format(args.arch))

    model = models.__dict__[args.arch](data=network_data).cuda()
    model = torch.nn.DataParallel(model).cuda()
    cudnn.benchmark = True

    # =========== creat optimizer, we use adam by default ==================
    assert(args.solver in ['adam', 'sgd'])
    print('=> setting {} solver'.format(args.solver))
    param_groups = [{'params': model.module.bias_parameters(), 'weight_decay': args.bias_decay},
                    {'params': model.module.weight_parameters(), 'weight_decay': args.weight_decay}]
    if args.solver == 'adam':
        optimizer = torch.optim.Adam(param_groups, args.lr,
                                     betas=(args.momentum, args.beta))
    elif args.solver == 'sgd':
        optimizer = torch.optim.SGD(param_groups, args.lr,
                                    momentum=args.momentum)

    # for continues training
    if args.pretrained and ('dataset' in network_data):
        if args.pretrained and args.dataset == network_data['dataset'] :
            optimizer.load_state_dict(network_data['optimizer'])
            best_EPE = network_data['best_EPE']
            args.start_epoch = network_data['epoch']
            save_path = os.path.dirname(args.pretrained)

    print('=> will save everything to {}'.format(save_path))
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    train_writer = SummaryWriter(os.path.join(save_path, 'train'))
    val_writer = SummaryWriter(os.path.join(save_path, 'val'))

    for epoch in range(args.start_epoch, args.epochs):
        # train for one epoch
        train_avg_loss, iteration = train(train_loader, model, optimizer, epoch,
                                                         train_writer)
        if epoch % args.record_freq == 0:
            train_writer.add_scalar('Mean avg_loss', train_avg_loss, epoch)

        # evaluate on validation set and save the module( and choose the best)
        with torch.no_grad():
            avg_loss = validate(val_loader, model, epoch, val_writer)
            if epoch % args.record_freq == 0:
                val_writer.add_scalar('Mean avg_loss', avg_loss, epoch)

        rec_dict = {
                'epoch': epoch + 1,
                'arch': args.arch,
                'state_dict': model.module.state_dict(),
                'best_EPE': best_EPE,
                'optimizer': optimizer.state_dict(),
                'dataset': args.dataset
            }

        if (iteration) >= (args.milestones[-1] + args.additional_step):
            save_checkpoint(rec_dict, is_best =False, filename='%d_step.tar' % iteration)
            print("Train finished!")
            break

        if best_EPE < 0:
            best_EPE = avg_loss
        is_best = avg_loss < best_EPE
        best_EPE = min(avg_loss, best_EPE)
        save_checkpoint(rec_dict, is_best)


def train(train_loader, model, optimizer, epoch, train_writer):
    global n_iter, args, intrinsic
    batch_time = AverageMeter()
    data_time = AverageMeter()

    total_loss = AverageMeter()

    epoch_size = len(train_loader) if args.epoch_size == 0 else min(len(train_loader), args.epoch_size)

    # switch to train mode
    model.train()
    end = time.time()
    iteration = 0

    for i, (input, label) in enumerate(train_loader):

        iteration = i + epoch * epoch_size

        # ========== adjust lr if necessary  ===============
        if (iteration + 1) in args.milestones:
            state_dict = optimizer.state_dict()
            for param_group in state_dict['param_groups']:
                param_group['lr'] = args.lr * ((0.5) ** (args.milestones.index(iteration + 1) + 1))
            optimizer.load_state_dict(state_dict)

        # ========== complete data loading ================
        input_gpu = input.to(device)
        label_gpu = label.to(device)
        torch.cuda.synchronize()
        data_time.update(time.time() - end)


        # ========== predict association map ============
        prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h = model(input_gpu)
        fast_loss = compute_labxy_loss(prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h,
                                                           label_gpu)

        # ========= back propagate ===============
        optimizer.zero_grad()
        fast_loss.backward()
        optimizer.step()

        # ========  measure batch time ===========
        torch.cuda.synchronize()
        batch_time.update(time.time() - end)
        end = time.time()

        # =========== record and display the loss ===========
        # record loss and EPE
        total_loss.update(fast_loss.item(), input_gpu.size(0))

        if i % args.print_freq == 0:
            print('train Epoch: [{0}][{1}/{2}]\t Time {3}\t Data {4}\t Total_loss {5}\t'
                  .format(epoch, i, epoch_size, batch_time, data_time, total_loss))

            train_writer.add_scalar('Train_loss', fast_loss.item(), i + epoch*epoch_size)
            train_writer.add_scalar('learning rate', optimizer.param_groups[0]['lr'], i + epoch * epoch_size)

        n_iter += 1
        if i >= epoch_size:
            break

        if (iteration) >= (args.milestones[-1] + args.additional_step):
            break

    # =========== write information to tensorboard ===========
    if epoch % args.record_freq == 0:
        train_writer.add_scalar('Train_loss_epoch', fast_loss.item(),  epoch )

        #save image
        mean_values = torch.tensor([0.411, 0.432, 0.45], dtype=input_gpu.dtype).view(3, 1, 1)
        input_l_save = (make_grid((input + mean_values).clamp(0, 1), nrow=args.batch_size))
        label_save = make_grid(args.label_factor * label)

        train_writer.add_image('Input', input_l_save, epoch)
        train_writer.add_image('label', label_save, epoch)

        curr_spixl_map = update_spixel_map(input_gpu, prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h)
        spixel_lab_save = make_grid(curr_spixl_map, nrow=args.batch_size)[0, :, :]
        spixel_viz, _ = get_spixel_image(input_l_save, spixel_lab_save)
        train_writer.add_image('Spixel viz', spixel_viz, epoch)

        if not isinstance(curr_spixl_map, np.ndarray):
            spix_index_np = curr_spixl_map[0, 0, :, :].detach().cpu().numpy().transpose(0, 1)
        else:
            spix_index_np = curr_spixl_map[0, 0, :, :]
        # save the unique maps as csv, uncomment it if needed
        if not os.path.isdir(os.path.join(save_path, 'map_csv')):
            os.makedirs(os.path.join(save_path, 'map_csv'))
        output_path = os.path.join(save_path, 'map_csv', str(i) + '.csv')
        np.savetxt(output_path, spix_index_np.astype(int), fmt='%i', delimiter=",")

        #save associ map,  --- for debug only
        # _, prob_idx = torch.max(output, dim=1, keepdim=True)
        # prob_map_save = make_grid(assign2uint8(prob_idx))
        # train_writer.add_image('assigment idx', prob_map_save, epoch)

        print('==> write train step %dth to tensorboard' % i)


    return total_loss.avg, iteration


def validate(val_loader, model, epoch, val_writer):
    global n_iter,   args,    intrinsic
    batch_time = AverageMeter()
    data_time = AverageMeter()

    total_loss = AverageMeter()

    # set the validation epoch-size, we only randomly val. 400 batches during training to save time
    epoch_size = min(len(val_loader), 400)

    # switch to train mode
    model.eval()
    end = time.time()

    for i, (input, label) in enumerate(val_loader):
        input_gpu = input.to(device)
        label_gpu = label.to(device)
        torch.cuda.synchronize()
        data_time.update(time.time() - end)

        # compute output
        with torch.no_grad():
            prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h = model(input_gpu)
            fast_loss = compute_labxy_loss(prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h,
                                                               label_gpu)

        # measure elapsed time
        torch.cuda.synchronize()
        batch_time.update(time.time() - end)
        end = time.time()

        # record loss and EPE
        total_loss.update(fast_loss.item(), input_gpu.size(0))

        if i % args.print_freq == 0:
            print('val Epoch: [{0}][{1}/{2}]\t Time {3}\t Data {4}\t Total_loss {5}\t'
                  .format(epoch, i, epoch_size, batch_time, data_time, total_loss))

        if i >= epoch_size:
            break

    # =============  write result to tensorboard ======================
    if epoch % args.record_freq == 0:
        val_writer.add_scalar('Train_loss_epoch', fast_loss.item(), epoch)

        mean_values = torch.tensor([0.411, 0.432, 0.45], dtype=input_gpu.dtype).view(3, 1, 1)
        input_l_save = (make_grid((input + mean_values).clamp(0, 1), nrow=args.batch_size))
        label_save = make_grid(args.label_factor * label)


        curr_spixl_map = update_spixel_map(input_gpu, prob0_v, prob0_h, prob1_v, prob1_h, prob2_v, prob2_h, prob3_v, prob3_h)
        spixel_lab_save = make_grid(curr_spixl_map, nrow=args.batch_size)[0, :, :]
        spixel_viz, _ = get_spixel_image(input_l_save, spixel_lab_save)


        val_writer.add_image('Input', input_l_save, epoch)
        val_writer.add_image('label', label_save, epoch)
        val_writer.add_image('Spixel viz', spixel_viz, epoch)

        # --- for debug
        #     _, prob_idx = torch.max(assign, dim=1, keepdim=True)
        #     prob_map_save = make_grid(assign2uint8(prob_idx))
        #     val_writer.add_image('assigment idx level %d' % j, prob_map_save, epoch)

        print('==> write val step %dth to tensorboard' % i)

    return total_loss.avg


def save_checkpoint(state, is_best, filename='checkpoint.tar'):
    torch.save(state, os.path.join(save_path,filename))
    if is_best:
        shutil.copyfile(os.path.join(save_path,filename), os.path.join(save_path,'model_best.tar'))


if __name__ == '__main__':
    main()
