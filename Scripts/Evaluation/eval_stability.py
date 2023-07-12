
#import argparse
import pandas as pd
import os.path

#parser = argparse.ArgumentParser()
#parser.add_argument('--input',help='The path of the source image')
#parser.add_argument('--output', help='The path of the labeled image')

#DATASETS = ['Birds', 'Sky', 'ECSSD', 'Insects']
METHODS = ['CRS', 'DAL-HERS', 'DISF', 'DRW', 'ERGC', 'ERS', 'ETPS', 'GMMSP', 'GRID', 'IBIS', 'ISF', 'LNSNet', 'LSC', 'ODISF', 'RSS', 'SCALP', 'SEEDS', 'SH', 'SLIC','SNIC', 'SICLE']
METRICS = ['BR', 'UE', 'EV', 'SIRS']
#SPXS = [25,50,75,100,200,300,400,500,600,700,800,900,1000]

DATASETS = ['test10']
SPXS = [50,100]

OUTPUT_PATH="../../RESULTS/Eval"
INUT_PATH="../../RESULTS/Eval"

data = {}

for dataset in DATASETS:
    for method in METHODS:
        for metric in METRICS:

          # get min / max / dp / avg
          data_spx, data_min, data_max, data_std = [], [], [], []

          for num_superpixel in SPXS:
            input_file = INUT_PATH + '/' +method+'/'+dataset+'/'+metric+'/'+method+'-'+num_superpixel+'.txt'

            data_imgs = pd.read_csv(input_file, sep=" ", names=["Image", "Superpixels", metric])
            
            data_spx.append(data_imgs["Superpixels"].mean())
            data_min.append(data_imgs[metric].min())
            data_std.append(data_imgs[metric].std())
            data_max.append(data_imgs[metric].max())

          key=metric+'-'+dataset+'-'+method+'-'
          data[key+'MIN'] = pd.DataFrame(data={'Superpixels':data_spx, metric:data_min})
          data[key+'MAX'] = pd.DataFrame(data={'Superpixels':data_spx, metric:data_max})
          data[key+'STD'] = pd.DataFrame(data={'Superpixels':data_spx, metric:data_std})

          for type in ['MIN', 'MAX', 'STD']:
            
            output_file = OUTPUT_PATH + '/' +method+'/'+dataset+'/stability/'+ type
            if not os.path.exists(output_file):
              os.makedirs(output_file)
            output_file = output_file+'/'+method+'-'+dataset+'-'+metric+'.txt'
            
            key=metric+'-'+dataset+'-'+method+'-'+type
            data[key].to_csv(output_file, header=[metric, type], index=None, sep=' ', mode='w')

          