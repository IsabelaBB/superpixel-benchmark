function wakeup (message)
  printf ("\a%s\n", message);
endfunction

function runSegms(x)
  #for i = 3 : length(list_files) # Avoid the dirs '.' and '..'
    #img_name = list_files(i).name;
  img_name = x.name;
  split_img_name = strsplit(img_name, '.');
  out_name = strcat(split_img_name(1), ".png"){1};
  
  #img = imread(strcat(dir_path, img_name));
  
  #for j = 1 : length(num_superpixels)
    #k = num_superpixels(j);
    #out_path = strcat(segm_path, int2str(k), '/', out_name);
    
    #[labels, numlabels] = slicmex(img,k,compact_factor);
    
    #imwrite(uint16(labels), "test.png");    
  #endfor  
  #endfor
endfunction


dataset_name = 'ECSSD';
num_superpixels = [20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200,...
                      225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475,...
                      500, 750, 1000];
compact_factor = 10;
                      
dir_path = strcat('/home/felipedkstro/Research/datasets/', dataset_name, '/images/');
res_path = strcat('/home/felipedkstro/Research/results/', dataset_name, '/SLIC/'); 
segm_path = strcat(res_path, 'segm/');
                 
list_files = dir(dir_path);

mkdir('/home/felipedkstro/Research/results/', dataset_name);
mkdir res_path;
mkdir segm_path;

for j = 1 : length(num_superpixels)
  mkdir(segm_path,int2str(num_superpixels(j)));
endfor

pararrayfun(-1, @runSegms,list_files);

#[labels, numlabels] = slicmex(img,500,20);%numlabels is the same as number of superpixels


