function edgesMap(imgPath,outPath)

addpath(genpath('models'));
addpath(genpath('toolbox'));

%% set opts for training (see edgesTrain.m)
opts=edgesTrain();                % default options (good settings)
opts.modelDir='models/';          % model will be in models/forest
opts.modelFnm='modelBsds';        % model name
opts.nPos=5e5; opts.nNeg=5e5;     % decrease to speedup training
opts.useParfor=0;                 % parallelize if sufficient memory

%% train edge detector (~20m/8Gb per tree, proportional to nPos/nNeg)
model=edgesTrain(opts); % will load model if already trained

%% set detection parameters (can set after training)
model.opts.multiscale=0;          % for top accuracy set multiscale=1
model.opts.sharpen=2;             % for top speed set sharpen=0
model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
model.opts.nThreads=4;            % max number threads for evaluation
model.opts.nms=0;                 % set to true to enable nms

%% iterate over images
img_names = dir(imgPath);

for i=1:length(img_names)

if img_names(i).isdir == 0

image_name = img_names(i).name;
[~, image_stem, ext] = fileparts(image_name);

%% detect edge and visualize results
img = imread([img_names(i).folder, '/', image_stem, ext]);

img_size = length(size(img));
if img_size == 3
    [rows,cols, aa] = size(img);
else
    [rows,cols] = size(img);
    img = repmat(reshape(img, [rows, cols, 1]), [1, 1, 3, 1]);
end

E=edgesDetect(img,model);
imwrite(E, [outPath, '/', image_stem, '.png'])

end

end

