function [ pb ] = calPbEdge( img )
%CALPBEDGE 此处显示有关此函数的摘要
%   此处显示详细说明

img = im2double(img);
pb = pbCGTG_nonmax(img); %superpixel segmentation-->probability img
[pb, ~] = max(pb, [], 3);

end

