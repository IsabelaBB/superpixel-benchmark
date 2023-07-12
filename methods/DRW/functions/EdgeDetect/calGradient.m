function [ grad ] = calGradient( img )
%CALGRADIENT 此处显示有关此函数的摘要
%   此处显示详细说明
g_img=rgb2gray(img); 
A=double(g_img);
[GX,GY]=gradient(A);
G=sqrt(GX.*GX+GY.*GY);
grad = G / 255.0;
end

