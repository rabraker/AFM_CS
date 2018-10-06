clear;
clc;
close all;

load('img_data.mat');

I = img_data.cs_im;
pixelifsampled = img_data.pixelifsampled;
maxiter  = round(0.008*256^2);


[ Ir ] = SMP_1D( img_data.cs_im,pixelifsampled,maxiter);












