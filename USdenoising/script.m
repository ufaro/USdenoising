clear all
close all
clc


%% Script performs Data-driven Kronecker-product wavelet-Fisz Despeckling 




x = 'DATA'; % .mat file first dimension (time), second and third (space)


smin=0; smax=255;

x=(x-min(x(:)))*smax/(max(x(:))-min(x(:)) )+smin;


sigma = 3; h = @(x) sigma^2*x ; % sigma is an estimate of the standard deviation.

Reconstructed  = ddhf_est_Iso_Kronecker(x, h, 1,1);
